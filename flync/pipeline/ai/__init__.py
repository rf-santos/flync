#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
import pickle
from typing import List, Dict, Optional
import logging
from dataclasses import dataclass
from scripts.progress_manager import ProgressManager, TaskProgress
import joblib

@dataclass
class AIPipelineConfig:
    """Configuration for AI pipeline stage"""
    workdir: Path
    model_path: Path
    threads: int

class MLPipeline:
    """Handles the AI/ML stage of the pipeline"""
    
    def __init__(self, config: AIPipelineConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.progress = ProgressManager(logger)
        
    def load_model(self) -> object:
        """Load the pre-trained ML model"""
        try:
            model = joblib.load(self.config.model_path)
            self.logger.info(f"Model loaded successfully from {self.config.model_path}")
            return model
        except Exception as e:
            self.progress.error("Failed to load ML model", e)
            raise

    def prepare_features(self, feature_matrix: pd.DataFrame) -> pd.DataFrame:
        """Prepare feature matrix for prediction"""
        try:
            # Required features based on the trained model
            required_features = [
                'length', 'bestTSS', 'bestTSS_inside', 'mean_gc',
                'mean_remap', 'cov_tfbs', 'cov_pol2', 'mean_pcons27',
                'mean_pPcons124'
            ]
            
            # Validate features
            missing_features = [f for f in required_features if f not in feature_matrix.columns]
            if missing_features:
                raise ValueError(f"Missing required features: {missing_features}")
                
            # Select and order features
            X = feature_matrix[required_features]
            
            # Fill missing values
            X = X.fillna(0)
            
            return X
            
        except Exception as e:
            self.progress.error("Failed to prepare features", e)
            raise

    def predict(self, features: pd.DataFrame) -> pd.DataFrame:
        """Run ML prediction on prepared features"""
        with self.progress.task_progress(TaskProgress(name="Running ML prediction")):
            try:
                # Load model
                model = self.load_model()
                
                # Prepare features
                X = self.prepare_features(features)
                
                # Get predictions and probabilities
                predictions = model.predict(X)
                probabilities = model.predict_proba(X)
                
                # Create results DataFrame
                results = pd.DataFrame({
                    'name': features['name'],
                    'lncRNA': predictions,
                    'Prob_False': probabilities[:, 0],
                    'Prob_True': probabilities[:, 1]
                })
                
                self.progress.status(f"Predictions complete: {sum(predictions)} lncRNAs identified")
                return results
                
            except Exception as e:
                self.progress.error("Prediction failed", e)
                raise

    def analyze_predictions(self, predictions: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
        """Analyze predictions and compile final results"""
        try:
            # Merge predictions with features
            results = pd.merge(predictions, features, on='name')
            
            # Calculate additional metrics if available
            if 'expression' in features.columns:
                results['log2_expression'] = np.log2(features['expression'] + 1)
                
            # Add genome browser links
            results['ucsc_link'] = results.apply(
                lambda x: f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=dm6&position={x['chr']}:{x['start']}-{x['end']}",
                axis=1
            )
            
            return results
            
        except Exception as e:
            self.progress.error("Failed to analyze predictions", e)
            raise

    def save_results(self, results: pd.DataFrame) -> None:
        """Save prediction results"""
        try:
            # Save complete results
            results_path = self.config.workdir / 'results' / 'lncrna_predictions.csv'
            results.to_csv(results_path, index=False)
            
            # Save positive predictions only
            pos_results = results[results['lncRNA']]
            pos_path = self.config.workdir / 'results' / 'lncrna_predictions_positive.csv'
            pos_results.to_csv(pos_path, index=False)
            
            # Save prediction summary
            summary = {
                'total_transcripts': len(results),
                'predicted_lncrnas': len(pos_results),
                'prediction_rate': len(pos_results) / len(results) * 100
            }
            
            summary_path = self.config.workdir / 'results' / 'prediction_summary.txt'
            with open(summary_path, 'w') as f:
                for key, value in summary.items():
                    f.write(f"{key}: {value}\n")
            
            self.logger.info(f"Results saved to {results_path}")
            self.logger.info(f"Positive predictions saved to {pos_path}")
            self.logger.info(f"Summary saved to {summary_path}")
            
        except Exception as e:
            self.progress.error("Failed to save results", e)
            raise

    def run_pipeline(self, features_path: Optional[Path] = None) -> None:
        """Run the complete AI pipeline"""
        try:
            # Load features
            if features_path is None:
                features_path = self.config.workdir / 'results' / 'features.csv'
                
            features = pd.read_csv(features_path)
            
            # Run prediction pipeline
            predictions = self.predict(features)
            results = self.analyze_predictions(predictions, features)
            self.save_results(results)
            
            # Print summary
            self.progress.print_summary()
            
        except Exception as e:
            self.progress.error("AI pipeline execution failed", e)
            raise