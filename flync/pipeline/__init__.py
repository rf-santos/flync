#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path
from enum import Enum, auto
from typing import Optional, List, Dict, Any
import logging
import json
from datetime import datetime
from dataclasses import asdict

from scripts.logger import setup_logger
from scripts.progress_manager import ProgressManager, TaskProgress, TaskStatus
from scripts.pipeline_bio import BioinformaticsPipeline, BioinformaticsConfig 
from scripts.pipeline_ai import MLPipeline, AIPipelineConfig
from scripts.config_manager import ConfigManager, PipelineConfig
from scripts.anime import ProgressAnimation, AnimationType

class PipelineStage(Enum):
    """Pipeline execution stages"""
    INIT = auto()
    PREPARE_GENOME = auto()
    MAP_READS = auto()
    ASSEMBLE = auto()
    EXTRACT_FEATURES = auto()
    PREDICT = auto()
    COMPLETE = auto()
    ERROR = auto()

class Pipeline:
    """Main pipeline controller"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.workdir = config.output
        self.logger = setup_logger(self.workdir)
        self.progress = ProgressManager(self.logger)
        self._stage = PipelineStage.INIT
        self._checkpoint_file = self.workdir / "pipeline_state.json"
        
    def _save_checkpoint(self) -> None:
        """Save pipeline state to checkpoint file"""
        checkpoint = {
            'stage': self._stage.name,
            'timestamp': datetime.now().isoformat(),
            'config': asdict(self.config)
        }
        
        try:
            with open(self._checkpoint_file, 'w') as f:
                json.dump(checkpoint, f, indent=2)
        except Exception as e:
            self.logger.warning(f"Failed to save checkpoint: {str(e)}")

    def _load_checkpoint(self) -> Optional[Dict[str, Any]]:
        """Load pipeline state from checkpoint file"""
        try:
            if self._checkpoint_file.exists():
                with open(self._checkpoint_file) as f:
                    return json.load(f)
        except Exception as e:
            self.logger.warning(f"Failed to load checkpoint: {str(e)}")
        return None

    def _setup_environment(self) -> None:
        """Setup conda environments and dependencies"""
        try:
            # Create necessary directories
            self.workdir.mkdir(parents=True, exist_ok=True)
            
            # Validate conda environments
            env_names = ['featureMod', 'mapMod', 'assembleMod', 'codMod', 'predictMod']
            for env in env_names:
                result = subprocess.run(['conda', 'env', 'list'], capture_output=True, text=True)
                if env not in result.stdout:
                    self.progress.error(f"Required conda environment '{env}' not found")
                    raise RuntimeError(f"Missing conda environment: {env}")
                    
            self.logger.info("Environment setup completed successfully")
            
        except Exception as e:
            self.progress.error("Environment setup failed", e)
            raise

    def _activate_conda_env(self, env_name: str) -> None:
        """Activate a conda environment"""
        try:
            subprocess.run(['conda', 'activate', env_name], check=True)
            self.logger.debug(f"Activated conda environment: {env_name}")
        except subprocess.CalledProcessError as e:
            self.progress.error(f"Failed to activate conda environment: {env_name}", e)
            raise

    def _run_bioinformatics_stage(self) -> None:
        """Run the bioinformatics pipeline stage"""
        bio_config = BioinformaticsConfig(
            workdir=self.workdir,
            threads=self.config.threads
        )
        
        with self.progress.task_progress(TaskProgress(
            name="Bioinformatics Pipeline",
            desc="Running bioinformatics analysis"
        )):
            pipeline = BioinformaticsPipeline(bio_config, self.logger)
            
            # Run stages based on checkpoint
            checkpoint = self._load_checkpoint()
            if checkpoint and checkpoint['stage'] != PipelineStage.INIT.name:
                self.logger.info(f"Resuming from checkpoint: {checkpoint['stage']}")
                self._stage = PipelineStage[checkpoint['stage']]
            
            try:
                if self._stage == PipelineStage.INIT:
                    pipeline.prepare_genome()
                    self._stage = PipelineStage.PREPARE_GENOME
                    self._save_checkpoint()
                
                if self._stage == PipelineStage.PREPARE_GENOME:
                    pipeline.process_reads(self.config)
                    self._stage = PipelineStage.MAP_READS
                    self._save_checkpoint()
                    
                if self._stage == PipelineStage.MAP_READS:
                    pipeline.assemble_transcriptome()
                    self._stage = PipelineStage.ASSEMBLE
                    self._save_checkpoint()
                    
                if self._stage == PipelineStage.ASSEMBLE:
                    pipeline.extract_features()
                    self._stage = PipelineStage.EXTRACT_FEATURES
                    self._save_checkpoint()
                    
            except Exception as e:
                self._stage = PipelineStage.ERROR
                self._save_checkpoint()
                raise

    def _run_ml_stage(self) -> None:
        """Run the machine learning pipeline stage"""
        model_path = Path(__file__).parent.parent / 'model' / 'rf_dm6_lncrna_classifier.model'
        
        ai_config = AIPipelineConfig(
            workdir=self.workdir,
            model_path=model_path,
            threads=self.config.threads
        )
        
        with self.progress.task_progress(TaskProgress(
            name="ML Pipeline",
            desc="Running ML prediction"
        )):
            pipeline = MLPipeline(ai_config, self.logger)
            
            try:
                pipeline.run_pipeline()
                self._stage = PipelineStage.COMPLETE
                self._save_checkpoint()
            except Exception as e:
                self._stage = PipelineStage.ERROR
                self._save_checkpoint()
                raise

    def run(self) -> None:
        """Run the complete pipeline"""
        try:
            with ProgressAnimation(AnimationType.BRAILLE):
                # Setup environment
                self._setup_environment()
                
                # Run pipeline stages
                self._run_bioinformatics_stage()
                self._run_ml_stage()
                
                self.progress.print_summary()
                
        except KeyboardInterrupt:
            self._stage = PipelineStage.ERROR
            self._save_checkpoint()
            self.progress.error("Pipeline interrupted by user")
            sys.exit(1)
        except Exception as e:
            self._stage = PipelineStage.ERROR
            self._save_checkpoint()
            self.progress.error(f"Pipeline failed at stage {self._stage.name}", e)
            raise

def main():
    """CLI entry point"""
    # Parse command line arguments
    import argparse
    parser = argparse.ArgumentParser(description="FLYNC Pipeline Controller")
    parser.add_argument('--config', '-c', type=Path, required=True,
                      help='Path to configuration file')
    args = parser.parse_args()
    
    try:
        # Load and validate configuration
        logger = logging.getLogger("flync")
        config_manager = ConfigManager(logger)
        config = config_manager.load_config(args.config)
        
        # Run pipeline
        pipeline = Pipeline(config)
        pipeline.run()
        
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()