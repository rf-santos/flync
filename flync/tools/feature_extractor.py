#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
from typing import Dict, List, Optional
from dataclasses import dataclass
import logging
from scripts.parallel_executor import ParallelExecutor, ExecutionConfig

@dataclass
class FeatureConfig:
    """Configuration for feature extraction"""
    paths_table: Path
    bed_path: Path
    outpath: Path
    prefix: str

class FeatureExtractor:
    """Handles extraction and processing of genomic features"""
    
    def __init__(self, config: FeatureConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.executor = ParallelExecutor(ExecutionConfig(), logger)
        
        # Column definitions
        self.bw_cols = ["name", "size", "covered_bases", "sum", "mean0", "mean", "min", "max"]
        self.bb_cols = ["name", "covered_percent", "mean", "min", "max"]
        
    def load_data(self) -> tuple:
        """Load and validate input data"""
        try:
            readin = pd.read_csv(self.config.paths_table, delimiter=',', header=None)
            readin.rename(columns={0: 'track', 1:'url'}, inplace=True)
            
            bed = pd.read_csv(self.config.bed_path, delimiter='\t', header=None)
            bed.columns = ["chr", "start", "end", "name", "score", "strand"]
            
            return readin, bed
            
        except Exception as e:
            self.logger.error(f"Error loading input data: {str(e)}")
            raise

    def _process_file(self, file_info: tuple) -> pd.DataFrame:
        """Process a single feature file"""
        track_name, file_path = file_info
        try:
            df = pd.read_csv(file_path, delimiter='\t', header=None)
            
            if track_name in ['GCcont', 'phastCons27', 'phyloP27', 'phyloP124', 'Pol2_S2', 'H3K4me3_S2']:
                df.columns = self.bw_cols
            else:
                df.columns = self.bb_cols
                
            df['name'] = df['name'].map(lambda x: x.rstrip('.'))
            return track_name, df
            
        except Exception as e:
            self.logger.error(f"Error processing {track_name}: {str(e)}")
            return track_name, pd.DataFrame()

    def process_feature_files(self, readin: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """Process all feature files in parallel"""
        file_info = list(zip(readin['track'], readin['url']))
        results = self.executor.parallel_io(self._process_file, file_info)
        
        return {name: df for name, df in results if not df.empty}

    def calculate_tss_features(self, features: pd.DataFrame, dfs: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Calculate TSS-related features"""
        try:
            # Process TSS peaks
            features = pd.merge(features, dfs['CAGE_pos']['name', 'startPosTSS'], on=['name'], how='outer')
            features = pd.merge(features, dfs['CAGE_neg']['name', 'endNegTSS'], on=['name'], how='outer')
            
            tss_peak = features[["startPosTSS", "endNegTSS"]]
            tss_peak = abs(tss_peak)
            tss_peak = tss_peak.max(axis=1)
            features["bestTSS"] = tss_peak
            features.drop(columns=["startPosTSS", "endNegTSS"], inplace=True)

            # Process whole transcript TSS features
            features = pd.merge(features, dfs['CAGE_pos_whole_trans'][['name', 'max']], on=['name'], how='outer')
            features = pd.merge(features, dfs['CAGE_neg_whole_trans'][['name', 'min']], on=['name'], how='outer')
            features.rename(columns={'max': 'PosTSS_inside', 'min': 'NegTSS_inside'}, inplace=True)

            tss_peak_ins = features[["PosTSS_inside", "NegTSS_inside"]]
            tss_peak_ins = abs(tss_peak_ins)
            tss_peak_ins = tss_peak_ins.max(axis=1)
            features["bestTSS_inside"] = tss_peak_ins
            features.drop(columns=["PosTSS_inside", "NegTSS_inside"], inplace=True)

            return features
            
        except Exception as e:
            self.logger.error(f"Error calculating TSS features: {str(e)}")
            raise

    def calculate_coverage_features(self, features: pd.DataFrame, dfs: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Calculate coverage-related features"""
        try:
            # Add metrics from bigWig/bigBed files
            features = pd.merge(features, dfs['GCcont'][['name', 'mean']], on=['name'], how='outer')
            features.rename(columns={'mean': 'mean_gc'}, inplace=True)

            features = pd.merge(features, dfs['ReMap'][['name', 'mean']], on=['name'], how='outer')
            features.rename(columns={'mean': 'mean_remap'}, inplace=True)

            # Calculate coverage features
            for track, new_col in [('H3K4me3_S2', 'cov_me3'), ('JASPAR_TF', 'cov_tfbs'), ('Pol2_S2', 'cov_pol2')]:
                if track in dfs:
                    features = pd.merge(features, dfs[track][['name', 'covered_bases']], on=['name'], how='outer')
                    if new_col in ['cov_me3', 'cov_pol2']:
                        features[new_col] = features['covered_bases'] / features['length']
                    else:
                        features.rename(columns={'covered_bases': new_col}, inplace=True)
                    features.drop(columns=['covered_bases'], inplace=True, errors='ignore')

            return features
            
        except Exception as e:
            self.logger.error(f"Error calculating coverage features: {str(e)}")
            raise

    def extract_features(self) -> pd.DataFrame:
        """Main method to extract all features"""
        self.logger.info("Starting feature extraction")
        
        try:
            # Load input data
            readin, bed = self.load_data()
            
            # Process bed file info
            features = bed[["name"]]
            features = pd.merge(features, bed[['name', 'start', 'end']], on=['name'])
            features["length"] = features["end"] - features["start"]
            
            # Process all feature files
            dfs = self.process_feature_files(readin)
            
            # Calculate features
            features = self.calculate_tss_features(features, dfs)
            features = self.calculate_coverage_features(features, dfs)
            
            # Handle missing values
            features = features.fillna(0)
            
            # Save results
            outfile = self.config.outpath / f"{self.config.prefix}.csv"
            features.to_csv(outfile, index=False)
            
            self.logger.info(f"Feature extraction complete. Results saved to {outfile}")
            return features
            
        except Exception as e:
            self.logger.error(f"Feature extraction failed: {str(e)}")
            raise

def main(appdir: Path, workdir: Path, bed_path: Optional[Path] = None):
    """CLI entry point"""
    try:
        # Setup paths
        paths_table = workdir / 'results/non-coding/features/paths.csv'
        outpath = workdir / 'results'
        
        if bed_path is None:
            bed_path = workdir / 'results/new-non-coding.bed'
            
        prefix = bed_path.stem
        
        # Configure and run feature extraction
        config = FeatureConfig(
            paths_table=paths_table,
            bed_path=bed_path,
            outpath=outpath,
            prefix=prefix
        )
        
        logger = logging.getLogger('flync.features')
        extractor = FeatureExtractor(config, logger)
        extractor.extract_features()
        
    except Exception as e:
        logger.error(f"Feature extraction failed: {str(e)}")
        raise

if __name__ == '__main__':
    import sys
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]) if len(sys.argv) > 3 else None)