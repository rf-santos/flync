#!/usr/bin/env python3

import logging
from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass
import pandas as pd
from .read_mapper import ReadMapper, ReadMappingConfig
from .transcript_assembler import TranscriptomeAssembler, AssemblyConfig
from .parallel_executor import ParallelExecutor, ExecutionConfig

@dataclass
class BioinformaticsConfig:
    """Configuration for bioinformatics pipeline stage"""
    workdir: Path
    threads: int
    genome_index: Path
    reference_gtf: Path
    reference_genome: Path
    splice_sites: Path
    min_coverage: float = 2.5
    min_fpkm: float = 0.1
    min_tpm: float = 0.1
    min_transcript_length: int = 200
    compressed: bool = True
    keep_intermediates: bool = False

class BioinformaticsPipeline:
    """Manages the bioinformatics stage of the pipeline"""
    
    def __init__(self, config: BioinformaticsConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.executor = ParallelExecutor(ExecutionConfig(max_threads=config.threads), logger)
        
        # Initialize pipeline components
        self.mapper = ReadMapper(
            ReadMappingConfig(
                workdir=config.workdir,
                appdir=config.workdir.parent,
                threads=config.threads,
                genome_index=config.genome_index,
                splice_sites=config.splice_sites,
                compressed=config.compressed,
                keep_intermediates=config.keep_intermediates
            ),
            logger
        )
        
        self.assembler = TranscriptomeAssembler(
            AssemblyConfig(
                workdir=config.workdir,
                threads=config.threads,
                reference_gtf=config.reference_gtf,
                reference_genome=config.reference_genome,
                min_coverage=config.min_coverage,
                min_fpkm=config.min_fpkm,
                min_tpm=config.min_tpm,
                min_transcript_length=config.min_transcript_length,
                keep_intermediates=config.keep_intermediates
            ),
            logger
        )
    
    def _get_sra_metadata(self, sra_list: Path) -> pd.DataFrame:
        """Load and validate SRA metadata"""
        try:
            metadata = pd.read_csv(
                self.config.workdir / "results" / "runinfo.csv",
                dtype={'Run': str}
            )
            return metadata
        except Exception as e:
            self.logger.error("Failed to load SRA metadata")
            raise

    def process_sra_data(self, sra_list: Path) -> List[Path]:
        """Process SRA accession data"""
        self.logger.info("Starting SRA data processing")
        
        try:
            # Load metadata
            metadata = self._get_sra_metadata(sra_list)
            
            # Process each SRA accession
            bam_files = []
            for _, row in metadata.iterrows():
                accession = row['Run']
                layout = row['LibraryLayout']
                try:
                    bam_file = self.mapper.process_sample(accession, layout)
                    bam_files.append(bam_file)
                except Exception as e:
                    self.logger.error(f"Failed to process {accession}: {str(e)}")
                    raise
                    
            return bam_files
            
        except Exception as e:
            self.logger.error("SRA processing failed")
            raise

    def process_fastq_data(self, fastq_dir: Path, paired: bool = False) -> List[Path]:
        """Process local FASTQ files"""
        self.logger.info("Starting FASTQ data processing")
        
        try:
            # Get FASTQ files
            pattern = "*_[12].fastq.gz" if paired else "*.fastq.gz"
            fastq_files = list(fastq_dir.glob(pattern))
            
            if not fastq_files:
                raise ValueError(f"No FASTQ files found in {fastq_dir}")
                
            # Process each FASTQ file
            bam_files = []
            for fastq in fastq_files:
                accession = fastq.stem.split('_')[0]  # Remove _1/_2 and extension
                layout = "PAIRED" if paired else "SINGLE"
                try:
                    bam_file = self.mapper.map_reads(fastq_dir, accession, layout)
                    bam_files.append(bam_file)
                except Exception as e:
                    self.logger.error(f"Failed to process {fastq}: {str(e)}")
                    raise
                    
            return bam_files
            
        except Exception as e:
            self.logger.error("FASTQ processing failed")
            raise

    def assemble_transcriptomes(self, bam_files: List[Path]) -> Path:
        """Assemble and process transcriptomes"""
        self.logger.info("Starting transcriptome assembly")
        
        try:
            # Assemble each sample
            assemblies = []
            for bam in bam_files:
                accession = bam.stem.split('.')[0]
                try:
                    gtf = self.assembler.assemble_sample(bam, accession)
                    assemblies.append(gtf)
                except Exception as e:
                    self.logger.error(f"Failed to assemble {accession}: {str(e)}")
                    raise
                    
            # Merge assemblies
            merged_gtf = self.assembler.merge_assemblies(assemblies)
            
            # Compare to reference
            self.assembler.compare_to_reference(merged_gtf)
            
            # Extract sequences
            self.assembler.extract_sequences(merged_gtf)
            self.assembler.extract_sequences(merged_gtf, novel_only=True)
            
            # Re-estimate abundances with merged assembly
            for bam in bam_files:
                accession = bam.stem.split('.')[0]
                try:
                    self.assembler.estimate_abundance(bam, merged_gtf, accession)
                except Exception as e:
                    self.logger.error(f"Failed abundance estimation for {accession}: {str(e)}")
                    raise
                    
            return merged_gtf
            
        except Exception as e:
            self.logger.error("Transcriptome assembly failed")
            raise

    def run_pipeline(self, input_data: Path, paired: bool = False) -> None:
        """Run the complete bioinformatics pipeline"""
        try:
            # Process input data
            if input_data.is_dir():
                bam_files = self.process_fastq_data(input_data, paired)
            else:
                bam_files = self.process_sra_data(input_data)
                
            # Assemble transcriptomes
            self.assemble_transcriptomes(bam_files)
            
            self.logger.info("Bioinformatics pipeline completed successfully")
            
        except Exception as e:
            self.logger.error("Pipeline execution failed")
            raise