#!/usr/bin/env python3

import subprocess
from pathlib import Path
import logging
from typing import Optional
from dataclasses import dataclass
import shutil
import os

@dataclass
class ReadMappingConfig:
    """Configuration for read mapping"""
    workdir: Path
    appdir: Path
    threads: int
    genome_index: Path
    splice_sites: Path
    compressed: bool = True
    keep_intermediates: bool = False

class ReadMapper:
    """Handles read mapping using HISAT2"""
    
    def __init__(self, config: ReadMappingConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self._validate_dependencies()
        
    def _validate_dependencies(self) -> None:
        """Check if required tools are available"""
        required_tools = ['hisat2', 'samtools', 'prefetch', 'fasterq-dump']
        for tool in required_tools:
            if not shutil.which(tool):
                raise RuntimeError(f"Required tool not found: {tool}")
        
    def _get_samtools_threads(self) -> int:
        """Calculate optimal thread count for samtools"""
        return max(1, self.config.threads // 2)
        
    def _run_command(self, cmd: list, check: bool = True) -> subprocess.CompletedProcess:
        """Run a shell command with proper error handling"""
        try:
            return subprocess.run(
                cmd,
                check=check,
                capture_output=True,
                text=True
            )
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command failed: {' '.join(cmd)}")
            self.logger.error(f"Error output: {e.stderr}")
            raise

    def download_sra(self, accession: str, layout: str) -> Path:
        """Download reads from SRA"""
        data_dir = self.config.workdir / "data" / accession
        data_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if files already exist
        if list(data_dir.glob("*.fastq.gz")):
            self.logger.info(f"Reads already downloaded for {accession}")
            return data_dir
            
        self.logger.info(f"Downloading reads for {accession}")
        
        try:
            # Download SRA data
            self._run_command([
                "prefetch",
                "-O", str(data_dir),
                accession
            ])
            
            # Convert to FASTQ
            self._run_command([
                "fasterq-dump",
                "-f",
                "-3",
                "-e", str(self.config.threads),
                "-O", str(data_dir),
                accession
            ])
            
            # Compress FASTQ files
            if self.config.compressed:
                for fastq in data_dir.glob("*.fastq"):
                    self._run_command(["gzip", str(fastq)])
                    
            return data_dir
            
        except Exception as e:
            self.logger.error(f"Failed to download reads for {accession}")
            raise
            
    def map_reads(self, input_dir: Path, accession: str, layout: str = "SINGLE") -> Path:
        """Map reads using HISAT2"""
        self.logger.info(f"Mapping reads for {accession}")
        
        # Prepare output paths
        sam_file = input_dir / f"{accession}.sam"
        bam_file = input_dir / f"{accession}.bam"
        sorted_bam = input_dir / f"{accession}.sorted.bam"
        
        # Skip if output already exists
        if sorted_bam.exists():
            self.logger.info(f"Mapped reads already exist for {accession}")
            return sorted_bam
            
        try:
            # Construct HISAT2 command
            hisat_cmd = [
                "hisat2",
                "-p", str(self.config.threads),
                "-x", str(self.config.genome_index),
                "--dta",
                "--dta-cufflinks",
                "--known-splicesite-infile", str(self.config.splice_sites),
                "-S", str(sam_file)
            ]
            
            # Add input files based on layout
            if layout == "SINGLE":
                fastq = next(input_dir.glob(f"{accession}*.fastq.gz"))
                hisat_cmd.extend(["-U", str(fastq)])
            else:
                fastq1 = next(input_dir.glob(f"{accession}_1.fastq.gz"))
                fastq2 = next(input_dir.glob(f"{accession}_2.fastq.gz"))
                hisat_cmd.extend(["-1", str(fastq1), "-2", str(fastq2)])
                
            # Run HISAT2
            self._run_command(hisat_cmd)
            
            # Convert SAM to BAM
            samthreads = self._get_samtools_threads()
            
            self._run_command([
                "samtools", "view",
                "-@", str(samthreads),
                "-b",
                "-o", str(bam_file),
                str(sam_file)
            ])
            
            # Sort BAM
            self._run_command([
                "samtools", "sort",
                "-@", str(samthreads),
                "-o", str(sorted_bam),
                str(bam_file)
            ])
            
            # Index BAM
            self._run_command([
                "samtools", "index",
                str(sorted_bam)
            ])
            
            # Clean up intermediate files
            if not self.config.keep_intermediates:
                sam_file.unlink(missing_ok=True)
                bam_file.unlink(missing_ok=True)
                if self.config.compressed:
                    for fastq in input_dir.glob(f"{accession}*.fastq.gz"):
                        fastq.unlink(missing_ok=True)
                        
            return sorted_bam
            
        except Exception as e:
            self.logger.error(f"Failed to map reads for {accession}")
            raise
            
    def process_sample(self, accession: str, layout: str) -> Path:
        """Process a single sample from download through mapping"""
        try:
            data_dir = self.download_sra(accession, layout)
            return self.map_reads(data_dir, accession, layout)
        except Exception as e:
            self.logger.error(f"Failed to process sample {accession}")
            raise