#!/usr/bin/env python3

import subprocess
from pathlib import Path
import logging
from typing import List, Optional, Dict, Tuple
from dataclasses import dataclass
import shutil

@dataclass
class AssemblyConfig:
    """Configuration for transcriptome assembly"""
    workdir: Path
    threads: int
    reference_gtf: Path
    reference_genome: Path
    min_coverage: float = 2.5
    min_fpkm: float = 0.1
    min_tpm: float = 0.1
    min_transcript_length: int = 200
    keep_intermediates: bool = False

class TranscriptomeAssembler:
    """Handles transcriptome assembly using StringTie"""
    
    def __init__(self, config: AssemblyConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self._assembly_dir = config.workdir / "assemblies" / "stringtie"
        self._cuffcompare_dir = config.workdir / "cuffcompare"
        self._validate_dependencies()
        
    def _validate_dependencies(self) -> None:
        """Check if required tools are available"""
        required_tools = ['stringtie', 'cuffcompare', 'gffread']
        for tool in required_tools:
            if not shutil.which(tool):
                raise RuntimeError(f"Required tool not found: {tool}")
            
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

    def assemble_sample(self, bam_file: Path, accession: str) -> Path:
        """Assemble transcripts for a single sample"""
        self._assembly_dir.mkdir(parents=True, exist_ok=True)
        output_gtf = self._assembly_dir / f"{accession}.rna.gtf"
        
        if output_gtf.exists():
            self.logger.info(f"Assembly already exists for {accession}")
            return output_gtf
            
        self.logger.info(f"Assembling transcripts for {accession}")
        
        try:
            # Run StringTie
            cmd = [
                "stringtie",
                str(bam_file),
                "-G", str(self.config.reference_gtf),
                "-o", str(output_gtf),
                "-p", str(self.config.threads),
                "-c", str(self.config.min_coverage),
                "-f", str(self.config.min_fpkm),
                "-T", str(self.config.min_tpm),
                "-m", str(self.config.min_transcript_length)
            ]
            
            self._run_command(cmd)
            return output_gtf
            
        except Exception as e:
            self.logger.error(f"Failed to assemble transcripts for {accession}")
            raise

    def merge_assemblies(self, assembly_files: List[Path], output_name: str = "merged") -> Path:
        """Merge multiple transcript assemblies"""
        if not assembly_files:
            raise ValueError("No assembly files provided for merging")
            
        merged_gtf = self._assembly_dir / f"{output_name}.gtf"
        
        if merged_gtf.exists():
            self.logger.info("Merged assembly already exists")
            return merged_gtf
            
        self.logger.info("Merging transcript assemblies")
        
        try:
            # Create manifest file for StringTie
            manifest = self._assembly_dir / "merge_manifest.txt"
            with open(manifest, "w") as f:
                for assembly in assembly_files:
                    f.write(f"{assembly}\n")
            
            # Run StringTie merge
            cmd = [
                "stringtie",
                "--merge",
                "-G", str(self.config.reference_gtf),
                "-o", str(merged_gtf),
                "-p", str(self.config.threads),
                str(manifest)
            ]
            
            self._run_command(cmd)
            
            # Clean up
            if not self.config.keep_intermediates:
                manifest.unlink()
                
            return merged_gtf
            
        except Exception as e:
            self.logger.error("Failed to merge assemblies")
            raise

    def estimate_abundance(self, bam_file: Path, merged_gtf: Path, accession: str) -> Path:
        """Re-estimate transcript abundances using merged assembly"""
        abundance_gtf = self._assembly_dir / f"{accession}.abundance.gtf"
        
        if abundance_gtf.exists():
            self.logger.info(f"Abundance estimates already exist for {accession}")
            return abundance_gtf
            
        self.logger.info(f"Estimating transcript abundances for {accession}")
        
        try:
            cmd = [
                "stringtie",
                str(bam_file),
                "-G", str(merged_gtf),
                "-o", str(abundance_gtf),
                "-e",  # Only estimate abundance
                "-p", str(self.config.threads),
                "-A", str(abundance_gtf).replace(".gtf", ".tab")
            ]
            
            self._run_command(cmd)
            return abundance_gtf
            
        except Exception as e:
            self.logger.error(f"Failed to estimate abundances for {accession}")
            raise

    def compare_to_reference(self, merged_gtf: Path) -> Tuple[Path, Path]:
        """Compare assembled transcripts to reference annotation"""
        self._cuffcompare_dir.mkdir(parents=True, exist_ok=True)
        output_prefix = self._cuffcompare_dir / "cuffcomp"
        
        if output_prefix.with_suffix('.gtf').exists():
            self.logger.info("Transcript comparison already exists")
            return output_prefix.with_suffix('.gtf'), output_prefix.with_suffix('.tracking')
            
        self.logger.info("Comparing assembly to reference annotation")
        
        try:
            cmd = [
                "cuffcompare",
                "-R",  # Show matching ref transcripts
                "-r", str(self.config.reference_gtf),
                str(merged_gtf),
                "-o", str(output_prefix)
            ]
            
            self._run_command(cmd)
            return output_prefix.with_suffix('.gtf'), output_prefix.with_suffix('.tracking')
            
        except Exception as e:
            self.logger.error("Failed to compare transcripts")
            raise

    def extract_sequences(self, merged_gtf: Path, novel_only: bool = False) -> Path:
        """Extract transcript sequences from assembly"""
        assemblies_dir = self.config.workdir / "assemblies"
        
        if novel_only:
            # Filter for novel transcripts (MSTRG prefix)
            filtered_gtf = assemblies_dir / "merged-new-transcripts.gtf"
            if not filtered_gtf.exists():
                self.logger.info("Filtering novel transcripts")
                with open(merged_gtf) as f_in, open(filtered_gtf, 'w') as f_out:
                    for line in f_in:
                        if 'MSTRG' in line:
                            f_out.write(line)
            output_fa = assemblies_dir / "assembled-new-transcripts.fa"
            source_gtf = filtered_gtf
        else:
            output_fa = assemblies_dir / "assembled-transcripts.fa"
            source_gtf = merged_gtf
            
        if output_fa.exists():
            self.logger.info(f"Sequences already extracted: {output_fa}")
            return output_fa
            
        self.logger.info(f"Extracting sequences from {source_gtf}")
        
        try:
            cmd = [
                "gffread",
                "-w", str(output_fa),
                "-g", str(self.config.reference_genome),
                str(source_gtf)
            ]
            
            self._run_command(cmd)
            return output_fa
            
        except Exception as e:
            self.logger.error("Failed to extract sequences")
            raise