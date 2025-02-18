#!/usr/bin/env python3

import yaml
from pathlib import Path
from typing import Dict, Optional, Union, Any
from dataclasses import dataclass
import logging

@dataclass
class MapperConfig:
    """Read mapping configuration"""
    aligner: str = 'hisat2'  # Default aligner
    threads: int = 1
    mismatch_tolerance: int = 2
    min_alignment_score: int = 20
    min_intron_length: int = 20
    max_intron_length: int = 500000

@dataclass 
class AssemblerConfig:
    """Transcript assembly configuration"""
    assembler: str = 'stringtie'  # Default assembler
    threads: int = 1
    min_coverage: float = 2.5
    min_isoform_fraction: float = 0.1
    min_transcript_length: int = 200

@dataclass
class LoggingConfig:
    """Logging configuration"""
    level: str = 'INFO'
    file: Optional[Path] = None
    console: bool = True
    
@dataclass
class CheckpointConfig:
    """Checkpointing configuration"""
    enabled: bool = True
    directory: Optional[Path] = None
    force_restart: bool = False

@dataclass
class PipelineConfig:
    """Pipeline configuration parameters"""
    output: Path
    threads: int
    metadata: Optional[Path] = None
    logging: LoggingConfig = LoggingConfig()
    checkpoint: CheckpointConfig = CheckpointConfig()
    fastq_active: bool = False
    fastq_path: Optional[Path] = None 
    fastq_paired: bool = False
    sra_path: Optional[Path] = None
    mapper: MapperConfig = MapperConfig()
    assembler: AssemblerConfig = AssemblerConfig()
    genome_release: str = 'Ensembl BDGP6.32'

class ConfigManager:
    """Manages configuration loading and validation"""
    
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        
    def _validate_path(self, path: Union[str, Path], required: bool = True, must_exist: bool = True) -> Optional[Path]:
        """Validate a path configuration value"""
        if path:
            try:
                path = Path(path).resolve()
                if must_exist and not path.exists():
                    raise ValueError(f"Path does not exist: {path}")
                return path
            except Exception as e:
                if required:
                    raise ValueError(f"Invalid path: {path}. Error: {str(e)}")
                self.logger.warning(f"Invalid path: {path}. Error: {str(e)}")
                return None
        elif required:
            raise ValueError("Required path not provided")
        return None

    def _validate_threads(self, threads: Any) -> int:
        """Validate threads configuration value"""
        try:
            threads = int(threads)
            if threads < 1:
                raise ValueError("Threads must be >= 1")
            return threads
        except ValueError:
            raise ValueError("Invalid threads value")

    def _validate_logging(self, config: Dict) -> LoggingConfig:
        """Validate logging configuration"""
        log_config = LoggingConfig()
        
        if 'logging' in config:
            logging_dict = config['logging']
            if 'level' in logging_dict:
                level = logging_dict['level'].upper()
                if level not in ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'):
                    raise ValueError(f"Invalid logging level: {level}")
                log_config.level = level
            
            if 'file' in logging_dict:
                log_config.file = self._validate_path(logging_dict['file'], required=False, must_exist=False)
            
            if 'console' in logging_dict:
                log_config.console = bool(logging_dict['console'])
                
        return log_config

    def _validate_checkpoint(self, config: Dict) -> CheckpointConfig:
        """Validate checkpoint configuration"""
        cp_config = CheckpointConfig()
        
        if 'checkpoint' in config:
            cp_dict = config['checkpoint']
            if 'enabled' in cp_dict:
                cp_config.enabled = bool(cp_dict['enabled'])
            
            if 'directory' in cp_dict:
                cp_config.directory = self._validate_path(cp_dict['directory'], required=False, must_exist=False)
                
            if 'force_restart' in cp_dict:
                cp_config.force_restart = bool(cp_dict['force_restart'])
                
        return cp_config

    def _validate_mapper(self, config: Dict) -> MapperConfig:
        """Validate read mapper configuration"""
        mapper_config = MapperConfig()
        
        if 'mapper' in config:
            mapper_dict = config['mapper']
            valid_aligners = ('hisat2', 'star')
            
            if 'aligner' in mapper_dict:
                aligner = mapper_dict['aligner'].lower()
                if aligner not in valid_aligners:
                    raise ValueError(f"Invalid aligner: {aligner}. Must be one of {valid_aligners}")
                mapper_config.aligner = aligner
                
            mapper_config.threads = mapper_dict.get('threads', mapper_config.threads)
            mapper_config.mismatch_tolerance = mapper_dict.get('mismatch_tolerance', mapper_config.mismatch_tolerance)
            mapper_config.min_alignment_score = mapper_dict.get('min_alignment_score', mapper_config.min_alignment_score)
            mapper_config.min_intron_length = mapper_dict.get('min_intron_length', mapper_config.min_intron_length)
            mapper_config.max_intron_length = mapper_dict.get('max_intron_length', mapper_config.max_intron_length)
            
        return mapper_config

    def _validate_assembler(self, config: Dict) -> AssemblerConfig:
        """Validate transcript assembler configuration"""
        assembler_config = AssemblerConfig()
        
        if 'assembler' in config:
            assembler_dict = config['assembler']
            valid_assemblers = ('stringtie', 'cufflinks')
            
            if 'assembler' in assembler_dict:
                assembler = assembler_dict['assembler'].lower()
                if assembler not in valid_assemblers:
                    raise ValueError(f"Invalid assembler: {assembler}. Must be one of {valid_assemblers}")
                assembler_config.assembler = assembler
                
            assembler_config.threads = assembler_dict.get('threads', assembler_config.threads)
            assembler_config.min_coverage = float(assembler_dict.get('min_coverage', assembler_config.min_coverage))
            assembler_config.min_isoform_fraction = float(assembler_dict.get('min_isoform_fraction', assembler_config.min_isoform_fraction))
            assembler_config.min_transcript_length = int(assembler_dict.get('min_transcript_length', assembler_config.min_transcript_length))
            
        return assembler_config

    def load_config(self, config_path: Path) -> PipelineConfig:
        """Load and validate configuration from YAML file"""
        self.logger.info(f"Loading configuration from {config_path}")
        
        try:
            with open(config_path) as f:
                config_data = yaml.safe_load(f)

            # Validate required fields
            output = self._validate_path(config_data.get('output'), required=True, must_exist=False)
            threads = self._validate_threads(config_data.get('threads', 1))

            # Create output directory
            output.mkdir(parents=True, exist_ok=True)

            # Validate optional fields
            metadata = self._validate_path(config_data.get('metadata'), required=False)
            fastq_active = bool(config_data.get('fastq_active', False))
            
            config = PipelineConfig(
                output=output,
                threads=threads,
                metadata=metadata,
                fastq_active=fastq_active,
                logging=self._validate_logging(config_data),
                checkpoint=self._validate_checkpoint(config_data),
                mapper=self._validate_mapper(config_data),
                assembler=self._validate_assembler(config_data)
            )

            # Validate pipeline-specific parameters
            if fastq_active:
                config.fastq_path = self._validate_path(config_data.get('fastq_path'), required=True)
                config.fastq_paired = bool(config_data.get('fastq_paired', False))
            else:
                config.sra_path = self._validate_path(config_data.get('sra'), required=True)

            self.logger.info("Configuration loaded successfully")
            return config

        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing YAML configuration: {str(e)}")
        except Exception as e:
            raise ValueError(f"Error loading configuration: {str(e)}")

    def generate_default_config(self, output_path: Path) -> None:
        """Generate a default configuration file"""
        default_config = {
            'output': 'results',
            'threads': 2,
            'metadata': None,
            'fastq_active': False,
            'fastq_path': None,
            'fastq_paired': False, 
            'sra': None,
            'logging': {
                'level': 'INFO',
                'file': 'flync.log',
                'console': True
            },
            'checkpoint': {
                'enabled': True,
                'directory': 'checkpoints',
                'force_restart': False
            },
            'mapper': {
                'aligner': 'hisat2',
                'threads': 2,
                'mismatch_tolerance': 2,
                'min_alignment_score': 20,
                'min_intron_length': 20,
                'max_intron_length': 500000
            },
            'assembler': {
                'assembler': 'stringtie',
                'threads': 2,
                'min_coverage': 2.5,
                'min_isoform_fraction': 0.1,
                'min_transcript_length': 200
            }
        }

        try:
            with open(output_path, 'w') as f:
                yaml.safe_dump(default_config, f, default_flow_style=False, sort_keys=False)
            self.logger.info(f"Default configuration written to {output_path}")
        except Exception as e:
            raise IOError(f"Error writing default configuration: {str(e)}")

def create_example_config(output_path: Path) -> None:
    """Create an example configuration file with comments"""
    example_config = '''# FLYNC Pipeline Configuration

# Output directory for results
output: results

# Number of threads to use for parallel processing
threads: 2

# Optional metadata file for differential expression analysis
metadata: null

# Pipeline mode configuration
fastq_active: false  # Set to true to use local FASTQ files instead of SRA

# FASTQ pipeline configuration (only used if fastq_active: true)
fastq_path: null     # Directory containing FASTQ files
fastq_paired: false  # Set to true for paired-end reads

# SRA pipeline configuration (only used if fastq_active: false)
sra: null  # File containing SRA accession numbers

# Logging configuration
logging:
  level: INFO       # DEBUG, INFO, WARNING, ERROR, or CRITICAL
  file: flync.log   # Log file path (optional)
  console: true     # Enable console output

# Checkpointing configuration
checkpoint:
  enabled: true           # Enable/disable checkpointing
  directory: checkpoints  # Directory to store checkpoints
  force_restart: false    # Force restart from beginning

# Read mapper configuration
mapper:
  aligner: hisat2           # Read aligner (hisat2 or star)
  threads: 2                # Threads for read mapping
  mismatch_tolerance: 2     # Number of allowed mismatches
  min_alignment_score: 20   # Minimum alignment score
  min_intron_length: 20     # Minimum intron length
  max_intron_length: 500000 # Maximum intron length

# Transcript assembler configuration  
assembler:
  assembler: stringtie          # Assembler (stringtie or cufflinks)
  threads: 2                    # Threads for assembly
  min_coverage: 2.5            # Minimum read coverage
  min_isoform_fraction: 0.1    # Minimum isoform fraction
  min_transcript_length: 200    # Minimum transcript length
'''
    
    try:
        with open(output_path, 'w') as f:
            f.write(example_config)
    except Exception as e:
        raise IOError(f"Error writing example configuration: {str(e)}")