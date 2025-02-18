#!/usr/bin/env python3

import subprocess
from pathlib import Path
import argparse
from argparse import RawTextHelpFormatter
import sys
from scripts.logger import setup_logger
from scripts.config_manager import ConfigManager, create_example_config
from scripts.pipeline_controller import Pipeline

# Constants
AUTHOR = 'Ricardo F. dos Santos'
MAIL = '<ricardo.santos@nms.unl.pt>'
VERSION = '2.0.0'

class PipelineError(Exception):
    """Base class for pipeline exceptions"""
    pass

class ConfigurationError(PipelineError):
    """Raised when there is an error in configuration"""
    pass

def main():
    """Command-line entry point"""
    parser = argparse.ArgumentParser(
        prog='flync',
        description=f'Fly Non-Coding RNA discovery and classification\nVersion: {VERSION}\nAuthor: {AUTHOR}\nContact: {MAIL}',
        formatter_class=RawTextHelpFormatter
    )
    
    parser.add_argument('--version', action='version', version=f'%(prog)s {VERSION}')
    
    subparsers = parser.add_subparsers(
        title='commands',
        description='Run flync <command> --help to see available options.',
        dest='command',
        required=True
    )
    
    # Run pipeline command
    run_parser = subparsers.add_parser(
        'run',
        help='Run pipeline using configuration file'
    )
    run_parser.add_argument(
        '--config',
        '-c',
        type=Path,
        default='config.yaml',
        help='Configuration file path (default: config.yaml)'
    )
    
    # Init command to create example config
    init_parser = subparsers.add_parser(
        'init',
        help='Initialize a new FLYNC project'
    )
    init_parser.add_argument(
        '--output',
        '-o',
        type=Path,
        default='config.yaml',
        help='Output path for example configuration (default: config.yaml)'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    try:
        if args.command == 'run':
            if not args.config.exists():
                raise ConfigurationError(f"Configuration file not found: {args.config}")
            
            pipeline = Pipeline(args.config)
            pipeline.run()
            
        elif args.command == 'init':
            create_example_config(args.output)
            print(f"Created example configuration at: {args.output}")
            
    except ConfigurationError as e:
        print(f"Configuration error: {str(e)}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()