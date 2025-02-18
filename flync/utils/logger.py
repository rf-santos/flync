#!/usr/bin/env python3

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path
import sys

def setup_logger(workdir: Path, name: str = "flync") -> logging.Logger:
    """
    Set up a logger with console and rotating file handlers.
    
    Args:
        workdir: Working directory path where log file will be stored
        name: Name of the logger (default: flync)
        
    Returns:
        logging.Logger: Configured logger instance
    """
    try:
        # Create logger
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)

        # Clear any existing handlers
        logger.handlers.clear()

        # Create formatters
        console_formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%H:%M:%S'
        )
        file_formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )

        # Create console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(console_formatter)

        # Create rotating file handler (10MB max size, keep 5 backups)
        log_file = workdir / "flync.log"
        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=10*1024*1024,  # 10MB
            backupCount=5,
            encoding='utf-8'
        )
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)

        # Add handlers to logger
        logger.addHandler(console_handler)
        logger.addHandler(file_handler)

        return logger

    except Exception as e:
        # Fallback to basic console logging if setup fails
        basic_logger = logging.getLogger(name)
        basic_logger.setLevel(logging.INFO)
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        basic_logger.addHandler(handler)
        basic_logger.error(f"Failed to setup logger: {str(e)}")
        return basic_logger