#!/usr/bin/env python3

import json
from pathlib import Path
from typing import Any, Dict, Optional
import logging
from dataclasses import dataclass
import hashlib
import time

@dataclass
class CheckpointData:
    """Represents a pipeline checkpoint"""
    stage: str
    data: Dict[str, Any]
    timestamp: float
    input_hash: str
    
class CheckpointManager:
    """Manages pipeline checkpoints"""
    
    def __init__(self, config: Dict[str, Any], logger: logging.Logger):
        """Initialize checkpoint manager"""
        self.enabled = config.get('enabled', True)
        self.directory = Path(config.get('directory', 'checkpoints'))
        self.force_restart = config.get('force_restart', False)
        self.logger = logger
        
        if self.enabled and not self.force_restart:
            self.directory.mkdir(parents=True, exist_ok=True)
            
    def _compute_input_hash(self, input_path: Path) -> str:
        """Compute hash of input file to detect changes"""
        hasher = hashlib.sha256()
        with open(input_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                hasher.update(chunk)
        return hasher.hexdigest()
        
    def _get_checkpoint_path(self, stage: str) -> Path:
        """Get path for checkpoint file"""
        return self.directory / f"{stage}.checkpoint"
        
    def save(self, stage: str, data: Dict[str, Any], input_path: Optional[Path] = None) -> None:
        """Save a checkpoint"""
        if not self.enabled or self.force_restart:
            return
            
        try:
            checkpoint = CheckpointData(
                stage=stage,
                data=data,
                timestamp=time.time(),
                input_hash=self._compute_input_hash(input_path) if input_path else ''
            )
            
            checkpoint_path = self._get_checkpoint_path(stage)
            with open(checkpoint_path, 'w') as f:
                # Convert checkpoint to dict for JSON serialization
                checkpoint_dict = {
                    'stage': checkpoint.stage,
                    'data': checkpoint.data,
                    'timestamp': checkpoint.timestamp,
                    'input_hash': checkpoint.input_hash
                }
                json.dump(checkpoint_dict, f, indent=2)
                
            self.logger.debug(f"Saved checkpoint for stage: {stage}")
            
        except Exception as e:
            self.logger.warning(f"Failed to save checkpoint for stage {stage}: {str(e)}")
            
    def load(self, stage: str, input_path: Optional[Path] = None) -> Optional[Dict[str, Any]]:
        """Load checkpoint data if valid"""
        if not self.enabled or self.force_restart:
            return None
            
        checkpoint_path = self._get_checkpoint_path(stage)
        if not checkpoint_path.exists():
            return None
            
        try:
            with open(checkpoint_path) as f:
                checkpoint_dict = json.load(f)
                
            # Verify input hasn't changed if path provided
            if input_path:
                current_hash = self._compute_input_hash(input_path)
                if current_hash != checkpoint_dict['input_hash']:
                    self.logger.debug(f"Input changed for stage {stage}, ignoring checkpoint")
                    return None
                    
            self.logger.info(f"Resuming from checkpoint: {stage}")
            return checkpoint_dict['data']
            
        except Exception as e:
            self.logger.warning(f"Failed to load checkpoint for stage {stage}: {str(e)}")
            return None
            
    def clear(self, stage: Optional[str] = None) -> None:
        """Clear checkpoints for specific stage or all stages"""
        if not self.enabled:
            return
            
        try:
            if stage:
                checkpoint_path = self._get_checkpoint_path(stage)
                if checkpoint_path.exists():
                    checkpoint_path.unlink()
                    self.logger.debug(f"Cleared checkpoint for stage: {stage}")
            else:
                # Clear all checkpoints
                for checkpoint_file in self.directory.glob("*.checkpoint"):
                    checkpoint_file.unlink()
                self.logger.debug("Cleared all checkpoints")
                
        except Exception as e:
            self.logger.warning(f"Failed to clear checkpoints: {str(e)}")
            
    def stage_completed(self, stage: str) -> bool:
        """Check if a stage has a valid checkpoint"""
        if not self.enabled or self.force_restart:
            return False
            
        checkpoint_path = self._get_checkpoint_path(stage)
        return checkpoint_path.exists()