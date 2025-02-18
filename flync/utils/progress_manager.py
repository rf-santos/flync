#!/usr/bin/env python3

import sys
import time
from typing import Optional, Any
from dataclasses import dataclass
from enum import Enum, auto
from contextlib import contextmanager
import logging
from tqdm import tqdm

class TaskStatus(Enum):
    """Task execution status"""
    PENDING = auto()
    RUNNING = auto()
    COMPLETED = auto()
    FAILED = auto()

@dataclass
class TaskProgress:
    """Represents the progress of a pipeline task"""
    name: str
    desc: Optional[str] = None
    total: int = 100
    unit: str = '%'
    status: TaskStatus = TaskStatus.PENDING

class ProgressManager:
    """Manages progress reporting and user feedback"""
    
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self._current_progress: Optional[tqdm] = None
        self._start_time = time.time()
        self._task_history = []
        self._current_task = None
        self._error_count = 0
        self._warning_count = 0
        
    def _format_time(self, seconds: float) -> str:
        """Format time duration in a human-readable format"""
        minutes, seconds = divmod(int(seconds), 60)
        hours, minutes = divmod(minutes, 60)
        return f"{hours:02d}:{minutes:02d}:{seconds:02d}"
        
    def start_task(self, task: TaskProgress) -> None:
        """Start tracking progress for a task"""
        self._current_task = task
        task.status = TaskStatus.RUNNING
        self._task_history.append(task)
        
        self.logger.info(f"Starting task: {task.name}")
        if task.desc:
            self.logger.info(f"Description: {task.desc}")
            
        self._current_progress = tqdm(
            total=task.total,
            desc=task.desc or task.name,
            unit=task.unit,
            file=sys.stdout
        )

    def update(self, amount: int = 1) -> None:
        """Update progress for current task"""
        if self._current_progress:
            self._current_progress.update(amount)

    def set_progress(self, value: int) -> None:
        """Set absolute progress value for current task"""
        if self._current_progress:
            self._current_progress.n = value
            self._current_progress.refresh()

    def complete_task(self, success: bool = True) -> None:
        """Mark current task as complete"""
        if self._current_task:
            self._current_task.status = TaskStatus.COMPLETED if success else TaskStatus.FAILED
            if self._current_progress:
                self._current_progress.close()
                self._current_progress = None
            self.logger.info(f"Task {self._current_task.name} {'completed' if success else 'failed'}")

    @contextmanager
    def task_progress(self, task: TaskProgress):
        """Context manager for task progress tracking"""
        try:
            self.start_task(task)
            yield self
        except Exception as e:
            self.complete_task(success=False)
            self.error(f"Task failed: {task.name}", e)
            raise
        else:
            self.complete_task(success=True)

    def print_summary(self) -> None:
        """Print execution summary"""
        total_time = time.time() - self._start_time
        
        # Calculate statistics
        total_tasks = len(self._task_history)
        completed_tasks = sum(1 for t in self._task_history if t.status == TaskStatus.COMPLETED)
        failed_tasks = sum(1 for t in self._task_history if t.status == TaskStatus.FAILED)
        
        # Print summary
        self.logger.info("-" * 40)
        self.logger.info("Pipeline Execution Summary")
        self.logger.info("-" * 40)
        self.logger.info(f"Total execution time: {self._format_time(total_time)}")
        self.logger.info(f"Tasks completed: {completed_tasks}/{total_tasks}")
        if failed_tasks > 0:
            self.logger.warning(f"Tasks failed: {failed_tasks}")
        self.logger.info(f"Errors encountered: {self._error_count}")
        self.logger.info(f"Warnings raised: {self._warning_count}")
        self.logger.info("-" * 40)

    def error(self, message: str, exception: Optional[Exception] = None) -> None:
        """Log error message with proper formatting"""
        self._error_count += 1
        self.logger.error(message)
        if exception:
            self.logger.error(f"Error details: {str(exception)}")
            if hasattr(exception, '__traceback__'):
                self.logger.debug("Stack trace:", exc_info=True)

    def warning(self, message: str) -> None:
        """Log warning message with proper formatting"""
        self._warning_count += 1
        self.logger.warning(message)

    def debug(self, message: str) -> None:
        """Log debug message"""
        self.logger.debug(message)

    def info(self, message: str) -> None:
        """Log info message"""
        self.logger.info(message)

    def status(self, message: str) -> None:
        """Print status message with timestamp"""
        self.logger.info(f"[{time.strftime('%H:%M:%S')}] {message}")

    def get_task_status(self, task_name: str) -> Optional[TaskStatus]:
        """Get the status of a specific task by name"""
        for task in self._task_history:
            if task.name == task_name:
                return task.status
        return None