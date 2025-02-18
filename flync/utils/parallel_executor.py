#!/usr/bin/env python3

import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from typing import Callable, List, Any, Dict, Union, Iterator
from pathlib import Path
from dataclasses import dataclass
import logging
import psutil
import threading
from queue import Queue
import time

@dataclass
class ExecutionConfig:
    """Configuration for parallel execution"""
    max_threads: int = 0  # 0 means auto-detect
    max_processes: int = 0  # 0 means auto-detect
    io_bound: bool = True  # True for I/O-bound tasks, False for CPU-bound
    chunk_size: int = 10  # Size of chunks for parallel processing

class ParallelExecutor:
    """Handles parallel execution of pipeline tasks"""
    
    def __init__(self, config: ExecutionConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self._thread_local = threading.local()
        
    def _get_available_threads(self) -> int:
        """Determine number of threads to use based on system resources"""
        if self.config.max_threads > 0:
            return min(self.config.max_threads, mp.cpu_count())
            
        cpu_count = mp.cpu_count()
        mem = psutil.virtual_memory()
        
        # Use memory info to avoid oversubscription
        mem_per_thread = 2 * 1024 * 1024 * 1024  # 2GB per thread
        mem_threads = max(1, int(mem.available / mem_per_thread))
        
        # Account for system load
        load = psutil.getloadavg()[0]
        load_threads = max(1, int(cpu_count - load))
        
        optimal_threads = min(cpu_count, mem_threads, load_threads)
        self.logger.debug(f"Optimal thread count: {optimal_threads} (CPU: {cpu_count}, Mem: {mem_threads}, Load: {load_threads})")
        
        return optimal_threads
        
    def _get_available_processes(self) -> int:
        """Determine number of processes to use based on system resources"""
        if self.config.max_processes > 0:
            return min(self.config.max_processes, mp.cpu_count())
            
        # For CPU-bound tasks, use N-1 processes by default
        return max(1, mp.cpu_count() - 1)

    def process_map(self, func: Callable, items: List[Any], **kwargs) -> List[Any]:
        """Execute function on items using process pool"""
        n_processes = self._get_available_processes()
        results = []
        errors = []
        
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            futures = [executor.submit(func, item, **kwargs) for item in items]
            
            for future in as_completed(futures):
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    self.logger.error(f"Process execution failed: {str(e)}")
                    errors.append(e)
                    
        if errors:
            self.logger.error(f"{len(errors)} tasks failed during parallel processing")
            raise RuntimeError("Parallel processing failed")
            
        return results

    def thread_map(self, func: Callable, items: List[Any], **kwargs) -> List[Any]:
        """Execute function on items using thread pool"""
        n_threads = self._get_available_threads()
        results = []
        errors = Queue()
        
        def worker(item):
            try:
                if not hasattr(self._thread_local, 'initialized'):
                    self._thread_local.initialized = True
                return func(item, **kwargs)
            except Exception as e:
                errors.put(e)
                return None
                
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            futures = [executor.submit(worker, item) for item in items]
            results = [f.result() for f in futures]
            
        if not errors.empty():
            error = errors.get()
            self.logger.error(f"Thread execution failed: {str(error)}")
            raise error
            
        return [r for r in results if r is not None]

    def parallel_map(self, func: Callable, items: List[Any], use_processes: bool = True, **kwargs) -> List[Any]:
        """Execute function on items in parallel using either processes or threads"""
        if not items:
            return []
            
        executor = self.process_map if use_processes else self.thread_map
        return executor(func, items, **kwargs)

    @staticmethod
    def chunk_list(items: List[Any], chunk_size: int) -> List[List[Any]]:
        """Split list into chunks for parallel processing"""
        return [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]

    def parallel_io(self, func: Callable, files: List[Path]) -> List[Any]:
        """Optimized parallel execution for I/O operations"""
        # Use threads for I/O operations since GIL impact is minimal
        chunks = self.chunk_list(files, self.config.chunk_size)
        results = []
        
        for chunk in chunks:
            chunk_results = self.thread_map(func, chunk)
            results.extend(chunk_results)
            
        return results

    def parallel_cpu(self, func: Callable, items: List[Any]) -> List[Any]:
        """Optimized parallel execution for CPU-intensive operations"""
        # Use processes for CPU-bound tasks to bypass GIL
        chunks = self.chunk_list(items, self.config.chunk_size)
        results = []
        
        for chunk in chunks:
            chunk_results = self.process_map(func, chunk)
            results.extend(chunk_results)
            
        return results