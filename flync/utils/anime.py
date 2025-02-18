#!/usr/bin/env python3

import sys
import time
import threading
from enum import Enum
from typing import List, Optional
import shutil

class AnimationType(Enum):
    """Available animation types"""
    DOTS = 'dots'
    BRAILLE = 'braille'
    ARROW = 'arrow'
    PULSE = 'pulse'
    METRO = 'metro'

class ProgressAnimation:
    """Handles animated progress indicators"""
    
    ANIMATIONS = {
        AnimationType.DOTS: ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏'],
        AnimationType.BRAILLE: ['⣾', '⣽', '⣻', '⢿', '⡿', '⣟', '⣯', '⣷'],
        AnimationType.ARROW: ['▹▹▹▹▹', '▸▹▹▹▹', '▹▸▹▹▹', '▹▹▸▹▹', '▹▹▹▸▹', '▹▹▹▹▸'],
        AnimationType.PULSE: ['█▁▁▁▁', '██▁▁▁', '███▁▁', '████▁', '█████', '████▁', '███▁▁', '██▁▁▁', '█▁▁▁▁'],
        AnimationType.METRO: ['[    ]', '[=   ]', '[==  ]', '[=== ]', '[ ===]', '[  ==]', '[   =]']
    }

    def __init__(self, animation_type: AnimationType = AnimationType.DOTS, delay: float = 0.1):
        self.animation = self.ANIMATIONS[animation_type]
        self.delay = delay
        self._running = False
        self._thread: Optional[threading.Thread] = None
        self._current_frame = 0
        self._terminal_width = shutil.get_terminal_size().columns
        
    def _animate(self):
        """Animation loop"""
        while self._running:
            frame = self.animation[self._current_frame]
            sys.stdout.write('\r' + frame)
            sys.stdout.flush()
            time.sleep(self.delay)
            self._current_frame = (self._current_frame + 1) % len(self.animation)
            
    def start(self):
        """Start the animation"""
        self._running = True
        self._thread = threading.Thread(target=self._animate)
        self._thread.daemon = True
        self._thread.start()
        
    def stop(self):
        """Stop the animation"""
        self._running = False
        if self._thread:
            self._thread.join()
        sys.stdout.write('\r' + ' ' * self._terminal_width + '\r')
        sys.stdout.flush()

    def __enter__(self):
        self.start()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()

class ProgressBar:
    """Progress bar with percentage and optional description"""
    
    def __init__(self, total: int, desc: str = '', width: int = 40):
        self.total = total
        self.desc = desc
        self.width = width
        self.current = 0
        self._start_time = time.time()
        
    def update(self, amount: int = 1):
        """Update progress by amount"""
        self.current += amount
        self._display()
        
    def _display(self):
        """Display the progress bar"""
        percentage = min(100, int(self.current / self.total * 100))
        filled = int(self.width * self.current / self.total)
        bar = '█' * filled + '▁' * (self.width - filled)
        
        elapsed = time.time() - self._start_time
        if self.current > 0:
            eta = elapsed * (self.total / self.current - 1)
            eta_str = f'ETA: {int(eta)}s'
        else:
            eta_str = 'ETA: --'
            
        print(f'\r{self.desc}: |{bar}| {percentage}% {eta_str}', end='')
        if self.current >= self.total:
            print()  # New line on completion