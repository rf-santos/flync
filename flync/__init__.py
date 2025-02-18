"""FLYNC v2 - FLY Long Non-Coding RNA Classifier"""

__version__ = "2.0.0"
__author__ = "Ricardo F. dos Santos"
__email__ = "ricardo.santos@nms.unl.pt"

from flync.pipeline import Pipeline
from flync.config import ConfigManager, PipelineConfig
from flync.utils.logger import setup_logger