"""Small public API for running AmPair from Python."""

from .api import (
    AmPairProject,
    FunctionalTestResult,
    PipelineRun,
    ResultPaths,
)

__all__ = [
    "AmPairProject",
    "FunctionalTestResult",
    "PipelineRun",
    "ResultPaths",
]
