"""Small public API for running AmPrime from Python."""

from .api import (
    AmPrimeProject,
    FunctionalTestResult,
    PipelineRun,
    ResultPaths,
)

__all__ = [
    "AmPrimeProject",
    "FunctionalTestResult",
    "PipelineRun",
    "ResultPaths",
]
