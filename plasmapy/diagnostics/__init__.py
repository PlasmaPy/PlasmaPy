"""The diagnostics subpackage contains tools for experimental research.
Currently, we have functionality for analyzing data from Langmuir probes.
"""
__all__ = ['AbstractProbe']

import abc


class AbstractProbe(abc.ABC):
    """Abstract class for defining probe characteristics/parameters."""
