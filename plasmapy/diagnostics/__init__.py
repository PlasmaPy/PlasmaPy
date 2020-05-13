"""The diagnostics subpackage contains tools for experimental research.
Currently, we have functionality for analyzing data from Langmuir probes.
"""
__all__ = [
    "AbstractProbe",
    "XAbstractDiagnostic",
    "xdiagnostics",
    "core",
]

from . import core
from .core import AbstractProbe, XAbstractDiagnostic, XDiagnostics

xdiagnostics = XDiagnostics()
"""
An instance of `~plasmapy.diagnostics.core.XDiagnostic` that allows users to only
load the xarray diagnostic accessors they want to use.

Examples
--------
>>> import xarray as xr
>>> from plasmapy.diagnostics import xdiagnostics
...
>>> xdiagnostics
<plasmapy.diagnostics.core.XDiagnostics>
Enabled   Available Diagnostic
    [ ]   swept_langmuir
...
>>> xdiagnostics.enable('swept_langmuir')
>>> xdiagnostics
<plasmapy.diagnostics.core.XDiagnostics>
Enabled   Available Diagnostic
    [x]   swept_langmuir
...
>>> hasattr(xr.Dataset, 'swept_langmuir')
True
"""
