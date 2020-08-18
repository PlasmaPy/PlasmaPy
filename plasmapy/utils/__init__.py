"""
Package of functions and classes used to develop clean, readable, and informative
code.
"""
__all__ = [
    "CouplingWarning",
    "PhysicsError",
    "PhysicsWarning",
    "PlasmaPyError",
    "PlasmaPyWarning",
    "RelativityError",
    "RelativityWarning",
    "call_string",
    "decorators",
    "exceptions",
    "roman",
]

from plasmapy.utils import (
    datatype_factory_base,
    decorators,
    exceptions,
    roman,
)
from plasmapy.utils.formatting import call_string
from plasmapy.utils.exceptions import (
    CouplingWarning,
    PhysicsError,
    PhysicsWarning,
    PlasmaPyError,
    PlasmaPyWarning,
    RelativityError,
    RelativityWarning,
)

from . import formatting

try:
    from plasmapy.utils import pytest_helpers
except ModuleNotFoundError:
    # pytest is not a hard dependency, so only import pytest_helpers is pytest
    # is installed
    pass
