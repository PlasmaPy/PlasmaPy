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
]

from plasmapy.utils import (
    decorators,
    datatype_factory_base,
    error_messages,
    exceptions,
    roman,
)
from plasmapy.utils.error_messages import call_string
from plasmapy.utils.exceptions import (
    CouplingWarning,
    PhysicsError,
    PhysicsWarning,
    PlasmaPyError,
    PlasmaPyWarning,
    RelativityError,
    RelativityWarning,
)

try:
    from plasmapy.utils import pytest_helpers
except ModuleNotFoundError:
    # pytest is not a hard dependency, so only import pytest_helpers is pytest
    # is installed
    pass
