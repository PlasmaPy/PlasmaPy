"""
Package of functions and classes used to develop clean, readable, and informative
code.
"""
from . import decorators
from . import roman

from plasmapy.utils.decorators.checks import (check_quantity,
                                              check_relativistic,
                                              _check_quantity,
                                              _check_relativistic)

from plasmapy.utils.decorators.converter import angular_freq_to_hz
from .decorators import preserve_signature
from .exceptions import (PlasmaPyError,
                         PhysicsError,
                         RelativityError,
                         PlasmaPyWarning,
                         PhysicsWarning,
                         CouplingWarning,
                         RelativityWarning)

from .pytest_helpers import (
    run_test,
    run_test_equivalent_calls,
    call_string,
    InconsistentTypeError,
    UnexpectedResultError,
    UnexpectedExceptionError,
    RunTestError,
    IncorrectResultError,
    MissingExceptionError,
    MissingWarningError,
    assert_can_handle_nparray,
)
