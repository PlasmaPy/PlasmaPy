import pytest
import warnings

from plasmapy.atomic.exceptions import AtomicError, MissingAtomicDataError, ChargeError, InvalidIonError, InvalidIsotopeError, InvalidElementError, \
    InvalidParticleError, AtomicWarning, MissingAtomicDataWarning
from plasmapy.physics.exceptions import *

plasmapy_exceptions = [
    PhysicsError,
    RelativityError,
    AtomicError,
    MissingAtomicDataError,
    ChargeError,
    InvalidIonError,
    InvalidIsotopeError,
    InvalidElementError,
    InvalidParticleError,
]


@pytest.mark.parametrize("exception", plasmapy_exceptions)
def test_exceptions(exception):
    r"""Test that custom PlasmaPy exceptions can be raised with an
    error message."""
    with pytest.raises(exception, message=f"Problem raising {exception}"):
        raise exception("What an exceptionally exceptional exception!")


plasmapy_warnings = [
    PhysicsWarning,
    RelativityWarning,
    AtomicWarning,
    MissingAtomicDataWarning,
]


@pytest.mark.parametrize("warning", plasmapy_warnings)
def test_warnings(warning):
    r"""Test that custom PlasmaPy warnings can be issued with a
    warning message."""
    with pytest.warns(warning, message=f"Problem issuing {warning}"):
        warnings.warn("Coverage decreased (-0.00002%)", warning)
