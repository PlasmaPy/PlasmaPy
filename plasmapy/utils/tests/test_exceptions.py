import pytest
import warnings

from plasmapy.atomic.exceptions import (
    AtomicError,
    MissingAtomicDataError,
    ChargeError,
    InvalidIonError,
    InvalidIsotopeError,
    InvalidElementError,
    InvalidParticleError,
    AtomicWarning,
    MissingAtomicDataWarning,
    )

from ..exceptions import *

plasmapy_exceptions = [
    PlasmaPyError,
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
    with pytest.raises(exception):
        raise exception("What an exceptionally exceptional exception!")
        pytest.fail(f"Problem raising {exception}")


@pytest.mark.parametrize("exception", plasmapy_exceptions)
def test_PlasmaPyError_subclassing(exception):
    r"""Test that each custom PlasmaPy exception can be caught
    as a PlasmaPyError."""
    with pytest.raises(PlasmaPyError):
        raise exception("I'm sorry, Dave.  I'm afraid I can't do that.")
        pytest.fail(f"Problem with subclassing of {exception}")


plasmapy_warnings = [
    PlasmaPyWarning,
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


@pytest.mark.parametrize("warning", plasmapy_warnings)
def test_PlasmaPyWarning_subclassing(warning):
    r"""Test that custom PlasmaPy warnings can all be caught
    as a PlasmaPyWarning."""
    with pytest.warns(PlasmaPyWarning, message=(
            f"Problem with subclassing of {warning}")):
        warnings.warn("Electrons are WEIRD.", warning)
