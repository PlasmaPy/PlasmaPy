import pytest
import warnings

from .. import (PlasmaPyError,
                PhysicsError,
                RelativityError,
                AtomicError,
                IonError,
                IsotopeError,
                ElementError,
                ParticleError)

from .. import (PlasmaPyWarning,
                PhysicsWarning,
                RelativityWarning,
                AtomicWarning)


plasmapy_exceptions = [
    PlasmaPyError,
    PhysicsError,
    RelativityError,
    AtomicError,
    IonError,
    IsotopeError,
    ElementError,
    ParticleError,
]

plasmapy_warnings = [
    PlasmaPyWarning,
    PhysicsWarning,
    RelativityWarning,
    AtomicWarning,
]


@pytest.mark.parametrize("exception", plasmapy_exceptions)
def test_exceptions(exception):
    r"""Test that custom PlasmaPy exceptions can be raised with an
    error message."""
    with pytest.raises(exception):
        raise exception("What an exceptionally exceptional exception!")


@pytest.mark.parametrize("warning", plasmapy_warnings)
def test_warnings(warning):
    r"""Test that custom PlasmaPy warnings can be issued with a
    warning message."""
    with pytest.warns(warning):
        warnings.warn("Coverage decreased (-0.00002%)", warning)


@pytest.mark.parametrize("exception", plasmapy_exceptions)
def test_PlasmaPyError_subclassing(exception):
    r"""Test that each custom PlasmaPy exception can be caught
    as a PlasmaPyError."""
    with pytest.raises(PlasmaPyError):
        raise exception("I'm sorry, Dave.  I'm afraid I can't do that.")


@pytest.mark.parametrize("warning", plasmapy_warnings)
def test_PlasmaPyWarning_subclassing(warning):
    r"""Test that custom PlasmaPy warnings can all be caught
    as a PlasmaPyWarning."""
    with pytest.warns(PlasmaPyWarning):
        warnings.warn("Electrons are WEIRD.", warning)
