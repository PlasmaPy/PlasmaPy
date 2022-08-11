"""Tests for functionality contained in `plasmapy.formulary.frequencies`."""
import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.densities import plasma_critical_density
from plasmapy.formulary.frequencies import plasma_frequency

n_i = 5e19 * u.m**-3
ion = "e-"


def test_plasma_critical_density():
    r"""Test the plasma_critical_density function in frequencies.py."""

    omega_p = plasma_frequency(n=n_i, particle=ion)
    n_c = plasma_critical_density(omega_p)

    assert n_c.unit.is_equivalent(u.m**-3)

    assert np.isclose(n_c, 5e19 * u.m**-3)

    with pytest.raises(TypeError):
        plasma_critical_density("Lorem Ipsum")

    with pytest.raises(u.UnitTypeError):
        plasma_critical_density(1 * u.Hz)

    with pytest.raises(ValueError):
        plasma_critical_density(-1 * u.rad / u.s)

    with pytest.warns(u.UnitsWarning):
        plasma_critical_density(1)
