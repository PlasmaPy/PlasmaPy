"""Tests for functionality contained in `plasmapy.formulary.frequencies`."""
import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.densities import plasma_critical_density
from plasmapy.formulary.frequencies import plasma_frequency


class TestPlasmaCriticalDensity:
    """Test the plasma_critical_density function in densities.py."""

    n_i = 5e19 * u.m**-3
    omega_p = plasma_frequency(n=n_i, particle="e-")

    @pytest.fixture()
    def get_return(self):
        """Get the critical density for the example frequency"""
        return plasma_critical_density(self.omega_p)

    def test_units(self, get_return):
        """Test the return units"""
        n_c = get_return

        assert n_c.unit.is_equivalent(u.m**-3)

    def test_value(self, get_return):
        """Compare the calculated critical density with the expected value"""
        n_c = get_return

        assert np.isclose(n_c.value, self.n_i.value)
