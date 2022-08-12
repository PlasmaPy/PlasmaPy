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
    def n_c(self):
        """Get the critical density for the example frequency"""
        return plasma_critical_density(self.omega_p)

    def test_units(self, n_c):
        """Test the return units"""

        assert n_c.unit.is_equivalent(u.m**-3)

    def test_value(self, n_c):
        """
        Compare the calculated critical density with the expected value.

        The plasma critical density is defined as the value at which the plasma frequency
        equals the frequency of radiation.
        """

        assert np.isclose(n_c.value, self.n_i.value)
