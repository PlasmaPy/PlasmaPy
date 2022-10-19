"""Tests for functionality contained in `plasmapy.formulary.frequencies`."""
import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.densities import critical_density
from plasmapy.formulary.frequencies import plasma_frequency


class TestCriticalDensity:
    """Test the plasma_critical_density function in densities.py."""

    n_i = 5e19 * u.m**-3
    omega = plasma_frequency(n=n_i, particle="e-")

    @pytest.fixture()
    def n_c(self):
        """Get the critical density for the example frequency"""
        return critical_density(self.omega)

    def test_units(self, n_c):
        """Test the return units"""

        assert n_c.unit.is_equivalent(u.m**-3)

    def test_value(self, n_c):
        """
        Compare the calculated critical density with the expected value.

        The plasma critical density for a given frequency of radiation is
        defined as the value at which the electron plasma frequency equals
        the frequency of the radiation.
        """

        assert np.isclose(n_c.value, self.n_i.value)
