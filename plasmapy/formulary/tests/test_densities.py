"""Tests for functionality contained in `plasmapy.formulary.frequencies`."""
import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.densities import critical_density, mass_density
from plasmapy.formulary.frequencies import plasma_frequency
from plasmapy.particles import Particle
from plasmapy.utils._pytest_helpers import assert_can_handle_nparray


class TestCriticalDensity:
    """Test the plasma_critical_density function in densities.py."""

    n_i = 5e19 * u.m**-3
    omega = plasma_frequency(n=n_i, particle="e-")

    @pytest.fixture()
    def n_c(self):
        """Get the critical density for the example frequency"""
        return critical_density(self.omega)

    def test_units(self, n_c) -> None:
        """Test the return units"""

        assert n_c.unit.is_equivalent(u.m**-3)

    def test_value(self, n_c) -> None:
        """
        Compare the calculated critical density with the expected value.

        The plasma critical density for a given frequency of radiation is
        defined as the value at which the electron plasma frequency equals
        the frequency of the radiation.
        """

        assert np.isclose(n_c.value, self.n_i.value)


class Test_mass_density:
    r"""Test the mass_density function in densities.py."""

    @pytest.mark.parametrize(
        ("args", "kwargs", "conditional"),
        [
            ((-1 * u.kg * u.m**-3, "He"), {}, pytest.raises(ValueError)),
            ((-1 * u.m**-3, "He"), {}, pytest.raises(ValueError)),
            (("not a Quantity", "He"), {}, pytest.raises(TypeError)),
            ((1 * u.m**-3,), {}, pytest.raises(TypeError)),
            ((1 * u.J, "He"), {}, pytest.raises(u.UnitTypeError)),
            ((1 * u.m**-3, None), {}, pytest.raises(TypeError)),
            (
                (1 * u.m**-3, "He"),
                {"z_ratio": "not a ratio"},
                pytest.raises(TypeError),
            ),
        ],
    )
    def test_raises(self, args, kwargs, conditional) -> None:
        with conditional:
            mass_density(*args, **kwargs)

    @pytest.mark.parametrize(
        ("args", "kwargs", "expected"),
        [
            ((1.0 * u.g * u.m**-3, ""), {}, 1.0e-3 * u.kg * u.m**-3),
            ((5.0e12 * u.cm**-3, "He"), {}, 3.32323849e-8 * u.kg * u.m**-3),
            (
                (5.0e12 * u.cm**-3, Particle("He")),
                {},
                3.32323849e-8 * u.kg * u.m**-3,
            ),
            (
                (5.0e12 * u.cm**-3, "He"),
                {"z_ratio": 0.5},
                1.66161925e-08 * u.kg * u.m**-3,
            ),
            (
                (5.0e12 * u.cm**-3, "He"),
                {"z_ratio": -0.5},
                1.66161925e-08 * u.kg * u.m**-3,
            ),
        ],
    )
    def test_values(self, args, kwargs, expected) -> None:
        assert np.isclose(mass_density(*args, **kwargs), expected)

    def test_handle_nparrays(self) -> None:
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(mass_density)
