"""Tests for functionality contained in `plasmapy.formulary.lengths`."""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.speeds import Alfven_speed, va_
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils.exceptions import RelativityError, RelativityWarning
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray


class TestAlfvenSpeed:
    """Test `~plasmapy.formulary.parameters.Alfven_speed`."""

    @pytest.mark.parametrize("alias", [va_])
    def test_aliases(self, alias):
        assert alias is Alfven_speed

    @pytest.mark.parametrize(
        "args, kwargs, _error",
        [
            # scenarios that raise RelativityError
            ((10 * u.T, 1.0e-10 * u.kg * u.m ** -3), {}, RelativityError),
            ((np.inf * u.T, 1 * u.m ** -3), {"ion": "p"}, RelativityError),
            ((-np.inf * u.T, 1 * u.m ** -3), {"ion": "p"}, RelativityError),
            #
            # scenarios that raise InvalidParticleError
            ((1 * u.T, 5e19 * u.m ** -3), {"ion": "spacecats"}, InvalidParticleError),
            #
            # scenarios that raise TypeError
            (("not a Bfield", 1.0e-10 * u.kg * u.m ** -3), {}, TypeError),
            ((10 * u.T, "not a density"), {}, TypeError),
            ((10 * u.T, 5), {"ion": "p"}, TypeError),
            ((1 * u.T, 1.0e18 * u.m ** -3), {"ion": ["He"]}, TypeError),
            ((1 * u.T, 1.0e18 * u.m ** -3), {"ion": "He", "z_mean": "nope"}, TypeError),
            #
            # scenarios that raise UnitTypeError
            ((1 * u.T, 1.0e18 * u.cm), {"ion": "He"}, u.UnitTypeError),
            ((1 * u.T, 5 * u.m ** -2), {"ion": "p"}, u.UnitTypeError),
            ((1 * u.cm, 1.0e18 * u.m ** -3), {"ion": "He"}, u.UnitTypeError),
            ((5 * u.A, 5e19 * u.m ** -3), {"ion": "p"}, u.UnitTypeError),
            #
            # scenarios that raise ValueError
            ((1 * u.T, -1.0e18 * u.m ** -3), {"ion": "He"}, ValueError),
            (
                (np.array([5, 6, 7]) * u.T, np.array([5, 6]) * u.m ** -3),
                {"ion": "p"},
                ValueError,
            ),
            (
                (np.array([0.001, 0.002]) * u.T, np.array([-5e19, 6e19]) * u.m ** -3),
                {"ion": "p"},
                ValueError,
            ),
        ],
    )
    def test_raises(self, args, kwargs, _error):
        """Test scenarios that raise exceptions or warnings."""
        with pytest.raises(_error):
            Alfven_speed(*args, **kwargs)

    @pytest.mark.parametrize(
        "args, kwargs, expected, isclose_kw, _warning",
        [
            # scenarios that issue RelativityWarning
            (
                (5 * u.T, 5e19 * u.m ** -3),
                {"ion": "H"},
                15413707.39,
                {},
                RelativityWarning,
            ),
            (
                (5 * u.T, 5e19 * u.m ** -3),
                {"ion": "H+"},
                15413707.39,
                {"rtol": 3.0e-4},
                RelativityWarning,
            ),
            (
                (5 * u.T, 5e19 * u.m ** -3),
                {"ion": "p"},
                15413707.39,
                {"rtol": 4.0e-4},
                RelativityWarning,
            ),
            #
            # scenarios that issue UnitsWarning
            ((0.5, 1.0e18 * u.m ** -3), {"ion": "He"}, 5470657.93, {}, u.UnitsWarning),
        ],
    )
    def test_warns(self, args, kwargs, expected, isclose_kw, _warning):
        """Test scenarios that issue warnings"""
        with pytest.warns(_warning):
            val = Alfven_speed(*args, **kwargs)
            assert isinstance(val, u.Quantity)
            assert val.unit == u.m / u.s
            assert np.isclose(val.value, expected, **isclose_kw)

    @pytest.mark.parametrize(
        "args, kwargs, expected, isclose_kw",
        [
            (
                (1 * u.T, 1e-8 * u.kg * u.m ** -3),
                {"ion": "p"},
                8920620.58 * u.m / u.s,
                {"rtol": 1e-6},
            ),
            (
                (1 * u.T, 1e-8 * u.kg * u.m ** -3),
                {},
                8920620.58 * u.m / u.s,
                {"rtol": 1e-6},
            ),
            (
                (0.05 * u.T, 1e18 * u.m ** -3),
                {"ion": "He"},
                Alfven_speed(0.05 * u.T, 6.64738793e-09 * u.kg * u.m ** -3),
                {},
            ),
            (
                (0.05 * u.T, 1e18 * u.m ** -3),
                {"ion": "He+"},
                Alfven_speed(0.05 * u.T, 1e18 * u.m ** -3, ion="He"),
                {"rtol": 7e-5},
            ),
            (
                (0.05 * u.T, 1e18 * u.m ** -3),
                {"ion": "He", "z_mean": 2},
                Alfven_speed(0.05 * u.T, 1e18 * u.m ** -3, ion="He +2"),
                {"rtol": 1.4e-4},
            ),
            (
                (0.05 * u.T, 1e18 * u.m ** -3),
                {"ion": Particle("He+")},
                Alfven_speed(0.05 * u.T, 1e18 * u.m ** -3, ion="He+"),
                {},
            ),
            (
                ([0.001, 0.002] * u.T, 5e-10 * u.kg * u.m ** -3),
                {},
                [
                    va_(0.001 * u.T, 5e-10 * u.kg * u.m ** -3).value,
                    va_(0.002 * u.T, 5e-10 * u.kg * u.m ** -3).value,
                ]
                * (u.m / u.s),
                {},
            ),
            (
                ([0.001, 0.002] * u.T, [5e-10, 2e-10] * u.kg * u.m ** -3),
                {},
                [
                    va_(0.001 * u.T, 5e-10 * u.kg * u.m ** -3).value,
                    va_(0.002 * u.T, 2e-10 * u.kg * u.m ** -3).value,
                ]
                * (u.m / u.s),
                {},
            ),
            (
                (0.001 * u.T, [1.0e18, 2e18] * u.m ** -3),
                {"ion": "p"},
                [
                    va_(0.001 * u.T, 1e18 * u.m ** -3, ion="p").value,
                    va_(0.001 * u.T, 2e18 * u.m ** -3, ion="p").value,
                ]
                * (u.m / u.s),
                {},
            ),
        ],
    )
    def test_values(self, args, kwargs, expected, isclose_kw):
        """Test expected values."""
        assert np.allclose(Alfven_speed(*args, **kwargs), expected, **isclose_kw)

    @pytest.mark.parametrize(
        "args, kwargs, nan_mask",
        [
            ((np.nan * u.T, 1 * u.kg * u.m ** -3), {}, []),
            ((0.001 * u.T, np.nan * u.kg * u.m ** -3), {}, []),
            (([np.nan, 0.001] * u.T, 1 * u.kg * u.m ** -3), {}, [True, False]),
            (
                (0.001 * u.T, [np.nan, 1.0, np.nan] * u.kg * u.m ** -3),
                {},
                [True, False, True],
            ),
            (([np.nan, 0.001] * u.T, [1, np.nan] * u.kg * u.m ** -3), {}, [True, True]),
            (
                (0.001 * u.T, [np.nan, 1e18, np.nan] * u.m ** -3),
                {"ion": "Ar+"},
                [True, False, True],
            ),
        ],
    )
    def test_nan_values(self, args, kwargs, nan_mask):
        """Input scenarios that leat to `numpy.nan` values being returned."""
        val = Alfven_speed(*args, **kwargs)
        if np.isscalar(val.value):
            assert np.isnan(val)
        else:
            nan_arr = np.isnan(val)
            assert np.all(nan_arr[nan_mask])
            assert np.all(np.logical_not(nan_arr[np.logical_not(nan_mask)]))

    def test_handle_nparrays(self):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(Alfven_speed)