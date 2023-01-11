"""Tests for functionality contained in `plasmapy.formulary.lengths`."""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.speeds import Alfven_speed, cs_, ion_sound_speed, va_
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils.exceptions import (
    PhysicsError,
    PhysicsWarning,
    RelativityError,
    RelativityWarning,
)
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray

Z = 1

k_1 = 3e1 * u.m**-1
k_2 = 3e7 * u.m**-1

n_e = Z * 5e19 * u.m**-3

T_e = 1e6 * u.K
T_i = 1e6 * u.K
T_nanarr = np.array([1e6, np.nan]) * u.K
T_negarr = np.array([1e6, -5151.0]) * u.K


@pytest.mark.parametrize(
    "alias, parent",
    [
        (va_, Alfven_speed),
        (cs_, ion_sound_speed),
    ],
)
def test_aliases(alias, parent):
    """Test all aliases defined in speeds.py"""
    assert alias is parent


class TestAlfvenSpeed:
    """Test `~plasmapy.formulary.speeds.Alfven_speed`."""

    @pytest.mark.parametrize(
        "args, kwargs, _error",
        [
            # scenarios that raise RelativityError
            ((10 * u.T, 1.0e-10 * u.kg * u.m**-3), {}, RelativityError),
            ((np.inf * u.T, 1 * u.m**-3), {"ion": "p"}, RelativityError),
            ((-np.inf * u.T, 1 * u.m**-3), {"ion": "p"}, RelativityError),
            #
            # scenarios that raise InvalidParticleError
            ((1 * u.T, 5e19 * u.m**-3), {"ion": "spacecats"}, InvalidParticleError),
            #
            # scenarios that raise TypeError
            (("not a Bfield", 1.0e-10 * u.kg * u.m**-3), {}, TypeError),
            ((10 * u.T, "not a density"), {}, TypeError),
            ((10 * u.T, 5), {"ion": "p"}, TypeError),
            ((1 * u.T, 1.0e18 * u.m**-3), {"ion": ["He"]}, TypeError),
            ((1 * u.T, 1.0e18 * u.m**-3), {"ion": "He", "z_mean": "nope"}, TypeError),
            #
            # scenarios that raise UnitTypeError
            ((1 * u.T, 1.0e18 * u.cm), {"ion": "He"}, u.UnitTypeError),
            ((1 * u.T, 5 * u.m**-2), {"ion": "p"}, u.UnitTypeError),
            ((1 * u.cm, 1.0e18 * u.m**-3), {"ion": "He"}, u.UnitTypeError),
            ((5 * u.A, 5e19 * u.m**-3), {"ion": "p"}, u.UnitTypeError),
            #
            # scenarios that raise ValueError
            ((1 * u.T, -1.0e18 * u.m**-3), {"ion": "He"}, ValueError),
            (
                (np.array([5, 6, 7]) * u.T, np.array([5, 6]) * u.m**-3),
                {"ion": "p"},
                ValueError,
            ),
            (
                (np.array([0.001, 0.002]) * u.T, np.array([-5e19, 6e19]) * u.m**-3),
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
                (5 * u.T, 5e19 * u.m**-3),
                {"ion": "H"},
                15413707.39,
                {},
                RelativityWarning,
            ),
            (
                (5 * u.T, 5e19 * u.m**-3),
                {"ion": "H+"},
                15413707.39,
                {"rtol": 3.0e-4},
                RelativityWarning,
            ),
            (
                (5 * u.T, 5e19 * u.m**-3),
                {"ion": "p"},
                15413707.39,
                {"rtol": 4.0e-4},
                RelativityWarning,
            ),
            #
            # scenarios that issue UnitsWarning
            ((0.5, 1.0e18 * u.m**-3), {"ion": "He"}, 5470657.93, {}, u.UnitsWarning),
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
                (1 * u.T, 1e-8 * u.kg * u.m**-3),
                {"ion": "p"},
                8920620.58 * u.m / u.s,
                {"rtol": 1e-6},
            ),
            (
                (1 * u.T, 1e-8 * u.kg * u.m**-3),
                {},
                8920620.58 * u.m / u.s,
                {"rtol": 1e-6},
            ),
            (
                (0.05 * u.T, 1e18 * u.m**-3),
                {"ion": "He"},
                Alfven_speed(0.05 * u.T, 6.64738793e-09 * u.kg * u.m**-3),
                {},
            ),
            (
                (0.05 * u.T, 1e18 * u.m**-3),
                {"ion": "He+"},
                Alfven_speed(0.05 * u.T, 1e18 * u.m**-3, ion="He"),
                {"rtol": 7e-5},
            ),
            (
                (0.05 * u.T, 1e18 * u.m**-3),
                {"ion": "He", "z_mean": 2},
                Alfven_speed(0.05 * u.T, 1e18 * u.m**-3, ion="He +2"),
                {"rtol": 1.4e-4},
            ),
            (
                (0.05 * u.T, 1e18 * u.m**-3),
                {"ion": Particle("He+")},
                Alfven_speed(0.05 * u.T, 1e18 * u.m**-3, ion="He+"),
                {},
            ),
            (
                ([0.001, 0.002] * u.T, 5e-10 * u.kg * u.m**-3),
                {},
                [
                    va_(0.001 * u.T, 5e-10 * u.kg * u.m**-3).value,
                    va_(0.002 * u.T, 5e-10 * u.kg * u.m**-3).value,
                ]
                * (u.m / u.s),
                {},
            ),
            (
                ([0.001, 0.002] * u.T, [5e-10, 2e-10] * u.kg * u.m**-3),
                {},
                [
                    va_(0.001 * u.T, 5e-10 * u.kg * u.m**-3).value,
                    va_(0.002 * u.T, 2e-10 * u.kg * u.m**-3).value,
                ]
                * (u.m / u.s),
                {},
            ),
            (
                (0.001 * u.T, [1.0e18, 2e18] * u.m**-3),
                {"ion": "p"},
                [
                    va_(0.001 * u.T, 1e18 * u.m**-3, ion="p").value,
                    va_(0.001 * u.T, 2e18 * u.m**-3, ion="p").value,
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
            ((np.nan * u.T, 1 * u.kg * u.m**-3), {}, []),
            ((0.001 * u.T, np.nan * u.kg * u.m**-3), {}, []),
            (([np.nan, 0.001] * u.T, 1 * u.kg * u.m**-3), {}, [True, False]),
            (
                (0.001 * u.T, [np.nan, 1.0, np.nan] * u.kg * u.m**-3),
                {},
                [True, False, True],
            ),
            (([np.nan, 0.001] * u.T, [1, np.nan] * u.kg * u.m**-3), {}, [True, True]),
            (
                (0.001 * u.T, [np.nan, 1e18, np.nan] * u.m**-3),
                {"ion": "Ar+"},
                [True, False, True],
            ),
        ],
    )
    def test_nan_values(self, args, kwargs, nan_mask):
        """Input scenarios that lead to `numpy.nan` values being returned."""
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


class Test_Ion_Sound_Speed:
    r"""Test the ion_sound_speed function in speeds.py."""

    @pytest.mark.parametrize(
        "args, kwargs, expected, isclose_kw",
        [
            (
                (),
                {
                    "T_i": 1.3232 * u.MK,
                    "T_e": 1.831 * u.MK,
                    "ion": "p",
                    "gamma_e": 1,
                    "gamma_i": 3,
                },
                218816.06086407552 * (u.m / u.s),
                {},
            ),
            (
                (1.831 * u.MK, 1.3232 * u.MK, "p"),
                {},
                218816.06086407552 * (u.m / u.s),
                {},
            ),  # Test that function call without keyword argument works correctly
            (
                (),
                {
                    "T_i": 1.3232 * u.MK,
                    "T_e": 1.831 * u.MK,
                    "n_e": n_e,
                    "k": k_1,
                    "ion": "p",
                    "gamma_e": 1,
                    "gamma_i": 3,
                },
                218816.06086407552 * (u.m / u.s),
                {},
            ),
            (
                (),
                {
                    "T_i": 1.3232 * u.MK,
                    "T_e": 1.831 * u.MK,
                    "n_e": n_e,
                    "k": k_2,
                    "ion": "p",
                    "gamma_e": 1,
                    "gamma_i": 3,
                },
                552.3212936293337 * (u.m / u.s),
                {},
            ),
            (
                (),
                {
                    "T_i": 0.88 * u.MK,
                    "T_e": 1.28 * u.MK,
                    "n_e": n_e,
                    "k": 0 * u.m**-1,
                    "ion": "p",
                    "gamma_e": 1.2,
                    "gamma_i": 3.4,
                },
                193328.52857788358 * (u.m / u.s),
                {},
            ),
            (
                (),
                {"T_i": T_i, "T_e": 0 * u.K, "n_e": n_e, "k": k_1, "ion": "p+"},
                ion_sound_speed(T_i=T_i, T_e=0 * u.K, n_e=n_e, k=k_1, ion="p+").value
                * (u.m / u.s),
                {},
            ),
            (
                (),
                {
                    "T_e": 1.2e6 * u.K,
                    "T_i": 0 * u.K,
                    "n_e": n_e,
                    "k": 0 * u.m**-1,
                    "z_mean": 0.8,
                    "ion": "p",
                },
                89018.09 * (u.m / u.s),
                {"atol": 0.0, "rtol": 1e-6},
            ),  # testing for user input z_mean
        ],
    )
    def test_values(self, args, kwargs, expected, isclose_kw):
        assert np.isclose(ion_sound_speed(*args, **kwargs), expected, **isclose_kw)

    # case when Z=1 is assumed
    # assert ion_sound_speed(T_i=T_i, T_e=T_e, ion='p+') == ion_sound_speed(T_i=T_i, T_e=T_e,
    # ion='H-1')

    @pytest.mark.parametrize(
        "kwargs1, kwargs2, _warning",
        [
            ({"T_i": T_i, "T_e": T_e, "n_e": n_e, "ion": "p"}, {}, PhysicsWarning),
            ({"T_i": T_i, "T_e": T_e, "k": k_1, "ion": "p"}, {}, PhysicsWarning),
            ({"T_i": 5e11 * u.K, "T_e": 0 * u.K, "ion": "p"}, {}, RelativityWarning),
            (
                {"T_e": 1.2e6, "T_i": 0 * u.K, "n_e": n_e, "k": k_1, "ion": "p"},
                {"T_e": 1.2e6 * u.K, "T_i": 0 * u.K, "n_e": n_e, "k": k_1, "ion": "p"},
                u.UnitsWarning,
            ),
            (
                {"T_i": 1.3e6, "T_e": 0 * u.K, "n_e": n_e, "k": k_1, "ion": "p"},
                {"T_i": 1.3e6 * u.K, "T_e": 0 * u.K, "n_e": n_e, "k": k_1, "ion": "p"},
                u.UnitsWarning,
            ),
        ],
    )
    def test_warns(self, kwargs1, kwargs2, _warning):
        with pytest.warns(_warning):
            val = ion_sound_speed(**kwargs1)
            if kwargs2 != {}:
                val == ion_sound_speed(**kwargs2)

    @pytest.mark.parametrize(
        "args, kwargs, _error",
        [
            (
                (),
                {
                    "T_i": T_i,
                    "T_e": T_e,
                    "n_e": n_e,
                    "k": k_1,
                    "ion": "p",
                    "gamma_i": np.inf,
                },
                RelativityError,
            ),
            (
                (),
                {
                    "T_i": np.array([5, 6, 5]) * u.K,
                    "T_e": np.array([3, 4]) * u.K,
                    "n_e": np.array([5, 6, 5]) * u.m**-3,
                    "k": np.array([3, 4]) * u.m**-3,
                    "ion": "p",
                },
                u.UnitTypeError,
            ),
            ((5 * u.T), {"ion": "p"}, TypeError),  # Is this test right??????
            ((), {"ion": "p"}, TypeError),
            (
                (),
                {"T_i": T_i, "T_e": 0 * u.K, "gamma_i": 0.9999, "ion": "p"},
                PhysicsError,
            ),
            (
                (),
                {"T_i": T_i, "T_e": 0 * u.K, "gamma_e": 0.9999, "ion": "p"},
                PhysicsError,
            ),
            (
                (),
                {"T_i": T_i, "T_e": 0 * u.K, "gamma_e": "sdjklsf", "ion": "p"},
                TypeError,
            ),
            (
                (),
                {"T_i": T_i, "T_e": 0 * u.K, "gamma_i": "fsdfas", "ion": "p"},
                TypeError,
            ),
            ((), {"T_i": T_i, "T_e": 0 * u.K, "ion": "cupcakes"}, InvalidParticleError),
            ((), {"T_i": -np.abs(T_i), "T_e": 0 * u.K, "ion": "p"}, ValueError),
            (
                (),
                {"T_i": T_i, "T_e": 0 * u.K, "n_e": -np.abs(n_e), "k": k_1, "ion": "p"},
                ValueError,
            ),
            (
                (),
                {"T_i": T_i, "T_e": 0 * u.K, "n_e": n_e, "k": -np.abs(k_1), "ion": "p"},
                ValueError,
            ),
            ((), {"T_i": 5e19 * u.K, "T_e": 0 * u.K, "ion": "p"}, RelativityError),
            (
                (),
                {"T_i": 5 * u.A, "T_e": 0 * u.K, "n_e": n_e, "k": k_1, "ion": "p"},
                u.UnitTypeError,
            ),
            (
                (),
                {"T_i": T_negarr, "T_e": 0 * u.K, "n_e": n_e, "k": k_1, "ion": "p"},
                ValueError,
            ),
            (
                (),
                {"T_e": T_negarr, "T_i": 0 * u.K, "n_e": n_e, "k": k_1, "ion": "p"},
                ValueError,
            ),
        ],
    )
    def test_raises(self, args, kwargs, _error):
        with pytest.raises(_error):
            ion_sound_speed(*args, **kwargs)

    @pytest.mark.parametrize(
        "kwargs",
        [
            ({"T_i": T_nanarr, "T_e": 0 * u.K, "n_e": n_e, "k": k_1, "ion": "p"}),
            ({"T_e": T_nanarr, "T_i": 0 * u.K, "n_e": n_e, "k": k_1, "ion": "p"}),
        ],
    )
    def test_nan_values(self, kwargs):
        np.isnan(ion_sound_speed(**kwargs)[1])

    def test_handle_nparrays(self):
        assert_can_handle_nparray(ion_sound_speed)
