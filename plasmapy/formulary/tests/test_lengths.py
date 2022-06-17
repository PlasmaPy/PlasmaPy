"""Tests for functionality contained in `plasmapy.formulary.lengths`."""
import astropy.units as u
import numpy as np
import pytest
import warnings

from astropy.constants import m_p
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary.lengths import (
    cwp_,
    Debye_length,
    gyroradius,
    inertial_length,
    lambdaD_,
    rc_,
    rhoc_,
)
from plasmapy.formulary.speeds import thermal_speed
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils.exceptions import PlasmaPyFutureWarning
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray

Z = 1
n_i = 5e19 * u.m**-3
n_e = Z * 5e19 * u.m**-3

B = 1.0 * u.T
B_arr = np.array([0.001, 0.002]) * u.T

T_e = 1e6 * u.K
T_i = 1e6 * u.K
T_arr = np.array([1e6, 2e6]) * u.K
T_nanarr = np.array([1e6, np.nan]) * u.K

V = 25.2 * u.m / u.s
V_arr = np.array([25, 50]) * u.m / u.s
V_nanarr = np.array([25, np.nan]) * u.m / u.s

mu = m_p.to(u.u).value


def test_Debye_length():
    r"""Test the Debye_length function in lengths.py."""

    assert Debye_length(T_e, n_e).unit.is_equivalent(u.m)

    assert np.isclose(Debye_length(1 * u.eV, 1 * u.cm**-3).value, 7.43, atol=0.005)

    with pytest.warns(u.UnitsWarning):
        Debye_length(5, 5 * u.m**-3)

    with pytest.raises(u.UnitTypeError):
        Debye_length(56 * u.kg, 5 * u.m**-3)

    with pytest.raises(ValueError):
        Debye_length(5 * u.eV, -5 * u.m**-3)

    with pytest.raises(ValueError):
        Debye_length(-45 * u.K, 5 * u.m**-3)

    Tarr2 = np.array([1, 2]) * u.K
    narr3 = np.array([1, 2, 3]) * u.m**-3
    with pytest.raises(ValueError):
        Debye_length(Tarr2, narr3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_length(2.0, 2.0) == Debye_length(2.0 * u.K, 2.0 * u.m**-3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_length(2.0 * u.K, 2.0) == Debye_length(2.0, 2.0 * u.m**-3)

    assert_can_handle_nparray(Debye_length)


class TestGyroradius:
    """Tests for `plasmapy.formulary.lengths.gyroradius`."""

    @pytest.mark.parametrize(
        "args, kwargs, _error",
        [
            ((u.T, "e-"), {}, TypeError),
            ((5 * u.A, "e-"), {"Vperp": 8 * u.m / u.s}, u.UnitTypeError),
            ((5 * u.T, "e-"), {"Vperp": 8 * u.m}, u.UnitTypeError),
            (
                (np.array([5, 6]) * u.T, "e-"),
                {"Vperp": np.array([5, 6, 7]) * u.m / u.s},
                ValueError,
            ),
            ((3.14159 * u.T, "e-"), {"T": -1 * u.K}, ValueError),
            ((1.1 * u.T, "e-"), {"Vperp": 1 * u.m / u.s, "T": 1.2 * u.K}, ValueError),
            ((1.1 * u.T, "e-"), {"Vperp": 1.1 * u.m, "T": 1.2 * u.K}, u.UnitTypeError),
            ((u.T,), {"particle": "p", "Vperp": 8 * u.m / u.s}, TypeError),
            ((B,), {"particle": "p", "T": -1 * u.K}, ValueError),
            (
                (1.1 * u.T,),
                {"particle": "p", "Vperp": 1 * u.m / u.s, "T": 1.2 * u.K},
                ValueError,
            ),
            (
                (1.1 * u.T,),
                {"particle": "p", "Vperp": 1.1 * u.m, "T": 1.2 * u.K},
                u.UnitTypeError,
            ),
            (
                (1.1 * u.T,),
                {"particle": "p", "Vperp": 1.2 * u.m, "T": 1.1 * u.K},
                u.UnitTypeError,
            ),
            ((B_arr, "e-"), {"Vperp": V, "T": T_arr}, ValueError),
            ((B_arr, "e-"), {"Vperp": V_arr, "T": T_i}, ValueError),
            ((B_arr, "e-"), {"Vperp": V, "T": T_nanarr}, ValueError),
            ((B_arr, "e-"), {"Vperp": V_nanarr, "T": T_i}, ValueError),
            ((0.4 * u.T, "e-"), {"T": 5 * u.eV, "T_i": 7 * u.eV}, ValueError),
        ],
    )
    def test_raises(self, args, kwargs, _error):
        """Test scenarios that raise an exception."""

        with warnings.catch_warnings(), pytest.raises(_error):
            # we don't care about warnings for these tests
            warnings.simplefilter("ignore")

            gyroradius(*args, **kwargs)

    @pytest.mark.parametrize(
        "args, kwargs, nan_mask",
        [
            ((np.nan * u.T,), {"particle": "e-", "T": 1 * u.K}, None),
            ((np.nan * u.T,), {"particle": "e-", "Vperp": 1 * u.m / u.s}, None),
            ((1 * u.T,), {"particle": "e-", "T": np.nan * u.K}, None),
            ((1 * u.T,), {"particle": "e-", "Vperp": np.nan * u.m / u.s}, None),
            (
                ([1, 2, np.nan] * u.T,),
                {"particle": "e-", "T": 1 * u.K},
                [False, False, True],
            ),
            (
                ([1, 2, np.nan] * u.T,),
                {"particle": "e-", "T": [np.nan, 1, 2] * u.K},
                [True, False, True],
            ),
            (
                ([1, np.nan, 2] * u.T,),
                {"particle": "e-", "Vperp": [np.nan, 1, 2] * u.m / u.s},
                [True, True, False],
            ),
        ],
    )
    def test_nan_values(self, args, kwargs, nan_mask):
        if nan_mask is None:
            assert np.all(np.isnan(gyroradius(*args, **kwargs)))
        else:
            rc_isnans = np.isnan(gyroradius(*args, **kwargs))
            assert np.all(rc_isnans[nan_mask])
            assert np.all(np.logical_not(rc_isnans[np.logical_not(nan_mask)]))

    @pytest.mark.parametrize(
        "args, kwargs, expected, atol",
        [
            (
                (1 * u.T,),
                {"particle": "e-", "Vperp": 1e6 * u.m / u.s},
                5.6856301e-06 * u.m,
                None,
            ),
            (
                (1 * u.T,),
                {"particle": "p", "Vperp": 1e6 * u.m / u.s},
                0.01043968 * u.m,
                None,
            ),
            (
                (1 * u.T,),
                {"particle": "positron", "T": 1 * u.MK},
                gyroradius(1 * u.T, particle="e-", T=1 * u.MK),
                None,
            ),
            (
                (B,),
                {"particle": "p", "T": T_i},
                gyroradius(B, particle="H+", T=T_i),
                1e-6,
            ),
            (
                (B,),
                {"particle": "p", "Vperp": V},
                gyroradius(B, particle="p", Vperp=-V),
                None,
            ),
            (
                (B_arr,),
                {"particle": "e-", "Vperp": V_arr},
                [1.42140753e-07, 1.42140753e-07] * u.m,
                None,
            ),
            (
                (B_arr,),
                {"particle": "e-", "T": T_arr},
                [0.03130334, 0.02213481] * u.m,
                None,
            ),
            (
                ([0.4, 0.6, 0.8] * u.T,),
                {"particle": "e-", "T": [6, 4, 2] * u.eV},
                [2.06499941e-05, 1.12404331e-05, 5.96113984e-06] * u.m,
                None,
            ),
            ((0.4 * u.T, "p"), {"T": 5800 * u.K}, 0.00025539 * u.m, None),
            (
                (0.4 * u.T, "p"),
                {"T": (5800 * u.K).to(u.eV, equivalencies=u.temperature_energy())},
                0.00025539 * u.m,
                None,
            ),
            #
            # If both Vperp or T are given, but only one of the two is valid
            # at each element, then the valid union is taken
            (
                ([0.001, 0.002] * u.T, "e-"),
                {"Vperp": [25, np.nan] * u.m / u.s, "T": [np.nan, 2e6] * u.K},
                [1.42140753e-07, 2.21348073e-02] * u.m,
                None,
            ),
            #
            # If either Vperp or T is a valid scalar and the other is a Qarray of
            # all nans, then the Qarray of nans should be ignored
            (
                ([0.001, 0.002] * u.T, "e-"),
                {"Vperp": 25.2 * u.m / u.s, "T": [np.nan, np.nan] * u.K},
                [1.43277879e-07, 7.16389393e-08] * u.m,
                None,
            ),
            (
                ([0.001, 0.002] * u.T, "e-"),
                {"Vperp": [np.nan, np.nan] * u.m / u.s, "T": 1e6 * u.K},
                [0.03130334, 0.01565167] * u.m,
                None,
            ),
        ],
    )
    def test_values(self, args, kwargs, expected, atol):
        if atol is None:
            atol = 1e-8

        # note allclose() checks values and units
        rc = gyroradius(*args, **kwargs)
        assert np.allclose(rc, expected, atol=atol)
        assert rc.unit == u.m

    @pytest.mark.parametrize(
        "args, kwargs, expected, _warns",
        [
            ((1.0, "e-"), {"Vperp": 1.0}, 5.6856301e-12 * u.m, u.UnitsWarning),
            ((1.0 * u.T, "e-"), {"Vperp": 1.0}, 5.6856301e-12 * u.m, u.UnitsWarning),
            (
                (1.0, "e-"),
                {"Vperp": 1.0 * u.m / u.s},
                5.6856301e-12 * u.m,
                u.UnitsWarning,
            ),
            ((1.1, "e-"), {"T": 1.2}, 3.11737236e-08 * u.m, u.UnitsWarning),
            ((1.1 * u.T, "e-"), {"T": 1.2}, 3.11737236e-08 * u.m, u.UnitsWarning),
            ((1.1, "e-"), {"T": 1.2 * u.K}, 3.11737236e-08 * u.m, u.UnitsWarning),
            #
            # Future warning for using T_i instead of T
            (
                (1.1 * u.T, "e-"),
                {"T_i": 1.2 * u.K},
                3.11737236e-08 * u.m,
                PlasmaPyFutureWarning,
            ),
        ],
    )
    def test_warns(self, args, kwargs, expected, _warns):
        with pytest.warns(_warns):
            rc = gyroradius(*args, **kwargs)
            if expected is not None:
                assert np.allclose(rc, expected)

    def test_keeps_arguments_unchanged(self):
        Vperp1 = u.Quantity([np.nan, 1], unit=u.m / u.s)
        Vperp2 = Vperp1.copy()
        T = u.Quantity([1, np.nan], unit=u.K)

        gyroradius(B_arr, "e-", Vperp=Vperp1, T=T)

        assert_quantity_allclose(Vperp1, Vperp2)

    def test_correct_thermal_speed_used(self):
        """
        Test the correct version of thermal_speed is used when
        temperature is given.
        """
        B = 123 * u.G
        T = 1.2 * u.MK
        particle = "alpha"

        vperp = thermal_speed(T, particle=particle, method="most_probable", ndim=3)

        assert gyroradius(B, particle=particle, T=T) == gyroradius(
            B, particle=particle, Vperp=vperp
        )


def test_inertial_length():
    r"""Test the inertial_length function in lengths.py."""

    assert inertial_length(n_i, particle="p").unit.is_equivalent(u.m)

    assert np.isclose(
        inertial_length(mu * u.cm**-3, particle="p").cgs.value, 2.28e7, rtol=0.01
    )

    inertial_length_electron_plus = inertial_length(5.351 * u.m**-3, particle="e+")
    assert inertial_length_electron_plus == inertial_length(
        5.351 * u.m**-3, particle="e"
    )

    assert inertial_length(n_i, particle="p") == inertial_length(n_i, particle="p")

    with pytest.warns(u.UnitsWarning):
        inertial_length(4, particle="p")

    with pytest.raises(u.UnitTypeError):
        inertial_length(4 * u.m**-2, particle="p")

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m**-3, particle="p")

    with pytest.raises(InvalidParticleError):
        inertial_length(n_i, particle=-135)

    with pytest.warns(u.UnitsWarning):
        inertial_length_no_units = inertial_length(1e19, particle="p")
        assert inertial_length_no_units == inertial_length(
            1e19 * u.m**-3, particle="p"
        )

    assert inertial_length(n_e, "e-").unit.is_equivalent(u.m)

    assert np.isclose(
        inertial_length(1 * u.cm**-3, "e-").cgs.value, 5.31e5, rtol=1e-3
    )

    with pytest.warns(u.UnitsWarning):
        inertial_length(5, "e-")

    with pytest.raises(u.UnitTypeError):
        inertial_length(5 * u.m, "e-")

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m**-3, "e-")

    with pytest.warns(u.UnitsWarning):
        assert inertial_length(1e19, "e-") == inertial_length(1e19 * u.m**-3, "e-")

    assert_can_handle_nparray(inertial_length)


@pytest.mark.parametrize(
    "alias, parent",
    [
        (cwp_, inertial_length),
        (lambdaD_, Debye_length),
        (rc_, gyroradius),
        (rhoc_, gyroradius),
    ],
)
def test_aliases(alias, parent):
    """Test all aliases defined in lengths.py"""
    assert alias is parent
