"""Tests for functions that calculate plasma parameters."""

import numpy as np
import pytest
import warnings

from astropy import units as u
from astropy.constants import m_e, m_p
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary.parameters import (
    Alfven_speed,
    betaH_,
    Bohm_diffusion,
    cs_,
    cwp_,
    DB_,
    Debye_length,
    Debye_number,
    gyrofrequency,
    gyroradius,
    Hall_parameter,
    inertial_length,
    ion_sound_speed,
    lambdaD_,
    lower_hybrid_frequency,
    magnetic_energy_density,
    magnetic_pressure,
    mass_density,
    nD_,
    oc_,
    plasma_frequency,
    pmag_,
    pth_,
    rc_,
    rho_,
    rhoc_,
    thermal_pressure,
    thermal_speed,
    ub_,
    upper_hybrid_frequency,
    va_,
    wc_,
    wlh_,
    wuh_,
)
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils.exceptions import (
    PhysicsError,
    PhysicsWarning,
    PlasmaPyFutureWarning,
    RelativityError,
    RelativityWarning,
)
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray

B = 1.0 * u.T
Z = 1
ion = "p"
m_i = m_p
n_i = 5e19 * u.m ** -3
n_e = Z * 5e19 * u.m ** -3
rho = n_i * m_i + n_e * m_e
T_e = 1e6 * u.K
T_i = 1e6 * u.K
k_1 = 3e1 * u.m ** -1
k_2 = 3e7 * u.m ** -1

B_arr = np.array([0.001, 0.002]) * u.T
B_nanarr = np.array([0.001, np.nan]) * u.T
B_allnanarr = np.array([np.nan, np.nan]) * u.T

rho_arr = np.array([5e-10, 2e-10]) * u.kg / u.m ** 3
rho_infarr = np.array([np.inf, 5e19]) * u.m ** -3
rho_negarr = np.array([-5e19, 6e19]) * u.m ** -3

T_arr = np.array([1e6, 2e6]) * u.K
T_nanarr = np.array([1e6, np.nan]) * u.K
T_nanarr2 = np.array([np.nan, 2e6]) * u.K
T_allnanarr = np.array([np.nan, np.nan]) * u.K
T_negarr = np.array([1e6, -5151.0]) * u.K

V = 25.2 * u.m / u.s
V_arr = np.array([25, 50]) * u.m / u.s
V_nanarr = np.array([25, np.nan]) * u.m / u.s
V_allnanarr = np.array([np.nan, np.nan]) * u.m / u.s

mu = m_p.to(u.u).value


class Test_mass_density:
    r"""Test the mass_density function in parameters.py."""

    @pytest.mark.parametrize(
        "args, kwargs, conditional",
        [
            ((-1 * u.kg * u.m ** -3, "He"), {}, pytest.raises(ValueError)),
            ((-1 * u.m ** -3, "He"), {}, pytest.raises(ValueError)),
            (("not a Quantity", "He"), {}, pytest.raises(TypeError)),
            ((1 * u.m ** -3,), {}, pytest.raises(TypeError)),
            ((1 * u.J, "He"), {}, pytest.raises(u.UnitTypeError)),
            ((1 * u.m ** -3, None), {}, pytest.raises(TypeError)),
            (
                (1 * u.m ** -3, "He"),
                {"z_ratio": "not a ratio"},
                pytest.raises(TypeError),
            ),
        ],
    )
    def test_raises(self, args, kwargs, conditional):
        with conditional:
            mass_density(*args, **kwargs)

    @pytest.mark.parametrize(
        "args, kwargs, expected",
        [
            ((1.0 * u.g * u.m ** -3, ""), {}, 1.0e-3 * u.kg * u.m ** -3),
            ((5.0e12 * u.cm ** -3, "He"), {}, 3.32323849e-8 * u.kg * u.m ** -3),
            (
                (5.0e12 * u.cm ** -3, Particle("He")),
                {},
                3.32323849e-8 * u.kg * u.m ** -3,
            ),
            (
                (5.0e12 * u.cm ** -3, "He"),
                {"z_ratio": 0.5},
                1.66161925e-08 * u.kg * u.m ** -3,
            ),
            (
                (5.0e12 * u.cm ** -3, "He"),
                {"z_ratio": -0.5},
                1.66161925e-08 * u.kg * u.m ** -3,
            ),
        ],
    )
    def test_values(self, args, kwargs, expected):
        assert np.isclose(mass_density(*args, **kwargs), expected)

    def test_handle_nparrays(self):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(mass_density)


# Assertions below that are in CGS units with 2-3 significant digits
# are generally from the NRL Plasma Formulary.


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


class Test_Ion_Sound_Speed:
    r"""Test the ion_sound_speed function in parameters.py."""

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
                    "k": 0 * u.m ** -1,
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
                    "k": 0 * u.m ** -1,
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
                    "n_e": np.array([5, 6, 5]) * u.m ** -3,
                    "k": np.array([3, 4]) * u.m ** -3,
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


def test_thermal_pressure():
    assert thermal_pressure(T_e, n_i).unit.is_equivalent(u.Pa)

    # TODO: may be array issues with arg "mass"
    assert_can_handle_nparray(thermal_pressure)


def test_gyrofrequency():
    r"""Test the gyrofrequency function in parameters.py."""

    assert gyrofrequency(B, "e-").unit.is_equivalent(u.rad / u.s)

    assert gyrofrequency(B, "e-", to_hz=True).unit.is_equivalent(u.Hz)

    assert np.isclose(gyrofrequency(1 * u.T, "e-").value, 175882008784.72018)

    assert np.isclose(gyrofrequency(2.4 * u.T, "e-").value, 422116821083.3284)

    assert np.isclose(
        gyrofrequency(1 * u.T, "e-", to_hz=True).value, 27992490076.528206
    )

    assert np.isclose(
        gyrofrequency(2.4 * u.T, "e-", signed=True).value, -422116821083.3284
    )

    assert np.isclose(gyrofrequency(1 * u.G, "e-").cgs.value, 1.76e7, rtol=1e-3)

    with pytest.raises(TypeError):
        with pytest.warns(u.UnitsWarning):
            gyrofrequency(u.m, "e-")

    with pytest.raises(u.UnitTypeError):
        gyrofrequency(u.m * 1, "e-")

    assert np.isnan(gyrofrequency(B_nanarr, "e-")[-1])

    # The following is a test to check that equivalencies from astropy
    # are working.
    omega_ce = gyrofrequency(2.2 * u.T, "e-")
    f_ce = (omega_ce / (2 * np.pi)) / u.rad
    f_ce_use_equiv = omega_ce.to(u.Hz, equivalencies=[(u.cy / u.s, u.Hz)])
    assert np.isclose(f_ce.value, f_ce_use_equiv.value)

    with pytest.warns(u.UnitsWarning):
        assert gyrofrequency(5.0, "e-") == gyrofrequency(5.0 * u.T, "e-")

    assert gyrofrequency(B, particle=ion).unit.is_equivalent(u.rad / u.s)

    assert np.isclose(gyrofrequency(1 * u.T, particle="p").value, 95788335.834874)

    assert np.isclose(gyrofrequency(2.4 * u.T, particle="p").value, 229892006.00369796)

    assert np.isclose(gyrofrequency(1 * u.G, particle="p").cgs.value, 9.58e3, rtol=2e-3)

    assert gyrofrequency(-5 * u.T, "p") == gyrofrequency(5 * u.T, "p")

    # Case when Z=1 is assumed
    # assert gyrofrequency(B, particle='p+') == gyrofrequency(B, particle='H-1')

    assert gyrofrequency(B, particle="e+") == gyrofrequency(B, "e-")

    with pytest.warns(u.UnitsWarning):
        gyrofrequency(8, "p")

    with pytest.raises(u.UnitTypeError):
        gyrofrequency(5 * u.m, "p")

    with pytest.raises(InvalidParticleError):
        gyrofrequency(8 * u.T, particle="asdfasd")

    with pytest.warns(u.UnitsWarning):
        # TODO this should be WARNS, not RAISES. and it's probably still raised
        assert gyrofrequency(5.0, "p") == gyrofrequency(5.0 * u.T, "p")

    gyrofrequency(1 * u.T, particle="p")
    # testing for user input Z
    testMeth1 = gyrofrequency(1 * u.T, particle="p", Z=0.8).si.value
    testTrue1 = 76630665.79318453
    errStr = f"gyrofrequency() gave {testMeth1}, should be {testTrue1}."
    assert np.isclose(testMeth1, testTrue1, atol=0.0, rtol=1e-5), errStr

    assert_can_handle_nparray(gyrofrequency, kwargs={"signed": True})

    assert_can_handle_nparray(gyrofrequency, kwargs={"signed": False})


class TestGyroradius:
    """Tests for `plasmapy.formulary.parameters.gyroradius`."""

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


def test_Debye_length():
    r"""Test the Debye_length function in parameters.py."""

    assert Debye_length(T_e, n_e).unit.is_equivalent(u.m)

    assert np.isclose(Debye_length(1 * u.eV, 1 * u.cm ** -3).value, 7.43, atol=0.005)

    with pytest.warns(u.UnitsWarning):
        Debye_length(5, 5 * u.m ** -3)

    with pytest.raises(u.UnitTypeError):
        Debye_length(56 * u.kg, 5 * u.m ** -3)

    with pytest.raises(ValueError):
        Debye_length(5 * u.eV, -5 * u.m ** -3)

    with pytest.raises(ValueError):
        Debye_length(-45 * u.K, 5 * u.m ** -3)

    Tarr2 = np.array([1, 2]) * u.K
    narr3 = np.array([1, 2, 3]) * u.m ** -3
    with pytest.raises(ValueError):
        Debye_length(Tarr2, narr3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_length(2.0, 2.0) == Debye_length(2.0 * u.K, 2.0 * u.m ** -3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_length(2.0 * u.K, 2.0) == Debye_length(2.0, 2.0 * u.m ** -3)

    assert_can_handle_nparray(Debye_length)


def test_Debye_number():
    r"""Test the Debye_number function in parameters.py."""

    assert Debye_number(T_e, n_e).unit.is_equivalent(u.dimensionless_unscaled)

    T_e_eV = T_e.to(u.eV, equivalencies=u.temperature_energy())
    assert np.isclose(Debye_number(T_e, n_e).value, Debye_number(T_e_eV, n_e).value)

    assert np.isclose(Debye_number(1 * u.eV, 1 * u.cm ** -3).value, 1720862385.43342)

    with pytest.warns(u.UnitsWarning):
        Debye_number(T_e, 4)

    with pytest.raises(ValueError):
        Debye_number(None, n_e)

    with pytest.raises(u.UnitTypeError):
        Debye_number(5 * u.m, 5 * u.m ** -3)

    with pytest.raises(u.UnitTypeError):
        Debye_number(5 * u.K, 5 * u.m ** 3)

    with pytest.raises(ValueError):
        Debye_number(5j * u.K, 5 * u.cm ** -3)

    Tarr2 = np.array([1, 2]) * u.K
    narr3 = np.array([1, 2, 3]) * u.m ** -3
    with pytest.raises(ValueError):
        Debye_number(Tarr2, narr3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_number(1.1, 1.1) == Debye_number(1.1 * u.K, 1.1 * u.m ** -3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_number(1.1 * u.K, 1.1) == Debye_number(1.1, 1.1 * u.m ** -3)

    assert_can_handle_nparray(Debye_number)


def test_inertial_length():
    r"""Test the inertial_length function in parameters.py."""

    assert inertial_length(n_i, particle="p").unit.is_equivalent(u.m)

    assert np.isclose(
        inertial_length(mu * u.cm ** -3, particle="p").cgs.value, 2.28e7, rtol=0.01
    )

    inertial_length_electron_plus = inertial_length(5.351 * u.m ** -3, particle="e+")
    assert inertial_length_electron_plus == inertial_length(
        5.351 * u.m ** -3, particle="e"
    )

    assert inertial_length(n_i, particle="p") == inertial_length(n_i, particle="p")

    with pytest.warns(u.UnitsWarning):
        inertial_length(4, particle="p")

    with pytest.raises(u.UnitTypeError):
        inertial_length(4 * u.m ** -2, particle="p")

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m ** -3, particle="p")

    with pytest.raises(InvalidParticleError):
        inertial_length(n_i, particle=-135)

    with pytest.warns(u.UnitsWarning):
        inertial_length_no_units = inertial_length(1e19, particle="p")
        assert inertial_length_no_units == inertial_length(
            1e19 * u.m ** -3, particle="p"
        )

    assert inertial_length(n_e, "e-").unit.is_equivalent(u.m)

    assert np.isclose(
        inertial_length(1 * u.cm ** -3, "e-").cgs.value, 5.31e5, rtol=1e-3
    )

    with pytest.warns(u.UnitsWarning):
        inertial_length(5, "e-")

    with pytest.raises(u.UnitTypeError):
        inertial_length(5 * u.m, "e-")

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m ** -3, "e-")

    with pytest.warns(u.UnitsWarning):
        assert inertial_length(1e19, "e-") == inertial_length(1e19 * u.m ** -3, "e-")

    assert_can_handle_nparray(inertial_length)


def test_magnetic_pressure():
    r"""Test the magnetic_pressure function in parameters.py."""

    assert magnetic_pressure(B_arr).unit.is_equivalent(u.Pa)

    assert magnetic_pressure(B).unit.is_equivalent(u.Pa)

    assert magnetic_pressure(B).unit.name == "Pa"

    assert magnetic_pressure(B).value == magnetic_energy_density(B).value

    assert magnetic_pressure(B) == magnetic_energy_density(B.to(u.G))

    assert np.isclose(magnetic_pressure(B).value, 397887.35772973835)

    with pytest.warns(u.UnitsWarning):
        magnetic_pressure(5)

    with pytest.raises(u.UnitTypeError):
        magnetic_pressure(5 * u.m)

    assert np.isnan(magnetic_pressure(np.nan * u.T))

    with pytest.raises(ValueError):
        magnetic_pressure(5j * u.T)

    assert np.isnan(magnetic_pressure(B_nanarr)[-1])

    with pytest.warns(u.UnitsWarning):
        assert magnetic_pressure(22.2) == magnetic_pressure(22.2 * u.T)

    assert_can_handle_nparray(magnetic_pressure)


def test_magnetic_energy_density():
    r"""Test the magnetic_energy_density function in parameters.py."""

    assert magnetic_energy_density(B_arr).unit.is_equivalent(u.J / u.m ** 3)

    assert magnetic_energy_density(B).unit.is_equivalent("J / m3")

    assert magnetic_energy_density(B).value == magnetic_pressure(B).value

    assert_quantity_allclose(
        magnetic_energy_density(2 * B), 4 * magnetic_energy_density(B)
    )

    assert_quantity_allclose(magnetic_energy_density(B).value, 397887.35772973835)

    assert_quantity_allclose(
        magnetic_energy_density(B), magnetic_energy_density(B.to(u.G))
    )

    assert isinstance(magnetic_energy_density(B_arr), u.Quantity)

    with pytest.warns(u.UnitsWarning):
        magnetic_energy_density(5)

    with pytest.raises(u.UnitTypeError):
        magnetic_energy_density(5 * u.m)

    assert np.isnan(magnetic_energy_density(np.nan * u.T))

    with pytest.raises(ValueError):
        magnetic_energy_density(5j * u.T)

    assert np.isnan(magnetic_energy_density(B_nanarr)[-1])

    with pytest.warns(u.UnitsWarning):
        assert magnetic_energy_density(22.2) == magnetic_energy_density(22.2 * u.T)

    assert_can_handle_nparray(magnetic_energy_density)


def test_upper_hybrid_frequency():
    r"""Test the upper_hybrid_frequency function in parameters.py."""

    omega_uh = upper_hybrid_frequency(B, n_e=n_e)
    omega_uh_hz = upper_hybrid_frequency(B, n_e=n_e, to_hz=True)
    omega_ce = gyrofrequency(B, "e-")
    omega_pe = plasma_frequency(n=n_e, particle="e-")
    assert omega_ce.unit.is_equivalent(u.rad / u.s)
    assert omega_pe.unit.is_equivalent(u.rad / u.s)
    assert omega_uh.unit.is_equivalent(u.rad / u.s)
    assert omega_uh_hz.unit.is_equivalent(u.Hz)
    left_hand_side = omega_uh ** 2
    right_hand_side = omega_ce ** 2 + omega_pe ** 2
    assert np.isclose(left_hand_side.value, right_hand_side.value)

    assert np.isclose(omega_uh_hz.value, 69385868857.90918)

    with pytest.raises(ValueError):
        upper_hybrid_frequency(5 * u.T, n_e=-1 * u.m ** -3)

    with pytest.warns(u.UnitsWarning):
        assert upper_hybrid_frequency(1.2, 1.3) == upper_hybrid_frequency(
            1.2 * u.T, 1.3 * u.m ** -3
        )

    with pytest.warns(u.UnitsWarning):
        assert upper_hybrid_frequency(1.4 * u.T, 1.3) == upper_hybrid_frequency(
            1.4, 1.3 * u.m ** -3
        )

    assert_can_handle_nparray(upper_hybrid_frequency)


def test_lower_hybrid_frequency():
    r"""Test the lower_hybrid_frequency function in parameters.py."""

    ion = "He-4 1+"
    omega_ci = gyrofrequency(B, particle=ion)
    omega_pi = plasma_frequency(n=n_i, particle=ion)
    omega_ce = gyrofrequency(B, "e-")
    omega_lh = lower_hybrid_frequency(B, n_i=n_i, ion=ion)
    omega_lh_hz = lower_hybrid_frequency(B, n_i=n_i, ion=ion, to_hz=True)
    assert omega_ci.unit.is_equivalent(u.rad / u.s)
    assert omega_pi.unit.is_equivalent(u.rad / u.s)
    assert omega_ce.unit.is_equivalent(u.rad / u.s)
    assert omega_lh.unit.is_equivalent(u.rad / u.s)
    left_hand_side = omega_lh ** -2
    right_hand_side = (
        1 / (omega_ci ** 2 + omega_pi ** 2) + omega_ci ** -1 * omega_ce ** -1
    )
    assert np.isclose(left_hand_side.value, right_hand_side.value)

    assert np.isclose(omega_lh_hz.value, 299878691.3223296)

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2 * u.T, n_i=5e19 * u.m ** -3, ion="asdfasd")

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2 * u.T, n_i=-5e19 * u.m ** -3, ion="asdfasd")

    with pytest.raises(ValueError):
        lower_hybrid_frequency(np.nan * u.T, n_i=-5e19 * u.m ** -3, ion="asdfasd")

    with pytest.warns(u.UnitsWarning):
        assert lower_hybrid_frequency(1.3, 1e19, "p+") == lower_hybrid_frequency(
            1.3 * u.T, 1e19 * u.m ** -3, "p+"
        )
    assert_can_handle_nparray(lower_hybrid_frequency)


def test_Bohm_diffusion():
    r"""Test Mag_Reynolds in dimensionless.py"""

    T_e = 5000 * u.K
    B = 10 * u.T

    assert (Bohm_diffusion(T_e, B)).unit == u.m ** 2 / u.s

    with pytest.warns(u.UnitsWarning):
        Bohm_diffusion(5000, B)

    with pytest.raises(u.UnitTypeError):
        Bohm_diffusion(2.2 * u.kg, B)


@pytest.mark.parametrize(
    "alias, parent",
    [
        (rho_, mass_density),
        (va_, Alfven_speed),
        (cs_, ion_sound_speed),
        (pth_, thermal_pressure),
        (betaH_, Hall_parameter),
        (oc_, gyrofrequency),
        (wc_, gyrofrequency),
        (rc_, gyroradius),
        (rhoc_, gyroradius),
        (lambdaD_, Debye_length),
        (nD_, Debye_number),
        (cwp_, inertial_length),
        (pmag_, magnetic_pressure),
        (ub_, magnetic_energy_density),
        (wuh_, upper_hybrid_frequency),
        (wlh_, lower_hybrid_frequency),
        (DB_, Bohm_diffusion),
    ],
)
def test_parameters_aliases(alias, parent):
    """Test all aliases defined in parameters.py"""
    assert alias is parent
