"""Tests for functions that calculate plasma parameters."""

import numpy as np
import pytest

from astropy import units as u
from astropy.constants import m_e, m_p
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary import dimensionless, frequencies, lengths, speeds
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
    kappa_thermal_speed,
    lambdaD_,
    lower_hybrid_frequency,
    magnetic_energy_density,
    magnetic_pressure,
    mass_density,
    nD_,
    oc_,
    plasma_frequency,
    plasma_frequency_lite,
    pmag_,
    pth_,
    rc_,
    rho_,
    rhoc_,
    thermal_pressure,
    thermal_speed,
    thermal_speed_coefficients,
    thermal_speed_lite,
    ub_,
    upper_hybrid_frequency,
    va_,
    vth_,
    vth_kappa_,
    wc_,
    wlh_,
    wp_,
    wuh_,
)
from plasmapy.particles import Particle
from plasmapy.utils.exceptions import PlasmaPyFutureWarning
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


def test_thermal_pressure():
    assert thermal_pressure(T_e, n_i).unit.is_equivalent(u.Pa)

    # TODO: may be array issues with arg "mass"
    assert_can_handle_nparray(thermal_pressure)


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
        (pth_, thermal_pressure),
        (pmag_, magnetic_pressure),
        (ub_, magnetic_energy_density),
        (DB_, Bohm_diffusion),
    ],
)
def test_parameters_aliases(alias, parent):
    """Test all aliases defined in parameters.py"""
    assert alias is parent


@pytest.mark.parametrize(
    "kwargs, deprecated_func, parent",
    [
        # dimensionless
        #
        (
            {"T_e": 58000 * u.K, "n_e": 1e18 * u.m ** -3},
            Debye_number,
            dimensionless.Debye_number,
        ),
        (
            {"T_e": 58000 * u.K, "n_e": 1e18 * u.m ** -3},
            nD_,
            dimensionless.nD_,
        ),
        (
            {
                "n": 1e18 * u.m ** -3,
                "T": 58000 * u.K,
                "B": 0.4 * u.T,
                "ion": "He+",
                "particle": "e-",
            },
            Hall_parameter,
            dimensionless.Hall_parameter,
        ),
        (
            {
                "n": 1e18 * u.m ** -3,
                "T": 58000 * u.K,
                "B": 0.4 * u.T,
                "ion": "He+",
                "particle": "e-",
            },
            betaH_,
            dimensionless.betaH_,
        ),
        #
        # frequencies
        #
        (
            {"B": 0.4 * u.T, "particle": "He+"},
            gyrofrequency,
            frequencies.gyrofrequency,
        ),
        ({"B": 0.4 * u.T, "particle": "He+"}, oc_, frequencies.oc_),
        ({"B": 0.4 * u.T, "particle": "He+"}, wc_, frequencies.wc_),
        (
            {"n": 1.0e18, "mass": 6.64556605e-27, "z_mean": 1.0},
            plasma_frequency_lite,
            frequencies.plasma_frequency_lite,
        ),
        (
            {"n": 1.0e18 * u.m ** -3, "particle": "He+"},
            plasma_frequency,
            frequencies.plasma_frequency,
        ),
        (
            {"n": 1.0e18 * u.m ** -3, "particle": "He+"},
            wp_,
            frequencies.wp_,
        ),
        (
            {"B": 0.4 * u.T, "n_i": 1.0e18 * u.m ** -3, "ion": "He+"},
            lower_hybrid_frequency,
            frequencies.lower_hybrid_frequency,
        ),
        (
            {"B": 0.4 * u.T, "n_i": 1.0e18 * u.m ** -3, "ion": "He+"},
            wlh_,
            frequencies.wlh_,
        ),
        (
            {"B": 0.4 * u.T, "n_e": 1.0e18 * u.m ** -3},
            upper_hybrid_frequency,
            frequencies.upper_hybrid_frequency,
        ),
        (
            {"B": 0.4 * u.T, "n_e": 1.0e18 * u.m ** -3},
            wuh_,
            frequencies.wuh_,
        ),
        #
        # lengths
        #
        (
            {"T_e": 58000 * u.K, "n_e": 1e18 * u.m ** -3},
            Debye_length,
            lengths.Debye_length,
        ),
        (
            {"T_e": 58000 * u.K, "n_e": 1e18 * u.m ** -3},
            lambdaD_,
            lengths.lambdaD_,
        ),
        (
            {"n": 1e18 * u.m ** -3, "particle": "p"},
            inertial_length,
            lengths.inertial_length,
        ),
        (
            {"n": 1e18 * u.m ** -3, "particle": "p"},
            cwp_,
            lengths.cwp_,
        ),
        (
            {"B": 0.4 * u.T, "particle": "He+", "T": 58000 * u.K},
            gyroradius,
            lengths.gyroradius,
        ),
        (
            {"B": 0.4 * u.T, "particle": "He+", "T": 58000 * u.K},
            rc_,
            lengths.rc_,
        ),
        (
            {"B": 0.4 * u.T, "particle": "He+", "T": 58000 * u.K},
            rhoc_,
            lengths.rhoc_,
        ),
        #
        # speeds
        #
        (
            {"B": 0.4 * u.T, "density": 1e18 * u.m ** -3, "ion": "He+"},
            Alfven_speed,
            speeds.Alfven_speed,
        ),
        ({"B": 0.4 * u.T, "density": 1e18 * u.m ** -3, "ion": "He+"}, va_, speeds.va_),
        (
            {"T_e": 58000 * u.K, "T_i": 12000 * u.K, "ion": "He+"},
            ion_sound_speed,
            speeds.ion_sound_speed,
        ),
        ({"T_e": 58000 * u.K, "T_i": 12000 * u.K, "ion": "He+"}, cs_, speeds.cs_),
        (
            {"method": "most_probable", "ndim": 3},
            thermal_speed_coefficients,
            speeds.thermal_speed_coefficients,
        ),
        (
            {"T": 58000, "mass": 6.64556605e-27, "coeff": 1.0},
            thermal_speed_lite,
            speeds.thermal_speed_lite,
        ),
        ({"T": 58000 * u.K, "particle": "He+"}, thermal_speed, speeds.thermal_speed),
        ({"T": 58000 * u.K, "particle": "He+"}, vth_, speeds.vth_),
        (
            {"T": 58000 * u.K, "particle": "He+", "kappa": 4.0},
            kappa_thermal_speed,
            speeds.kappa_thermal_speed,
        ),
        (
            {"T": 58000 * u.K, "particle": "He+", "kappa": 4.0},
            vth_kappa_,
            speeds.vth_kappa_,
        ),
    ],
)
def test_deprecated(kwargs, deprecated_func, parent):
    assert hasattr(deprecated_func, "__wrapped__")
    assert deprecated_func.__wrapped__ is parent

    with pytest.warns(PlasmaPyFutureWarning):
        deprecated_func(**kwargs)
