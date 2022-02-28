"""Tests for functions that calculate plasma parameters."""

import numpy as np
import pytest

from astropy import units as u
from astropy.constants import m_e, m_p

from plasmapy.formulary import dimensionless, frequencies, lengths, misc, speeds
from plasmapy.formulary.parameters import (
    _grab_charge,
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
from plasmapy.utils.exceptions import PlasmaPyFutureWarning

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
        # misc
        #
        ({"ion": "He+", "z_mean": 1.32}, _grab_charge, misc._grab_charge),
        ({"T_e": 58000 * u.K, "B": 0.4 * u.T}, Bohm_diffusion, misc.Bohm_diffusion),
        ({"T_e": 58000 * u.K, "B": 0.4 * u.T}, DB_, misc.DB_),
        ({"B": 0.4 * u.T}, magnetic_energy_density, misc.magnetic_energy_density),
        ({"B": 0.4 * u.T}, ub_, misc.ub_),
        ({"B": 0.4 * u.T}, magnetic_pressure, misc.magnetic_pressure),
        ({"B": 0.4 * u.T}, pmag_, misc.pmag_),
        (
            {"density": 1e18 * u.m ** -3, "particle": "He+"},
            mass_density,
            misc.mass_density,
        ),
        ({"density": 1e18 * u.m ** -3, "particle": "He+"}, rho_, misc.rho_),
        (
            {"T": 5800 * u.K, "n": 1e18 * u.m ** -3},
            thermal_pressure,
            misc.thermal_pressure,
        ),
        ({"T": 5800 * u.K, "n": 1e18 * u.m ** -3}, pth_, misc.pth_),
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
