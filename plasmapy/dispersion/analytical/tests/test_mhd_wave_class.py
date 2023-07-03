import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.analytical.mhd_wave_class import mhd_waves
from plasmapy.formulary.speeds import Alfven_speed, ion_sound_speed

kwargs_wave_limits = {
    "k": 0.01 * u.rad / u.m,
    "theta": [0, 90] * u.deg,
}


class TestMHDWave:
    _kwargs_plasma_cold = {
        "B": 8.3e-9 * u.T,
        "ion": "p+",
        "n_i": 5e6 * u.m**-3,
    }
    _kwargs_plasma_hydro = {
        "B": 0 * u.T,
        "ion": "p+",
        "n_i": 5e6 * u.m**-3,
    }
    _kwargs_plasma_temp = {
        "T_e": 1.6e6 * u.K,
        "T_i": 4.0e5 * u.K,
    }

    _v_A = Alfven_speed(
        B=_kwargs_plasma_cold["B"],
        ion=_kwargs_plasma_cold["ion"],
        density=_kwargs_plasma_cold["n_i"],
    )
    _c_s = ion_sound_speed(
        **_kwargs_plasma_temp, ion=_kwargs_plasma_cold["ion"], gamma_e=1, gamma_i=3
    )
    _c_m = np.sqrt(_v_A**2 + _c_s**2)

    @pytest.mark.parametrize(
        ("kwargs", "expected"),
        [
            (  # Cold plasma limit
                {**_kwargs_plasma_cold},
                {
                    "alfven": [_v_A, 0 * u.m / u.s],
                    "fast": [_v_A, _v_A],
                    "slow": [0, 0] * u.m / u.s,
                },
            ),
            (  # No magnetic field
                {**_kwargs_plasma_hydro, **_kwargs_plasma_temp},
                {
                    "alfven": [0, 0] * u.m / u.s,
                    "fast": [_c_s, _c_s],
                    "slow": [0, 0] * u.m / u.s,
                },
            ),
            (  # Finite B and temperature with plasma beta > 1
                {**_kwargs_plasma_cold, **_kwargs_plasma_temp},
                {
                    "alfven": [_v_A, 0 * u.m / u.s],
                    "fast": [_c_s, _c_m],
                    "slow": [_v_A, 0 * u.m / u.s],
                },
            ),
        ],
    )
    def test_angular_frequency_limiting_vals(self, kwargs, expected):
        waves = mhd_waves(**kwargs)
        for mode in waves:
            omega = waves[mode].angular_frequency(**kwargs_wave_limits)
            assert np.all(np.isclose(omega / kwargs_wave_limits["k"], expected[mode]))

    @pytest.mark.parametrize(
        ("kwargs", "expected"),
        [
            ({"k": 0.01 * u.rad / u.m, "theta": 0 * u.deg}, ()),
            ({"k": [0.01, 0.02, 0.03, 0.04] * u.rad / u.m, "theta": 0 * u.deg}, (4,)),
            ({"k": 0.01 * u.rad / u.m, "theta": [0, 45, 90] * u.deg}, (3,)),
            (
                {
                    "k": [0.01, 0.02, 0.03, 0.04] * u.rad / u.m,
                    "theta": [0, 45, 90] * u.deg,
                },
                (4, 3),
            ),
        ],
    )
    def test_angular_frequency_return_structure(self, kwargs, expected):
        waves = mhd_waves(8.3e-9 * u.T, "p+", 5e6 * u.m**-3)

        assert isinstance(waves, dict)
        assert {"alfven", "fast", "slow"} == set(waves.keys())

        for mode in waves:
            omega = waves[mode].angular_frequency(**kwargs)
            assert isinstance(omega, u.Quantity)
            assert omega.unit == u.rad / u.s
            assert omega.shape == expected
