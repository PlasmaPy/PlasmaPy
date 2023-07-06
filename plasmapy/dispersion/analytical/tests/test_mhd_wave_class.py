import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.analytical.mhd_wave_class import AlfvenWave, mhd_waves

kwargs_wave_limits = {
    "k": 0.01 * u.rad / u.m,
    "theta": [0, 90] * u.deg,
}


class TestMHDWave:
    _kwargs_plasma_cold = {
        "B": 8.3e-9 * u.T,
        "density": 5e6 * u.m**-3,
        "ion": "p+",
    }
    _kwargs_plasma_hydro = {
        "B": 0 * u.T,
        "density": 5e6 * u.m**-3,
        "ion": "p+",
    }
    _T = 1.6e6 * u.K

    # get speeds calculated by an instance
    _test_wave = AlfvenWave(**_kwargs_plasma_cold, T=_T)
    _v_a = _test_wave.alfven_speed
    _c_s = _test_wave.sound_speed
    _c_ms = _test_wave.magnetosonic_speed

    @pytest.mark.parametrize(
        ("kwargs", "expected"),
        [
            (  # Cold plasma limit
                {**_kwargs_plasma_cold},
                {
                    "alfven": [_v_a, 0 * u.m / u.s],
                    "fast": [_v_a, _v_a],
                    "slow": [0, 0] * u.m / u.s,
                },
            ),
            (  # No magnetic field
                {**_kwargs_plasma_hydro, "T": _T},
                {
                    "alfven": [0, 0] * u.m / u.s,
                    "fast": [_c_s, _c_s],
                    "slow": [0, 0] * u.m / u.s,
                },
            ),
            (  # Finite B and temperature with plasma beta > 1
                {**_kwargs_plasma_cold, "T": _T},
                {
                    "alfven": [_v_a, 0 * u.m / u.s],
                    "fast": [_c_s, _c_ms],
                    "slow": [_v_a, 0 * u.m / u.s],
                },
            ),
        ],
    )
    def test_angular_frequency_limiting_vals(self, kwargs, expected):
        waves = mhd_waves(**kwargs)
        for mode in waves:
            omega = waves[mode].angular_frequency(**kwargs_wave_limits)
            assert np.allclose(omega / kwargs_wave_limits["k"], expected[mode], atol=1e-04)

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
        waves = mhd_waves(8.3e-9 * u.T, 5e6 * u.m**-3, "p+")

        assert isinstance(waves, dict)
        assert {"alfven", "fast", "slow"} == set(waves.keys())

        for mode in waves:
            omega = waves[mode].angular_frequency(**kwargs)
            assert isinstance(omega, u.Quantity)
            assert omega.unit == u.rad / u.s
            assert omega.shape == expected
