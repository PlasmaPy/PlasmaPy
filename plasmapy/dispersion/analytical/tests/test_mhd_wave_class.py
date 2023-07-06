import astropy.units.core
import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.analytical.mhd_wave_class import AlfvenWave, mhd_waves
from plasmapy.particles.exceptions import InvalidIonError

kwargs_wave_limits = {
    "k": 0.01 * u.rad / u.m,
    "theta": [0, 90] * u.deg,
}


class TestMHDWave:
    kwargs_wave_limits = {
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
    _test_wave = AlfvenWave(**kwargs_wave_limits, T=_T)
    _v_a = _test_wave.alfven_speed
    _c_s = _test_wave.sound_speed
    _c_ms = _test_wave.magnetosonic_speed

    @pytest.mark.parametrize(
        ("kwargs", "error"),
        [
            ({**kwargs_wave_limits, "B": "wrong type"}, TypeError),
            ({**kwargs_wave_limits, "B": [8e-9, 8.5e-9] * u.T}, ValueError),
            ({**kwargs_wave_limits, "B": -1 * u.T}, ValueError),
            ({**kwargs_wave_limits, "B": 5 * u.m}, u.UnitTypeError),
            ({**kwargs_wave_limits, "ion": {"not": "a particle"}}, TypeError),
            ({**kwargs_wave_limits, "ion": "e-"}, InvalidIonError),
            ({**kwargs_wave_limits, "ion": "He", "Z": "wrong type"}, TypeError),
            ({**kwargs_wave_limits, "density": "wrong type"}, TypeError),
            ({**kwargs_wave_limits, "density": [5e6, 6e6] * u.m**-3}, ValueError),
            ({**kwargs_wave_limits, "density": -5e6 * u.m**-3}, ValueError),
            ({**kwargs_wave_limits, "density": 2 * u.s}, u.UnitTypeError),
            ({**kwargs_wave_limits, "T": "wrong type"}, TypeError),
            ({**kwargs_wave_limits, "T": [1.4e6, 1.7e6] * u.K}, ValueError),
            ({**kwargs_wave_limits, "T": -10 * u.eV}, ValueError),
            ({**kwargs_wave_limits, "T": 2 * u.s}, u.UnitTypeError),
            ({**kwargs_wave_limits, "gamma": "wrong type"}, TypeError),
        ],
    )
    def test_raises_init(self, kwargs, error):
        """Test scenarios that raise an `Exception`."""
        with pytest.raises(error):
            AlfvenWave(**kwargs)

    @pytest.mark.parametrize(
        ("kwargs", "error"),
        [
            ({"k": np.ones((3, 2)) * u.rad / u.m, "theta": 45 * u.deg}, ValueError),
            ({"k": 0 * u.rad / u.m, "theta": 45 * u.deg}, ValueError),
            ({"k": -1.0 * u.rad / u.m, "theta": 45 * u.deg}, ValueError),
            ({"k": 0.01 * u.eV, "theta": 45 * u.deg}, astropy.units.core.UnitTypeError),
            ({"k": 0.01 * u.rad / u.m, "theta": np.ones((3, 2)) * u.deg}, ValueError),
            ({"k": 0.01 * u.rad / u.m, "theta": 5 * u.eV}, u.UnitTypeError),
        ],
    )
    def test_raises_angular_frequency(self, kwargs, error):
        """Test scenarios that raise an `Exception`."""
        _test_wave = AlfvenWave(8.3e-9 * u.T, 5e6 * u.m**-3, "p+")
        with pytest.raises(error):
            _test_wave.angular_frequency(**kwargs)

    @pytest.mark.parametrize(
        ("kwargs", "expected"),
        [
            (  # Cold plasma limit
                {**kwargs_wave_limits},
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
                {**kwargs_wave_limits, "T": _T},
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
            phase_velocity = waves[mode].phase_velocity(**kwargs_wave_limits)
            assert np.allclose(
                omega / kwargs_wave_limits["k"], expected[mode], atol=1e-02
            )
            # test phase_velocity
            assert np.allclose(omega / kwargs_wave_limits["k"], phase_velocity)

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
