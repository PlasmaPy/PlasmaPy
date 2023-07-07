import astropy.units.core
import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.analytical.mhd_wave_class import AlfvenWave, mhd_waves
from plasmapy.particles.exceptions import InvalidIonError
from plasmapy.utils.exceptions import PhysicsWarning

kwargs_wave_limits = {
    "k": 1e-5 * u.rad / u.m,
    "theta": [0, 90] * u.deg,
}


class TestMHDWave:
    _kwargs_plasma_cold = {
        "B": 1e-3 * u.T,
        "density": 1e16 * u.m**-3,
        "ion": "p+",
    }
    _kwargs_plasma_hydro = {
        "B": 0 * u.T,
        "density": 1e16 * u.m**-3,
        "ion": "p+",
    }
    _T = 2.5e6 * u.K

    # get speeds calculated by an instance
    _test_wave = AlfvenWave(**_kwargs_plasma_cold, T=_T)
    _v_a = _test_wave.alfven_speed
    _c_s = _test_wave.sound_speed
    _c_ms = _test_wave.magnetosonic_speed

    @pytest.mark.parametrize(
        ("kwargs", "error"),
        [
            ({**_kwargs_plasma_cold, "B": "wrong type"}, TypeError),
            ({**_kwargs_plasma_cold, "B": [8e-9, 8.5e-9] * u.T}, ValueError),
            ({**_kwargs_plasma_cold, "B": -1 * u.T}, ValueError),
            ({**_kwargs_plasma_cold, "B": 5 * u.m}, u.UnitTypeError),
            ({**_kwargs_plasma_cold, "ion": {"not": "a particle"}}, TypeError),
            ({**_kwargs_plasma_cold, "ion": "e-"}, InvalidIonError),
            ({**_kwargs_plasma_cold, "ion": "He", "Z": "wrong type"}, TypeError),
            ({**_kwargs_plasma_cold, "density": "wrong type"}, TypeError),
            ({**_kwargs_plasma_cold, "density": [5e6, 6e6] * u.m**-3}, ValueError),
            ({**_kwargs_plasma_cold, "density": -5e6 * u.m**-3}, ValueError),
            ({**_kwargs_plasma_cold, "density": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_plasma_cold, "T": "wrong type"}, TypeError),
            ({**_kwargs_plasma_cold, "T": [1.4e6, 1.7e6] * u.K}, ValueError),
            ({**_kwargs_plasma_cold, "T": -10 * u.eV}, ValueError),
            ({**_kwargs_plasma_cold, "T": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_plasma_cold, "gamma": "wrong type"}, TypeError),
        ],
    )
    def test_raises_init(self, kwargs, error):
        """Test scenarios that raise an `Exception`."""
        with pytest.raises(error):
            AlfvenWave(**kwargs)

    @pytest.mark.parametrize(
        ("kwargs", "error"),
        [
            (
                {"k": 1e-5 * np.ones((3, 2)) * u.rad / u.m, "theta": 45 * u.deg},
                ValueError,
            ),
            ({"k": 0 * u.rad / u.m, "theta": 45 * u.deg}, ValueError),
            ({"k": -1.0 * u.rad / u.m, "theta": 45 * u.deg}, ValueError),
            ({"k": 1e-5 * u.eV, "theta": 45 * u.deg}, astropy.units.core.UnitTypeError),
            ({"k": 1e-5 * u.rad / u.m, "theta": np.ones((3, 2)) * u.deg}, ValueError),
            ({"k": 1e-5 * u.rad / u.m, "theta": 5 * u.eV}, u.UnitTypeError),
        ],
    )
    def test_raises_angular_frequency(self, kwargs, error):
        """Test scenarios that raise an `Exception`."""
        _test_wave = AlfvenWave(1e-3 * u.T, 1e16 * u.m**-3, "p+")
        with pytest.raises(error):
            _test_wave.angular_frequency(**kwargs)

    @pytest.mark.parametrize(
        ("kwargs_plasma", "kwargs_wave", "error"),
        [
            (
                {**_kwargs_plasma_cold},
                {"k": 1 * u.rad / u.m, "theta": 0 * u.deg},
                PhysicsWarning,
            ),
            (
                {**_kwargs_plasma_cold},
                {"k": [1e-5, 1] * (u.rad / u.m), "theta": 0 * u.deg},
                PhysicsWarning,
            ),
        ],
    )
    def test_warns(self, kwargs_plasma, kwargs_wave, error):
        """Test scenarios the issue a `Warning`."""
        with pytest.warns(error):
            AlfvenWave(**kwargs_plasma).angular_frequency(**kwargs_wave)

    @pytest.mark.parametrize(
        ("kwargs", "expected"),
        [
            (  # Cold plasma limit
                {**_kwargs_plasma_cold},
                {
                    "alfven": [_v_a, 0 * u.m / u.s],
                    "fast": [_v_a, _v_a],
                    "slow": [0, 0] * (u.m / u.s),
                },
            ),
            (  # No magnetic field
                {**_kwargs_plasma_hydro, "T": _T},
                {
                    "alfven": [0, 0] * (u.m / u.s),
                    "fast": [_c_s, _c_s],
                    "slow": [0, 0] * (u.m / u.s),
                },
            ),
            (  # Finite B and temperature with plasma beta < 1
                {**_kwargs_plasma_cold, "T": _T},
                {
                    "alfven": [_v_a, 0 * u.m / u.s],
                    "fast": [_v_a, _c_ms],
                    "slow": [_c_s, 0 * u.m / u.s],
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
            ({"k": 1e-5 * u.rad / u.m, "theta": 0 * u.deg}, ()),
            ({"k": [1e-5, 2e-5, 3e-5, 4e-5] * (u.rad / u.m), "theta": 0 * u.deg}, (4,)),
            ({"k": 1e-5 * u.rad / u.m, "theta": [0, 45, 90] * u.deg}, (3,)),
            (
                {
                    "k": [1e-5, 2e-5, 3e-5, 4e-5] * (u.rad / u.m),
                    "theta": [0, 45, 90] * u.deg,
                },
                (4, 3),
            ),
        ],
    )
    def test_angular_frequency_return_structure(self, kwargs, expected):
        waves = mhd_waves(1e-3 * u.T, 1e16 * u.m**-3, "p+")

        assert isinstance(waves, dict)
        assert {"alfven", "fast", "slow"} == set(waves.keys())

        for mode in waves:
            omega = waves[mode].angular_frequency(**kwargs)
            assert isinstance(omega, u.Quantity)
            assert omega.unit == u.rad / u.s
            assert omega.shape == expected
