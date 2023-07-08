import astropy.units.core
import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.analytical.mhd_wave_class import mhd_waves
from plasmapy.particles.exceptions import InvalidIonError
from plasmapy.utils.exceptions import PhysicsWarning

kwargs_plasma_cold = {
    "B": 1e-3 * u.T,
    "density": 1e16 * u.m**-3,
    "ion": "p+",
}
kwargs_wave_limits = {
    "k": 1e-5 * u.rad / u.m,
    "theta": [0, 90] * u.deg,
}

sample_waves = mhd_waves(**kwargs_plasma_cold, T=1e6 * u.K)


class TestMHDWave:
    @pytest.mark.parametrize(
        ("kwargs", "error"),
        [
            ({**kwargs_plasma_cold, "B": "wrong type"}, TypeError),
            ({**kwargs_plasma_cold, "B": [8e-9, 8.5e-9] * u.T}, ValueError),
            ({**kwargs_plasma_cold, "B": -1 * u.T}, ValueError),
            ({**kwargs_plasma_cold, "B": 5 * u.m}, u.UnitTypeError),
            ({**kwargs_plasma_cold, "ion": {"not": "a particle"}}, TypeError),
            ({**kwargs_plasma_cold, "ion": "e-"}, InvalidIonError),
            ({**kwargs_plasma_cold, "ion": "He", "Z": "wrong type"}, TypeError),
            ({**kwargs_plasma_cold, "density": "wrong type"}, TypeError),
            ({**kwargs_plasma_cold, "density": [5e6, 6e6] * u.m**-3}, ValueError),
            ({**kwargs_plasma_cold, "density": -5e6 * u.m**-3}, ValueError),
            ({**kwargs_plasma_cold, "density": 2 * u.s}, u.UnitTypeError),
            ({**kwargs_plasma_cold, "T": "wrong type"}, TypeError),
            ({**kwargs_plasma_cold, "T": [1.4e6, 1.7e6] * u.K}, ValueError),
            ({**kwargs_plasma_cold, "T": -10 * u.eV}, ValueError),
            ({**kwargs_plasma_cold, "T": 2 * u.s}, u.UnitTypeError),
            ({**kwargs_plasma_cold, "gamma": "wrong type"}, TypeError),
        ],
    )
    def test_raises_init(self, kwargs, error):
        """Test scenarios that raise an `Exception`."""
        with pytest.raises(error):
            mhd_waves(**kwargs)

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
        for mode in sample_waves:
            with pytest.raises(error):
                sample_waves[mode].angular_frequency(**kwargs)

    @pytest.mark.parametrize(
        ("kwargs_wave", "error"),
        [
            (
                {"k": 1 * u.rad / u.m, "theta": 0 * u.deg},
                PhysicsWarning,
            ),
            (
                {"k": [1e-5, 1] * (u.rad / u.m), "theta": 0 * u.deg},
                PhysicsWarning,
            ),
        ],
    )
    def test_warns(self, kwargs_wave, error):
        """Test scenarios the issue a `Warning`."""
        for mode in sample_waves:
            with pytest.warns(error):
                sample_waves[mode].angular_frequency(**kwargs_wave)

    @pytest.mark.parametrize("B", [1e-3, 1e-2])
    @pytest.mark.parametrize("density", [1e16, 1e-11 * u.kg])
    @pytest.mark.parametrize("T", [0, 1e5, 1e6])
    def test_angular_frequency_limiting_vals(self, B, density, T):
        """Test limiting values of the angular frequencies and phase velocities"""
        waves = mhd_waves(B * u.T, density * u.m**-3, "p+", T=T * u.K)
        v_a = waves["alfven"].alfven_speed
        c_s = waves["alfven"].sound_speed
        c_ms = waves["alfven"].magnetosonic_speed
        expected = {
            "alfven": [v_a, 0 * u.m / u.s],
            "fast": [max(v_a, c_s), c_ms],
            "slow": [min(v_a, c_s), 0 * u.m / u.s],
        }

        for mode in waves:
            omega = waves[mode].angular_frequency(**kwargs_wave_limits)
            v_ph = waves[mode].phase_velocity(**kwargs_wave_limits)
            assert np.allclose(omega / kwargs_wave_limits["k"], expected[mode])
            assert np.allclose(omega / kwargs_wave_limits["k"], v_ph)

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
        assert isinstance(sample_waves, dict)
        assert {"alfven", "fast", "slow"} == set(sample_waves)

        for mode in sample_waves:
            omega = sample_waves[mode].angular_frequency(**kwargs)
            assert isinstance(omega, u.Quantity)
            assert omega.unit == u.rad / u.s
            assert omega.shape == expected
