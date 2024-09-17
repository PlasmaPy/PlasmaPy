import astropy.units as u
import numpy as np
import pytest

from plasmapy.dispersion.analytical.mhd_waves_ import mhd_waves
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
    def test_raises_init(self, kwargs, error) -> None:
        """Test scenarios that raise an exception."""
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
            ({"k": 1e-5 * u.eV, "theta": 45 * u.deg}, u.UnitTypeError),
            ({"k": 1e-5 * u.rad / u.m, "theta": np.ones((3, 2)) * u.deg}, ValueError),
            ({"k": 1e-5 * u.rad / u.m, "theta": 5 * u.eV}, u.UnitTypeError),
        ],
    )
    @pytest.mark.parametrize("mode", range(3))
    def test_raises_angular_frequency(self, kwargs, error, mode) -> None:
        """Test scenarios that raise an exception."""
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
    @pytest.mark.parametrize("mode", range(3))
    def test_warns(self, kwargs_wave, error, mode) -> None:
        """Test scenarios the issue a `Warning`."""
        with pytest.warns(error):
            sample_waves[mode].angular_frequency(**kwargs_wave)

    @pytest.mark.parametrize("B", [1e-3, 1e-2])
    @pytest.mark.parametrize("density", [1e16, 1e-11 * u.kg])
    @pytest.mark.parametrize("T", [0, 1e5, 1e6])
    def test_angular_frequency_limiting_vals(self, B, density, T) -> None:
        """Test limiting values of the angular frequencies and phase velocities"""
        waves = mhd_waves(B * u.T, density * u.m**-3, "p+", T=T * u.K)
        v_a = waves[0].alfven_speed
        c_s = waves[0].sound_speed
        c_ms = waves[0].magnetosonic_speed
        expected = [
            [v_a, 0 * u.m / u.s],
            [max(v_a, c_s), c_ms],
            [min(v_a, c_s), 0 * u.m / u.s],
        ]

        for mode in range(3):
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
    def test_angular_frequency_return_structure(self, kwargs, expected) -> None:
        assert isinstance(sample_waves, tuple)

        for mode in range(3):
            omega = sample_waves[mode].angular_frequency(**kwargs)
            assert isinstance(omega, u.Quantity)
            assert omega.unit == u.rad / u.s
            assert omega.shape == expected

    @pytest.mark.parametrize("B", [1e-3, 1e-2])
    @pytest.mark.parametrize("density", [1e16, 1e-11 * u.kg])
    @pytest.mark.parametrize("T", [0, 1e5, 1e6])
    def test_group_velocity_vals(self, B, density, T) -> None:
        """Test limiting values of the group velocities"""
        waves = mhd_waves(B * u.T, density * u.m**-3, "p+", T=T * u.K)

        k = 1e-5 * u.rad / u.m
        theta = np.linspace(0, 2 * np.pi, 10) * u.rad
        dt = 1e-3 * u.rad

        for mode in range(3):
            waves[mode].phase_velocity(k, theta)
            group_velocity_k, group_velocity_theta = waves[mode].group_velocity(
                k, theta
            )

            phase_velocity = waves[mode].phase_velocity(k, theta)
            phase_velocity_p = waves[mode].phase_velocity(k, theta + dt)
            phase_velocity_m = waves[mode].phase_velocity(k, theta - dt)
            # symmetric difference quotient
            dv_dtheta = (phase_velocity_p - phase_velocity_m) / (2 * dt / u.rad)

            assert np.allclose(group_velocity_k, phase_velocity)
            assert np.allclose(group_velocity_theta, dv_dtheta)
