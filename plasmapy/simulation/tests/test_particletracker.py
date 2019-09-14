import numpy as np
import pytest
from astropy import units as u
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit

from plasmapy.simulation.particletracker import ParticleTracker
from plasmapy.classes.sources import AnalyticalPlasma

# def test_basic_particletracker_functionality():
#     plasma = uniform_magnetic_field()

#     s = ParticleTracker(plasma=plasma, dt=1e-3 * u.s, nt=1)
#     assert np.isclose(s.kinetic_energy, 0 * u.J, atol=1e-4 * u.J)

#     # this should crash as neither `dt` nor `NT` are not provided
#     with pytest.raises(ValueError):
#         ParticleTracker(plasma=plasma)


def fit_sine_curve(position, t, expected_gyrofrequency, phase=0):
    def sine(t, amplitude, omega, phase, mean):
        return amplitude * np.sin(omega * t + phase) + mean

    mean = position.mean().si.value
    amplitude = 3 * position.std().si.value
    omega = expected_gyrofrequency.si.value
    params, covariances = curve_fit(sine, position.si.value, t.si.value,
                                    p0=(amplitude, omega, phase, mean))
    stds = np.sqrt(np.diag(covariances))
    return params, stds


def test_particle_uniform_magnetic():
    r"""
        Tests the particle stepper for a uniform magnetic field motion.
    """
    def magnetic_field(r):
        return u.Quantity([[0, 0, 1]]*len(r), u.T)
    def electric_field(r):
        return np.zeros_like(r).value * (u.V / u.m)
    test_plasma = AnalyticalPlasma(magnetic_field, electric_field)

    particle_type = 'N-14++'
    s = ParticleTracker(test_plasma, particle_type=particle_type, dt=1e-2 * u.s,
                nt=int(1e2))

    perp_speed = 0.01 * u.m / u.s
    parallel_speed = 1e-5 * u.m / u.s
    mean_B = 1 * u.T
    expected_gyrofrequency = (s.q * mean_B / s.m).to(1 / u.s)
    expected_gyroradius = perp_speed / expected_gyrofrequency
    expected_gyroperiod = 2 * np.pi / expected_gyrofrequency

    dt = expected_gyroperiod / 100

    s = ParticleTracker(test_plasma, particle_type=particle_type, dt=dt, nt=int(1e4))
    s.v[:, 1] = perp_speed

    s.v[:, 2] = parallel_speed
    s.run()

    x = s.position_history[:, 0, 0]
    z = s.position_history[:, 0, 2]

    try:
        params, stds = fit_sine_curve(x, s.t, expected_gyrofrequency)
    except RuntimeError as e:
        print(s)
        raise e
    estimated_gyroradius = np.abs(params[0]) * u.m
    estimated_gyroradius_std = np.abs(stds[0]) * u.m
    estimated_gyrofrequency = np.abs(params[1]) / u.s
    estimated_gyrofrequency_std = np.abs(stds[1]) / u.s

    assert np.isclose(expected_gyroradius, estimated_gyroradius,
                      atol=estimated_gyroradius_std), \
        "Gyroradii don't match!"

    assert np.isclose(expected_gyrofrequency, estimated_gyrofrequency,
                      atol=estimated_gyrofrequency_std), \
        "Gyrofrequencies don't match!"

    p_init = models.Polynomial1D(degree=1)
    fit_p = fitting.LinearLSQFitter()
    p = fit_p(p_init, s.t, z)

    assert np.allclose(z, p(s.t), atol=1e-4 * u.m), \
        "z-velocity doesn't stay constant!"

    s.test_kinetic_energy()
    # s.plot_trajectories()


def test_particle_exb_drift():
    r"""
        Tests the particle stepper for a field with magnetic field in the Z
        direction, electric field in the y direction. This should produce a
        drift in the negative X direction, with the drift velocity

        v_e = ExB / B^2

        which is independent of ion charge.
    """
    def magnetic_field(r):
        return u.Quantity([[0, 0, 1]]*len(r), u.T)
    def electric_field(r):
        return u.Quantity([[0, 1, 0]]*len(r), u.V/u.m)
    test_plasma =  AnalyticalPlasma(magnetic_field, electric_field)

    expected_drift_velocity = -1 * u.m / u.s

    s = ParticleTracker(test_plasma, 'p', 5, dt=1e-10 * u.s, nt=int(5e3))
    s.v[:, 2] += np.random.normal(size=s.N) * u.m / u.s

    s.run()

    p_init = models.Polynomial1D(degree=1)
    for x in s.position_history[:, :, 0].T:
        fit_p = fitting.LinearLSQFitter()
        p = fit_p(p_init, s.t, x)
        fit_velocity = p.parameters[1] * u.m / u.s

        assert np.allclose(x, p(s.t), atol=1e-3 * u.m), \
            "x position doesn't follow linear fit!"

        assert np.isclose(expected_drift_velocity, fit_velocity,
                          atol=1e-3 * u.m / u.s), \
            "x velocity doesn't agree with expected drift velocity!"

    s.test_kinetic_energy()

    # s.plot_trajectories()


# def test_particle_exb_nonuniform_drift():
#     """
#         Tests the particle stepper for a field with magnetic field in the Z
#         direction, electric field in the y direction. This should produce a
#         drift in the negative X direction, with the drift velocity

#         v_e = ExB / B^2

#         which is independent of ion charge.
#     """
#     def magnetic_field(r):
#         return u.Quantity([[0, 0, 1]]*len(r), u.T)
#     def electric_field(r):
#         c1 = 10 * u.V / u.m**3
#         E_val = u.Quantity([[0, 0, 0]]*len(r), u.V/u.m)
#         E_val[:, 1] = c1 * r[:,1]**2 / 2
#         return E_val

#     test_plasma =  AnalyticalPlasma(magnetic_field, electric_field)
#     perp_speed = 1e3 * u.m / u.s
#     parallel_speed = 1e-5 * u.m / u.s


#     s = ParticleTracker(test_plasma, particle_type='p', dt=1e-10 * u.s, nt=int(1e4))
#     s.v[:, 1] = parallel_speed
#     s.v[:, 2] = perp_speed

#     expected_gyrofrequency = (s.q * 1 * u.T / s.m).to(1 / u.s)
#     expected_gyroradius = perp_speed / expected_gyrofrequency
#     s.x[:, 0] = expected_gyroradius / 2

#     # expected_drift_velocity = -(test_plasma.electric_field_strength /
#     #                             test_plasma.magnetic_field_strength).to(u.m/u.s)
#     # expected_drift_velocity += 0.25 * expected_gyroradius**2 * c1 / u.T
#     # expected_drift_velocity = expected_drift_velocity.mean().to(u.m / u.s)

#     s.run()

#     x = s.position_history[:, 0, 0]

#     s.plot_time_trajectories("x")

    # p_init = models.Polynomial1D(degree=1)
    # fit_p = fitting.LinearLSQFitter()
    # p = fit_p(p_init, s.t, x)
    # fit_velocity = p.parameters[1] * u.m / u.s
    #
    # # s.plot_trajectories()

    # assert np.allclose(x, p(s.t), atol=1e-3 * u.m), \
    #     "x position doesn't follow linear fit!"
    #

    # assert np.isclose(expected_drift_velocity, fit_velocity,
    #                   atol=1e-3 * u.m / u.s), \
    #     "x velocity doesn't agree with expected drift velocity!"
    #
    # s.test_kinetic_energy()


# def test_particle_nonuniform_grid():
#     '''
#         Test the particle stepper when the spatial domain dimensions
#         are unequal
#     '''
#     x = np.linspace(0, 1, 10)*u.m
#     y = np.linspace(0, 1, 20)*u.m
#     z = np.linspace(0, 1, 30)*u.m

#     plasma = Plasma3D(x, y, z)

#     ParticleTracker(plasma, 'e', dt=1e-14*u.s, nt=2).run()

if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
