import numpy as np
import pytest
from astropy import units as u
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit

from plasmapy.simulation.particletracker import ParticleTracker
from plasmapy.classes.sources import AnalyticalPlasma
from plasmapy.utils.exceptions import PhysicsError

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


# precalculating unit for efficiency
E_unit = u.V / u.m

def test_run_no_fields():
    def magnetic_field(r):
        return u.Quantity(np.zeros(3), u.T)

    def electric_field(r):
        return u.Quantity(np.zeros(3), E_unit)
    test_plasma = AnalyticalPlasma(magnetic_field, electric_field)

    s = ParticleTracker(test_plasma, 'p', 50, dt=1e-10 * u.s, nt=int(2))
    s.run()
    assert np.isfinite(s._x).all()
    assert np.isfinite(s._v).all()

def test_boris_push_no_fields_no_movement():
    x = np.zeros((2, 3), dtype=float)
    x_copy = x.copy()
    v = np.zeros_like(x)
    b = np.zeros_like(x)
    e = np.zeros_like(x)
    q = 1
    m = 1
    dt = 1e-3
    hqmdt = 0.5 * q * m * dt
    ParticleTracker._boris_push(x, v, b, e, hqmdt, dt)
    assert np.isfinite(x).all()
    assert np.isfinite(v).all()
    assert (x == x_copy).all()

def test_boris_push_no_fields_just_velocity():
    x = np.zeros((2, 3), dtype=float)
    x_copy = x.copy()
    v = np.ones_like(x)
    b = np.zeros_like(x)
    e = np.zeros_like(x)
    q = 1
    m = 1
    dt = 1e-3
    hqmdt = 0.5 * q * m * dt
    ParticleTracker._boris_push(x, v, b, e, hqmdt, dt)
    assert np.isfinite(x).all()
    assert np.isfinite(v).all()
    assert ((x - v * dt) == x_copy).all()

def test_boris_push_electric_field():
    x = np.zeros((2, 3), dtype=float)
    x_copy = x.copy()
    v = np.zeros_like(x)
    b = np.zeros_like(x)
    e = np.zeros_like(x)
    e[:,0] = 1

    q = 1
    m = 1
    dt = 1e-3
    hqmdt = 0.5 * q * m * dt
    ParticleTracker._boris_push(x, v, b, e, hqmdt, dt)
    assert np.isfinite(x).all()
    assert np.isfinite(v).all()
    assert ((x - e * dt**2 / 2 - v * dt / 2) == x_copy).all()

@pytest.mark.xfail
def test_set_particle_velocity():
    test_plasma = AnalyticalPlasma(lambda r: None, lambda r: None)
    particle_type = 'N-14++'
    s = ParticleTracker(test_plasma, 'p', dt=1 * u.s, nt=5)
    s.v[0,0] = 5 * u.m/u.s
    assert s._v[0,0] == 5

def test_set_particle_velocity_by_value():
    test_plasma = AnalyticalPlasma(lambda r: np.zeros(3), lambda r: np.zeros(3))
    particle_type = 'N-14++'
    s = ParticleTracker(test_plasma, 'p', dt=1 * u.s, nt=5)
    s.v = np.array([[5, 0, 0]]) * u.m/u.s
    assert s._v[0,0] == 5

def test_particle_uniform_magnetic():
    r"""
        Tests the particle stepper for a uniform magnetic field motion.
    """
    np.random.seed(0)
    def magnetic_field(r):
        return u.Quantity([0, 0, 1], u.T)

    def electric_field(r):
        return u.Quantity(np.zeros(3), E_unit)

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
    s._v[:, 1] = perp_speed.si.value

    s._v[:, 2] = parallel_speed.si.value
    assert s.v[0, 1] == perp_speed
    assert s.v[0, 2] == parallel_speed
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

    assert u.isclose(expected_gyroradius, estimated_gyroradius,
                      atol=estimated_gyroradius_std), \
        "Gyroradii don't match!"

    assert u.isclose(expected_gyrofrequency, estimated_gyrofrequency,
                      atol=estimated_gyrofrequency_std), \
        "Gyrofrequencies don't match!"

    p_init = models.Polynomial1D(degree=1)
    fit_p = fitting.LinearLSQFitter()
    p = fit_p(p_init, s.t, z)

    assert u.allclose(z, p(s.t), atol=1e-4 * u.m), \
        "z-velocity doesn't stay constant!"

    # s.plot_trajectories()
    s.test_kinetic_energy()


def test_particle_exb_drift():
    r"""
        Tests the particle stepper for a field with magnetic field in the Z
        direction, electric field in the y direction. This should produce a
        drift in the negative X direction, with the drift velocity

        v_e = ExB / B^2

        which is independent of ion charge.
    """
    np.random.seed(0)
    def magnetic_field(r):
        return u.Quantity([0, 0, 1], u.T)

    def electric_field(r):
        return u.Quantity([0, 1, 0], E_unit)
    test_plasma = AnalyticalPlasma(magnetic_field, electric_field)

    expected_drift_velocity = -1 * u.m / u.s

    s = ParticleTracker(test_plasma, 'p', 50, dt=1e-10 * u.s, nt=int(5e3))
    s._v[:, 2] += np.random.normal(size=s.N)
    assert np.isfinite(s._v).all()
    assert np.isfinite(s._x).all()
    s.run()
    assert np.isfinite(s.position_history).all()
    assert np.isfinite(s.velocity_history).all()

    p_init = models.Polynomial1D(degree=1)
    for x in s.position_history[:, :, 0].T:
        fit_p = fitting.LinearLSQFitter()
        p = fit_p(p_init, s.t, x)
        fit_velocity = p.parameters[1] * u.m / u.s

        assert u.allclose(x, p(s.t), atol=1e-3 * u.m), \
            "x position doesn't follow linear fit!"

    assert u.isclose(expected_drift_velocity, fit_velocity,
                      atol=1e-3 * u.m / u.s), \
        "x velocity doesn't agree with expected drift velocity!"

    # s.plot_trajectories()
    with pytest.raises(PhysicsError):
        s.test_kinetic_energy()

if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-s"])
