import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit

from plasmapy import formulary
from plasmapy.simulation.particletracker import ParticleTracker, _boris_push
from plasmapy.classes.sources import AnalyticalFields
from plasmapy.utils.exceptions import PhysicsError

def fit_sine_curve(position, t, expected_gyrofrequency, phase=0):
    def sine(t, amplitude, omega, phase, mean):
        return amplitude * np.sin(omega * t + phase) + mean

    mean = position.mean()
    amplitude = 3 * position.std()
    omega = expected_gyrofrequency.si.value
    params, covariances = curve_fit(sine, position, t,
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
    test_plasma = AnalyticalFields(magnetic_field, electric_field)

    s = ParticleTracker(test_plasma)
    sol = s.run(2e-10*u.s, dt=1e-10 * u.s)
    assert np.isfinite(sol.data.position).all()
    assert np.isfinite(sol.data.velocity).all()

def test_adjust_position_velocity():
    def magnetic_field(r):
        return u.Quantity(np.zeros(3), u.T)

    def electric_field(r):
        return u.Quantity(np.zeros(3), E_unit)
    test_plasma = AnalyticalFields(magnetic_field, electric_field)

    s = ParticleTracker(test_plasma)
    s.x = s.x + 2 * u.m
    s.v = s.v + 2 * u.m/u.s

    assert_quantity_allclose(s.x, 2 * u.m)
    assert_quantity_allclose(s.v, 2 * u.m / u.s)

    s.v = u.Quantity(np.zeros_like(s.v), u.m/u.s)
    assert_quantity_allclose(s.kinetic_energy(), 0 * u.J)

    r = repr(s)
    assert 'N = 1' in r
    assert 'AnalyticalFields' in r
    assert 'particle_type=p' in r



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
    _boris_push(x, v, b, e, hqmdt, dt)
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
    _boris_push(x, v, b, e, hqmdt, dt)
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
    _boris_push(x, v, b, e, hqmdt, dt)
    assert np.isfinite(x).all()
    assert np.isfinite(v).all()
    assert ((x - e * dt**2 / 2 - v * dt / 2) == x_copy).all()

def test_particle_uniform_magnetic():
    r"""
        Tests the particle stepper for a uniform magnetic field motion.
    """
    np.random.seed(0)
    def magnetic_field(r):
        return u.Quantity([0, 0, 1], u.T)

    def electric_field(r):
        return u.Quantity(np.zeros(3), E_unit)

    test_plasma = AnalyticalFields(magnetic_field, electric_field)

    particle_type = 'N-14++'
    perp_speed = 0.01 * u.m / u.s
    parallel_speed = 1e-5 * u.m / u.s
    mean_B = 1 * u.T
    expected_gyrofrequency = formulary.gyrofrequency(mean_B, particle_type, to_hz=True)
    expected_gyroradius = formulary.gyroradius(mean_B, particle_type, Vperp = perp_speed)
    expected_gyroperiod = 1 / expected_gyrofrequency

    dt = expected_gyroperiod / 100
    v = u.Quantity([0 * u.m/u.s, perp_speed, parallel_speed]).reshape((1,3))
    s = ParticleTracker(test_plasma, particle_type=particle_type,  v = v)
    sol = s.run(1e4 * dt, dt)

    x = sol.data.position.sel(particle=0, dimension='x')
    z = sol.data.position.sel(particle=0, dimension='z')

    try:
        params, stds = fit_sine_curve(x, sol.data.time, expected_gyrofrequency)
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
    p = fit_p(p_init, sol.data.time, z)

    assert np.allclose(z, p(sol.data.time), atol=1e-4), \
        "z-velocity doesn't stay constant!"

    # s.plot_trajectories()
    sol.test_kinetic_energy()


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
    test_plasma = AnalyticalFields(magnetic_field, electric_field)

    expected_drift_velocity = -1 * u.m / u.s

    v = np.zeros((50, 3))
    v[:, 2] += np.random.normal(size=50)
    s = ParticleTracker(test_plasma, v = v * u.m / u.s)
    assert np.isfinite(s._v).all()
    assert np.isfinite(s._x).all()
    sol = s.run(5e-7 * u.s, dt=1e-10 * u.s)
    assert np.isfinite(sol.data.position).all()
    assert np.isfinite(sol.data.velocity).all()

    p_init = models.Polynomial1D(degree=1)
    for x in sol.data.position.sel(dimension='x').T:
        fit_p = fitting.LinearLSQFitter()
        p = fit_p(p_init, sol.data.time, x)
        fit_velocity = p.parameters[1] * u.m / u.s

        assert np.allclose(x, p(sol.data.time), atol=1e-3), \
            "x position doesn't follow linear fit!"

    assert np.isclose(expected_drift_velocity, fit_velocity,
                      atol=1e-3), \
        "x velocity doesn't agree with expected drift velocity!"

    # s.plot_trajectories()
    with pytest.raises(PhysicsError):   # Kinetic energy is not conserved here due to the electric field
        sol.test_kinetic_energy()

if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-s"])
