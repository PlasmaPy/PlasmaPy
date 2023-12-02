"""
Tests for particle_tracker.py
"""
import astropy.units as u
import numpy as np
import pytest

from hypothesis import given, settings
from hypothesis import strategies as st

from plasmapy.formulary.lengths import gyroradius
from plasmapy.particles import CustomParticle
from plasmapy.plasma import Plasma
from plasmapy.plasma.grids import CartesianGrid
from plasmapy.simulation.particle_tracker import (
    IntervalSaveRoutine,
    NoParticlesOnGridsStoppingCondition,
    ParticleTracker,
    TimeElapsedStopCondition,
)

rng = np.random.default_rng()

# (([CartesianGrid(-1 * u.m, 1 * u.m), CartesianGrid(-10 * u.m, 10 * u.m)]), ),


@pytest.mark.parametrize(
    ("constructor_args", "expected_exception"),
    [
        # Deprecation error
        (
            (
                Plasma(
                    domain_x=np.linspace(-1, 1, 10) * u.m,
                    domain_y=np.linspace(-1, 1, 10) * u.m,
                    domain_z=np.linspace(-1, 1, 10) * u.m,
                ),
            ),
            TypeError,
        ),
        # Unrecognized grid type
        ((42,), TypeError),
    ],
)
def test_particle_tracker_constructor_errors(constructor_args, expected_exception):
    with pytest.raises(expected_exception):
        ParticleTracker(*constructor_args)


@pytest.fixture()
def no_particles_on_grids_instantiated():
    return NoParticlesOnGridsStoppingCondition()


@pytest.fixture()
def time_elapsed_stop_condition_instantiated():
    return TimeElapsedStopCondition(1 * u.s)


@pytest.fixture()
def disk_interval_save_routine_instantiated(tmp_path):
    return IntervalSaveRoutine(1 * u.s, output_directory=tmp_path)


@pytest.fixture()
def memory_interval_save_routine_instantiated():
    return IntervalSaveRoutine(1 * u.s)


@pytest.mark.parametrize(
    ("stop_condition", "save_routine"),
    [
        (
            "no_particles_on_grids_instantiated",
            "memory_interval_save_routine_instantiated",
        ),
        (
            "time_elapsed_stop_condition_instantiated",
            "memory_interval_save_routine_instantiated",
        ),
        (
            "no_particles_on_grids_instantiated",
            "disk_interval_save_routine_instantiated",
        ),
        (
            "time_elapsed_stop_condition_instantiated",
            "disk_interval_save_routine_instantiated",
        ),
    ],
)
def test_interval_save_routine(request, stop_condition, save_routine):
    x = [[0, 0, 0]] * u.m
    v = [[0, 1, 0]] * u.m / u.s
    point_particle = CustomParticle(1 * u.kg, 1 * u.C)

    L = 1 * u.m
    num = 2
    grid = CartesianGrid(-L, L, num=num)
    grid_shape = (num,) * 3

    Ex = np.full(grid_shape, 1) * u.V / u.m
    grid.add_quantities(E_x=Ex)

    simulation = ParticleTracker(grid)
    simulation.load_particles(x, v, point_particle)

    stop_condition = request.getfixturevalue(stop_condition)
    save_routine = request.getfixturevalue(save_routine)

    simulation.run(stop_condition, save_routine)


class TestParticleTrackerGyroradius:
    v_x = rng.integers(1, 10, size=100) * u.m / u.s

    v = np.array([[v_x_element.value, 0, 0] for v_x_element in v_x]) * u.m / u.s

    B_strength = rng.integers(1, 10) * u.T

    point_particle = CustomParticle(1 * u.kg, 1 * u.C)

    # Set the initial position to the gyroradius
    # This means the particle will orbit the origin
    R_L = gyroradius(B_strength, point_particle, Vperp=v_x)
    x = np.array([[0, R_L_element.value, 0] for R_L_element in R_L]) * u.m

    L = 1e6 * u.km
    num = 2
    grid = CartesianGrid(-L, L, num=num)
    grid_shape = (num,) * 3

    Bz = np.full(grid_shape, B_strength) * u.T
    grid.add_quantities(B_z=Bz)

    simulation = ParticleTracker(grid)
    simulation.load_particles(x, v, point_particle)

    stop_condition = TimeElapsedStopCondition(6 * u.s)
    save_routine = IntervalSaveRoutine(0.1 * u.s)

    simulation.run(stop_condition, save_routine, dt=1e-2 * u.s)

    def test_gyroradius(self):
        """Test to ensure particles maintain their gyroradius over time"""
        positions = np.asarray(self.save_routine.r_all) * u.m
        distances = np.linalg.norm(positions, axis=-1)

        assert np.isclose(distances, self.R_L, rtol=5e-2).all()

    def test_kinetic_energy(self):
        """Test to ensure particles maintain their gyroradius over time"""

        initial_kinetic_energies = 0.5 * self.point_particle.mass * self.v_x**2

        velocities = np.asarray(self.save_routine.v_all) * u.m / u.s
        speeds = np.linalg.norm(velocities, axis=-1)
        simulation_kinetic_energies = 0.5 * self.point_particle.mass * speeds**2

        assert np.isclose(initial_kinetic_energies, simulation_kinetic_energies).all()


@given(st.integers(1, 10), st.integers(1, 10), st.integers(1, 10), st.integers(1, 10))
@settings(deadline=2e4, max_examples=10)
def test_particle_tracker_potential_difference(request, E_strength, L, mass, charge):
    # Apply appropriate units to the random inputs
    E_strength = E_strength * u.V / u.m
    L = L * u.m
    mass = mass * u.kg
    charge = charge * u.C

    num = 2
    dt = 1e-2 * u.s

    grid = CartesianGrid(-L, L, num=num)
    grid_shape = (num,) * 3

    Ex = np.full(grid_shape, E_strength) * u.V / u.m
    grid.add_quantities(E_x=Ex)

    point_particle = CustomParticle(mass, charge)

    x = [[0, 0, 0]] * u.m
    v = [[0, 0, 0]] * u.m / u.s

    simulation = ParticleTracker(grid)
    simulation.load_particles(x, v, point_particle)

    stop_condition = request.getfixturevalue("no_particles_on_grids_instantiated")
    save_routine = request.getfixturevalue("memory_interval_save_routine_instantiated")

    simulation.run(stop_condition, save_routine, dt=dt)

    velocities = np.asarray(save_routine.v_all)[:, 0] * u.m / u.s
    speeds = np.linalg.norm(velocities, axis=-1)

    # Final energy is given by the product of the charge and potential difference
    final_expected_energy = (E_strength * L * point_particle.charge).to(u.J)
    final_simulated_energy = (0.5 * point_particle.mass * speeds[-1] ** 2).to(u.J)

    assert np.isclose(
        final_expected_energy, final_simulated_energy, atol=0.5, rtol=5e-2
    )


def test_particle_tracker_stop_particles(request):
    E_strength = 1 * u.V / u.m
    L = 1 * u.m
    mass = 1 * u.kg
    charge = 1 * u.C

    num = 2
    dt = 1e-2 * u.s

    grid = CartesianGrid(-L, L, num=num)
    grid_shape = (num,) * 3

    Ex = np.full(grid_shape, E_strength) * u.V / u.m
    grid.add_quantities(E_x=Ex)

    point_particle = CustomParticle(mass, charge)

    x = [[0, 0, 0]] * u.m
    v = [[0, 0, 0]] * u.m / u.s

    stop_condition = request.getfixturevalue("no_particles_on_grids_instantiated")
    save_routine = request.getfixturevalue("memory_interval_save_routine_instantiated")

    simulation = ParticleTracker(grid)
    simulation.load_particles(x, v, point_particle)

    simulation.run(stop_condition, save_routine, dt=dt)
    simulation._stop_particles([True])

    assert np.isnan(simulation.v[0, :]).all()
    assert not np.isnan(simulation.x[0, :]).all()

    simulation._remove_particles([True])
    assert np.isnan(simulation.x[0, :]).all()
