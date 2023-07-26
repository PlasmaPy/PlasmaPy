"""
Tests for particle_tracker.py
"""
import astropy.units as u
import numpy as np
import pytest

from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

from plasmapy.formulary.lengths import gyroradius
from plasmapy.particles import CustomParticle
from plasmapy.plasma.grids import CartesianGrid
from plasmapy.simulation.particle_tracker import (
    MemoryIntervalSaveRoutine,
    NoFieldsStoppingCondition,
    ParticleTracker,
    TimeElapsedStopCondition,
)


class TestParticleTrackerGyroradius:
    point_particle = CustomParticle(1 * u.kg, 1 * u.C)

    def instantiate_simulation(self, B_strength=1 * u.T):
        L = 1e6 * u.km
        num = 2
        grid = CartesianGrid(-L, L, num=num)
        grid_shape = (num,) * 3

        Bz = np.full(grid_shape, B_strength) * u.T
        grid.add_quantities(B_z=Bz)

        return ParticleTracker(
            grid, req_quantities=["E_x", "E_y", "E_z", "B_x", "B_y", "B_z"]
        )

    @pytest.mark.slow()
    @given(arrays(int, (5,), elements=st.integers(1, 100)), st.integers(1, 100))
    @settings(deadline=5e4, max_examples=3)
    def test_gyroradius(self, v_x, B_strength):
        """Test to ensure particles maintain their gyroradius over time"""
        v_x = v_x * u.m / u.s
        v = np.array([[v_x_element.value, 0, 0] for v_x_element in v_x]) * u.m / u.s

        # Set the initial position to the gyroradius
        # This means the particle will orbit the origin
        R_L = gyroradius(B_strength, self.point_particle, Vperp=v_x)
        x = np.array([[0, R_L_element.value, 0] for R_L_element in R_L]) * u.m

        simulation = self.instantiate_simulation(B_strength=B_strength)
        simulation.load_particles(x, v, self.point_particle)

        stop_condition = TimeElapsedStopCondition(10 * u.s)
        save_routine = MemoryIntervalSaveRoutine(0.1 * u.s)

        simulation.run(stop_condition, save_routine, dt=1e-3 * u.s)

        positions = np.asarray(save_routine.x_all) * u.m
        distances = np.linalg.norm(positions, axis=-1)

        assert np.isclose(distances, R_L, rtol=1e-2).all()

    @pytest.mark.slow()
    @given(arrays(int, (5,), elements=st.integers(1, 100)))
    @settings(deadline=1e5, max_examples=3)
    def test_kinetic_energy(self, v_x):
        """Test to ensure particles maintain their gyroradius over time"""
        v_x = v_x * u.m / u.s
        v = np.array([[v_x_element.value, 0, 0] for v_x_element in v_x]) * u.m / u.s
        initial_kinetic_energies = 0.5 * self.point_particle.mass * v_x**2

        # Set the initial position to the gyroradius
        # This means the particle will orbit the origin
        R_L = gyroradius(self.B_strength, self.point_particle, Vperp=v_x)
        x = np.array([[0, R_L_element.value, 0] for R_L_element in R_L]) * u.m

        simulation = self.instantiate_simulation()
        simulation.load_particles(x, v, self.point_particle)

        stop_condition = TimeElapsedStopCondition(10 * u.s)
        save_routine = MemoryIntervalSaveRoutine(0.1 * u.s)

        simulation.run(stop_condition, save_routine, dt=1e-3 * u.s)

        velocities = np.asarray(save_routine.v_all) * u.m / u.s
        speeds = np.linalg.norm(velocities, axis=-1)
        simulation_kinetic_energies = 0.5 * self.point_particle.mass * speeds**2

        assert np.isclose(initial_kinetic_energies, simulation_kinetic_energies).all()


@pytest.mark.slow()
@given(
    st.integers(1, 100), st.integers(1, 100), st.integers(1, 100), st.integers(1, 100)
)
@settings(deadline=2e4, max_examples=10)
def test_particle_tracker_potential_difference(E_strength, L, mass, charge):
    # Apply appropriate units to the random inputs
    E_strength = E_strength * u.V / u.m
    L = L * u.m
    mass = mass * u.kg
    charge = charge * u.C

    num = 2
    dt = 1e-3 * u.s

    grid = CartesianGrid(-L, L, num=num)
    grid_shape = (num,) * 3

    Ex = np.full(grid_shape, E_strength) * u.V / u.m
    grid.add_quantities(E_x=Ex)

    point_particle = CustomParticle(mass, charge)

    x = [[0, 0, 0]] * u.m
    v = [[0, 0, 0]] * u.m / u.s

    simulation = ParticleTracker(
        grid, req_quantities=["E_x", "E_y", "E_z", "B_x", "B_y", "B_z"]
    )
    simulation.load_particles(x, v, point_particle)

    stop_condition = NoFieldsStoppingCondition()
    save_routine = MemoryIntervalSaveRoutine(dt)

    simulation.run(stop_condition, save_routine, dt=dt)

    velocities = np.asarray(save_routine.v_all)[:, 0] * u.m / u.s
    speeds = np.linalg.norm(velocities, axis=-1)

    # Final energy is given by the product of the charge and potential difference
    final_expected_energy = (E_strength * L * point_particle.charge).to(u.J)
    final_simulated_energy = (0.5 * point_particle.mass * speeds[-1] ** 2).to(u.J)

    assert np.isclose(final_expected_energy, final_simulated_energy, rtol=1e-2)
