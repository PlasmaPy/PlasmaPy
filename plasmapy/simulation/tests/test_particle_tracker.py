"""
Tests for particle_tracker.py
"""
import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.lengths import gyroradius
from plasmapy.particles import CustomParticle
from plasmapy.plasma.grids import CartesianGrid
from plasmapy.simulation.particle_tracker import (
    MemoryIntervalSaveRoutine,
    ParticleTracker,
    TimeElapsedStopCondition,
)

rng = np.random.default_rng()


class TestParticleTracker:
    L = 1e6 * u.km
    num = 2
    grid = CartesianGrid(-L, L, num=num)
    grid_shape = (num,) * 3

    B_strength = 1 * u.T
    Bz = np.full(grid_shape, B_strength) * u.T
    grid.add_quantities(B_z=Bz)

    point_particle = CustomParticle(1 * u.kg, 1 * u.C)

    @pytest.mark.parametrize("v_x", [rng.integers(0, 100, 10) * u.m / u.s])
    def test_gyroradius(self, v_x):
        """Test to ensure particles maintain their gyroradius over time"""
        simulation = ParticleTracker(
            self.grid, req_quantities=["E_x", "E_y", "E_z", "B_x", "B_y", "B_z"]
        )

        v = np.array([[v_x_element.value, 0, 0] for v_x_element in v_x]) * u.m / u.s

        # Set the initial position to the gyroradius
        # This means the particle will orbit the origin
        R_L = gyroradius(self.B_strength, self.point_particle, Vperp=v_x)
        x = np.array([[0, R_L_element.value, 0] for R_L_element in R_L]) * u.m

        simulation.load_particles(x, v, self.point_particle)

        stop_condition = TimeElapsedStopCondition(10 * u.s)
        save_routine = MemoryIntervalSaveRoutine(0.1 * u.s)

        simulation.run(stop_condition, save_routine, dt=1e-3 * u.s)

        positions = np.asarray(save_routine.x_all) * u.m
        distances = np.linalg.norm(positions, axis=-1)

        assert np.isclose(distances, R_L, rtol=1e-3).all()
