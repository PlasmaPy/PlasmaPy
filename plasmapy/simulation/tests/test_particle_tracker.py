"""
Tests for particle_tracker.py
"""
import astropy.units as u
import numpy as np

from plasmapy.formulary.lengths import gyroradius
from plasmapy.particles import CustomParticle
from plasmapy.plasma.grids import CartesianGrid
from plasmapy.simulation.particle_tracker import (
    MemoryIntervalSaveRoutine,
    ParticleTracker,
    TimeElapsedStopCondition,
)


class TestParticleTracker:
    def test_gyroradius(self):
        """Test to ensure particles maintain their gyroradius over time"""
        L = 1 * u.km
        num = 2

        grid = CartesianGrid(-L, L, num=num)

        grid_shape = (num,) * 3
        B_strength = 1 * u.T
        Bz = np.full(grid_shape, B_strength) * u.T

        grid.add_quantities(B_z=Bz)

        simulation = ParticleTracker(
            grid,
            req_quantities=["E_x", "E_y", "E_z", "B_x", "B_y", "B_z"],
            verbose=False,
        )

        x = np.array([[1, 0, 0], [2, 0, 0], [3, 0, 0]]) * u.m
        v = np.array([[0, -1, 0], [0, -2, 0], [0, -3, 0]]) * u.m / u.s
        point_particle = CustomParticle(1 * u.kg, 1 * u.C)

        simulation.load_particles(x, v, point_particle)

        stop_condition = TimeElapsedStopCondition(10 * u.s)
        save_routine = MemoryIntervalSaveRoutine(0.1 * u.s)

        simulation.run(stop_condition, save_routine, dt=1e-3 * u.s)

        positions = np.asarray(save_routine.x_all) * u.m
        distances = np.linalg.norm(positions, axis=-1)

        gyroradii = gyroradius(B_strength, point_particle, Vperp=v[:, 1])

        assert np.isclose(distances, gyroradii, rtol=1e-3).all()
