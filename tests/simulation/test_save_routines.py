"""
Tests for save_routines.py
"""

import warnings

import astropy.units as u
import numpy as np
import pytest

from plasmapy.particles import CustomParticle
from plasmapy.plasma.grids import CartesianGrid
from plasmapy.simulation.particle_tracker.particle_tracker import ParticleTracker
from plasmapy.simulation.particle_tracker.save_routines import IntervalSaveRoutine
from plasmapy.simulation.particle_tracker.termination_conditions import (
    NoParticlesOnGridsTerminationCondition,
    TimeElapsedTerminationCondition,
)


@pytest.fixture
def no_particles_on_grids_instantiated():
    return NoParticlesOnGridsTerminationCondition()


@pytest.fixture
def time_elapsed_termination_condition_instantiated():
    return TimeElapsedTerminationCondition(1 * u.s)


@pytest.fixture
def disk_interval_save_routine_instantiated(tmp_path):
    return IntervalSaveRoutine(
        1 * u.s, output_directory=tmp_path, output_basename="test_name"
    )


@pytest.fixture
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
            "time_elapsed_termination_condition_instantiated",
            "memory_interval_save_routine_instantiated",
        ),
        (
            "no_particles_on_grids_instantiated",
            "disk_interval_save_routine_instantiated",
        ),
        (
            "time_elapsed_termination_condition_instantiated",
            "disk_interval_save_routine_instantiated",
        ),
    ],
)
def test_interval_save_routine(request, stop_condition, save_routine) -> None:
    x = [[0, 0, 0]] * u.m
    v = [[0, 1, 0]] * u.m / u.s
    point_particle = CustomParticle(1 * u.kg, 1 * u.C)

    L = 1 * u.m
    num = 2
    grid = CartesianGrid(-L, L, num=num)
    grid_shape = (num,) * 3

    Ex = np.full(grid_shape, 1) * u.V / u.m
    grid.add_quantities(E_x=Ex)

    termination_condition = request.getfixturevalue(stop_condition)
    save_routine = request.getfixturevalue(save_routine)

    simulation = ParticleTracker(grid, termination_condition, save_routine)
    simulation.load_particles(x, v, point_particle)

    with warnings.catch_warnings():
        # Ignore warning raises by this special small grid
        warnings.filterwarnings("ignore", message="Quantities should go to zero")
        simulation.run()
