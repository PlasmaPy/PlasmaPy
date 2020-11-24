"""
Tests for grids.py
"""

import astropy.units as u
import numpy as np
import pytest
import warnings

from plasmapy.plasma import grids as grids


def test_AbstractGrid():
    grid = grids.AbstractGrid(-1 * u.cm, 1 * u.cm)

    array = grid.grid
    units = grid.units

    pts0, pts1, pts2 = grid.pts0, grid.pts1, grid.pts2


def test_AbstractGrid_number_of_positionals():
    with pytest.raises(TypeError):
        grid = grids.AbstractGrid(1 * u.cm, 1 * u.eV, 1)


def test_AbstractGrid_incompatible_units():
    with pytest.raises(ValueError):
        grid = grids.AbstractGrid(1 * u.cm, 1 * u.eV)


def test_AbstractGrid_invalid_lengths_of_arguments():
    with pytest.raises(ValueError):
        grid = grids.AbstractGrid(-1, [2, 3], units=[u.m, u.m])


def test_CartesianGrid():

    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm)

    array = grid.grid
    x_arr, y_arr, z_arr = grid.x_arr, grid.y_arr, grid.z_arr
    radius = grid.distance_from_origin
    x_axis, y_axis, z_axis = grid.x_axis, grid.y_axis, grid.z_axis
    d_x, d_y, d_z = grid.d_x, grid.d_y, grid.d_z
    is_uniform_grid = grid.is_uniform_grid
    shape = grid.shape
    unit = grid.units

    # Grid should be uniform
    assert grid.is_uniform_grid == True

    # Test initializing with a provided grid
    grid2 = grids.CartesianGrid(
        grid.grid[..., 0] * grid.unit0,
        grid.grid[..., 1] * grid.unit1,
        grid.grid[..., 2] * grid.unit2,
    )

    # Units not all consistent
    with pytest.raises(ValueError):
        grid = grids.CartesianGrid(-1, 1, units=[u.m, u.rad, u.rad])

    with pytest.raises(ValueError):
        grid = grids.CartesianGrid(-1, 1, units=[u.m, u.rad, u.m])


def test_interpolators():
    # Test interpolator

    # Test interpolator
    # Create grid
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=100)
    # Add some data to the grid
    data = grid.add_quantity("positions", grid.grid[..., 0] * grid.unit)

    pos = np.array([0.1, -0.3, 0]) * u.cm
    pos2 = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm

    # Interpolate indices
    i = grid.interpolate_indices(pos)[0]

    pout = grid.grid[i[0], i[1], i[2], :] * grid.unit

    # Assert that nearest grid cell was found
    assert np.allclose(pos, pout, atol=0.03)

    # Interpolate grid values using nearest-neighbor interpolator
    pout = grid.nearest_neighbor_interpolator(pos, "pts0")
    assert np.allclose(pos[0], pout, atol=0.03)

    # Interpolate grid values using volume-weighted interpolator
    pout = grid.volume_averaged_interpolator(pos, "pts0")
    assert np.allclose(pos[0], pout, atol=0.03)

    # Test using multiple arguments
    pout, pout2 = grid.nearest_neighbor_interpolator(pos, "pts0", "pts1")

    # Test multiple grid points
    pos = np.array([[0, 0, 0], [0.5, 0.2, 0]]) * u.cm
    i = grid.interpolate_indices(pos)
    pout = grid.nearest_neighbor_interpolator(pos.value, "pts0")
    pout = grid.volume_averaged_interpolator(pos, "pts0")


def test_NonUniformCartesianGrid():
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm)

    # Grid should be non-uniform
    assert grid.is_uniform_grid == False
