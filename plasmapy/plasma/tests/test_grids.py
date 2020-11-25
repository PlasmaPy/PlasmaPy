"""
Tests for grids.py
"""

import astropy.units as u
import numpy as np
import pytest
import warnings

from plasmapy.plasma import grids as grids


def test_AbstractGrid():
    grid = grids.AbstractGrid(-1 * u.cm, 1 * u.cm, num=(10, 20, 5))

    array = grid.grid
    units = grid.units

    pts0, pts1, pts2 = grid.pts0, grid.pts1, grid.pts2

    # Test wrong number of positional arguments: 1 or more than 3
    with pytest.raises(TypeError):
        grid = grids.AbstractGrid(1 * u.cm, num=10)
    with pytest.raises(TypeError):
        grid = grids.AbstractGrid(-1 * u.cm, 1 * u.cm, 1 * u.cm, 1 * u.cm)

    # Test unequal lengths of arguments raises error
    with pytest.raises(ValueError):
        grid = grids.AbstractGrid(-1 * u.m, [2 * u.m, 3 * u.m], num=10)

    with pytest.raises(ValueError):
        grid = grids.AbstractGrid(
            np.random.randn(2, 5, 3) * u.m,
            np.random.randn(2, 5, 3) * u.m,
            np.random.randn(2, 5, 4) * u.m,
        )

    # Test incompatible units
    with pytest.raises(ValueError):
        grid = grids.AbstractGrid(1 * u.cm, 1 * u.eV, num=10)

    # Test adding a quantity
    q1 = np.random.randn(10, 20, 5) * u.kg
    grid.add_quantity("test quantity", q1)

    # Test adding a quantity of incompatible size
    q2 = np.random.randn(5, 20, 5) * u.kg
    with pytest.raises(ValueError):
        grid.add_quantity("test quantity2", q2)


def test_CartesianGrid():

    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=10)

    x_arr, y_arr, z_arr = grid.grids
    x_axis, y_axis, z_axis = grid.ax0, grid.ax1, grid.ax2
    d_x, d_y, d_z = grid.dax0, grid.dax1, grid.dax2
    is_uniform_grid = grid.is_uniform_grid
    shape = grid.shape
    unit = grid.units

    # Grid should be uniform
    assert grid.is_uniform_grid == True

    # Test initializing with a provided grid
    grid2 = grids.CartesianGrid(grid.grids[0], grid.grids[1], grid.grids[2],)

    # Units not all consistent
    with pytest.raises(ValueError):
        grid = grids.CartesianGrid(
            [-1 * u.m, -1 * u.rad, -1 * u.m], [1 * u.m, 1 * u.rad, 1 * u.m]
        )


def test_interpolators():
    # Test interpolator

    # Test interpolator
    # Create grid
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=25)
    # Add some data to the grid
    data = grid.add_quantity("pts0", grid.grids[0])
    data = grid.add_quantity("pts1", grid.grids[1])

    pos = np.array([0.1, -0.3, 0]) * u.cm
    pos2 = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm

    # Interpolate indices
    i = grid.interpolate_indices(pos)[0]

    pout = grid.grid[i[0], i[1], i[2], :] * grid.unit

    # Assert that nearest grid cell was found
    assert np.allclose(pos, pout, atol=0.1)

    # Interpolate grid values using nearest-neighbor interpolator
    pout = grid.nearest_neighbor_interpolator(pos, "pts0")
    assert np.allclose(pos[0], pout, atol=0.1)

    # Interpolate grid values using volume-weighted interpolator
    pout = grid.volume_averaged_interpolator(pos, "pts0")
    assert np.allclose(pos[0], pout, atol=0.1)

    # Test using multiple arguments
    pout, pout2 = grid.nearest_neighbor_interpolator(pos, "pts0", "pts1")

    # Test multiple grid points
    pos = np.array([[0, 0, 0], [0.5, 0.2, 0]]) * u.cm
    i = grid.interpolate_indices(pos)
    pout = grid.nearest_neighbor_interpolator(pos.value, "pts0")
    pout = grid.volume_averaged_interpolator(pos, "pts0")


def test_NonUniformCartesianGrid():
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=10)

    pts0, pts1, pts2 = grid.grids
    shape = grid.shape
    units = grid.units

    # Grid should be non-uniform
    assert grid.is_uniform_grid == False

    # Test assigning a quantity
    q1 = np.random.randn(10, 10, 10) * u.kg
    grid.add_quantity("test_quantity", q1)
