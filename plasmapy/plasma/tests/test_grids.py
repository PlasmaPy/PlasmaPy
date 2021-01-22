"""
Tests for grids.py
"""

import astropy.units as u
import numpy as np
import pytest
import warnings

from plasmapy.plasma import grids as grids


def test_AbstractGrid():

    # Create grid with single u.Quantity args
    grid = grids.AbstractGrid(-1 * u.cm, 1 * u.cm, num=10)
    assert grid.is_uniform_grid

    # Create grid with lists of  u.Quantity args
    grid = grids.AbstractGrid(
        [-1 * u.cm, -1 * u.cm, -1 * u.cm],
        [1 * u.cm, 1 * u.cm, 1 * u.cm],
        num=[10, 10, 10],
    )
    assert grid.is_uniform_grid

    # Create grid with arrays of u.quantities
    grid = grids.AbstractGrid(
        np.array([-1] * 3) * u.cm, np.array([1] * 3) * u.cm, num=[10, 10, 10]
    )
    assert grid.is_uniform_grid

    print(grid)

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
    q = np.random.randn(10, 10, 10) * u.T
    grid.add_quantities(B_x=q)

    # Test accessing a quantity using __getitem__ or directly
    Bx = grid.ds["B_x"]
    Bx = grid["B_x"]

    # Test adding a quantity with wrong units
    q = np.random.randn(10, 10, 10) * u.kg
    with pytest.raises(ValueError):
        grid.add_quantities(B_x=q)

    # Testing adding a quantity with an unrecognized key name
    with pytest.warns(UserWarning):
        grid.add_quantities(not_a_recognized_key=q)

    # Test adding a quantity of incompatible size
    q = np.random.randn(5, 20, 5) * u.T
    with pytest.raises(ValueError):
        grid.add_quantities(B_x=q)

    # Test adding multiple quantites at once
    q = np.random.randn(10, 10, 10) * u.T
    grid.add_quantities(B_x=q, B_y=q, B_z=q)

    print(grid)


def test_CartesianGrid():

    grid = grids.CartesianGrid(
        np.array([-1, -1, -1]) * u.cm, np.array([1, 1, 1]) * u.cm, num=(10, 10, 10)
    )

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


def test_grid_methods():
    grid = grids.CartesianGrid(
        np.array([-1, -1, -1]) * u.cm, np.array([1, 1, 1]) * u.cm, num=(10, 10, 10)
    )


    # Test on-grid
    pos = np.array([[0.1, -0.3, 0], [3, -0.3, 0]]) * u.cm
    out = grid.on_grid(pos)
    assert np.all(out == np.array([True, False]))


    # Test vector_intersects
    # This vector passes through the grid
    p1, p2 = np.array([0,-5,0])*u.cm, np.array([0,5,0])*u.cm
    assert grid.vector_intersects(p1, p2) == True
    # This one doesn't
    p1, p2 = np.array([0,-5,0])*u.cm, np.array([0,-5,10])*u.cm
    assert grid.vector_intersects(p1, p2) == False



def test_interpolate_indices():
    # Create grid
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=25)

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    i = grid.interpolate_indices(pos)[0]
    # Assert that nearest grid cell was found
    pout = grid.grid[int(i[0]), int(i[1]), int(i[2])] * grid.unit
    assert np.allclose(pos, pout, atol=0.1)

    # Two positions
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    i = grid.interpolate_indices(pos)[0]

    # Contains out-of-bounds values (index array should contain NaNs)
    pos = np.array([5, -0.3, 0]) * u.cm
    i = grid.interpolate_indices(pos)[0]
    assert np.sum(np.isnan(i)) > 0

    # ***********************************************************************

    # Create a non-uniform grid
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=100)

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    i = grid.interpolate_indices(pos)[0]

    # Assert that nearest grid cell was found
    pout = grid.grid[int(i)] * grid.unit
    assert np.allclose(pos, pout, atol=0.5)


def test_nearest_neighbor_interpolator():
    # Create grid
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=25)
    # Add some data to the grid
    grid.add_quantities(x=grid.grids[0])
    grid.add_quantities(y=grid.grids[1])

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x")
    assert np.allclose(pos[0], pout, atol=0.1)

    # Test quantity key not present in dataset
    with pytest.raises(KeyError):
        pout = grid.nearest_neighbor_interpolator(pos, "B_x")

    # Two positions, two quantities
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x", "y")

    # Contains out-of-bounds values (must handle NaNs correctly)
    pos = np.array([5, -0.3, 0]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x")
    assert np.allclose(pout, 0 * u.cm, atol=0.1)

    # ***********************************************************************

    # Create a non-uniform grid
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=100)
    grid.add_quantities(x=grid.grids[0])

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x")
    assert np.allclose(pos[0], pout, atol=0.5)


def test_volume_averaged_interpolator():
    # Create grid
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=25)
    # Add some data to the grid
    grid.add_quantities(x=grid.grids[0])
    grid.add_quantities(y=grid.grids[1])

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    pout = grid.volume_averaged_interpolator(pos, "x")
    assert np.allclose(pos[0], pout, atol=0.1)

    # Test quantity key not present in dataset
    with pytest.raises(KeyError):
        pout = grid.volume_averaged_interpolator(pos, "B_x")

    # Two positions, two quantities
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    pout = grid.volume_averaged_interpolator(pos, "x", "y")

    # Contains out-of-bounds values (must handle NaNs correctly)
    pos = np.array([5, -0.3, 0]) * u.cm
    pout = grid.volume_averaged_interpolator(pos, "x")
    assert np.allclose(pout, 0 * u.cm, atol=0.1)


def test_NonUniformCartesianGrid():
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=10)

    pts0, pts1, pts2 = grid.grids
    shape = grid.shape
    units = grid.units

    # Grid should be non-uniform
    assert grid.is_uniform_grid == False

    # Test assigning a quantity
    q1 = np.random.randn(10, 10, 10) * u.kg / u.cm ** 3
    grid.add_quantities(rho=q1)


if __name__ == "__main__":
    #test_AbstractGrid()
    #test_CartesianGrid()
    test_grid_methods()
    # test_interpolate_indices()
    #test_nearest_neighbor_interpolator()
    #test_volume_averaged_interpolator()
    # test_NonUniformCartesianGrid()
    pass
