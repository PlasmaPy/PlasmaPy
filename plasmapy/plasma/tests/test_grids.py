"""
Tests for grids.py
"""

import astropy.units as u
import numpy as np
import pytest
import warnings

from numpy.testing import assert_allclose

from plasmapy.plasma import grids as grids


def test_AbstractGrid():

    # Create grid with single u.Quantity args
    grid = grids.AbstractGrid(-1 * u.cm, 1 * u.cm, num=10)
    assert grid.is_uniform

    # Create grid with lists of  u.Quantity args
    grid = grids.AbstractGrid(
        [-1 * u.cm, -1 * u.cm, -1 * u.cm],
        [1 * u.cm, 1 * u.cm, 1 * u.cm],
        num=[10, 10, 10],
    )
    assert grid.is_uniform

    # Create grid with arrays of u.quantities
    grid = grids.AbstractGrid(
        np.array([-1] * 3) * u.cm, np.array([1] * 3) * u.cm, num=[10, 10, 10]
    )
    assert grid.is_uniform

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

    # Test setting a subset of a quantity array
    grid["B_x"][0, 0, 0] = 21 * u.T
    assert grid["B_x"][0, 0, 0] == 21 * u.T

    # Test accessing a quantity using __getitem__ or directly
    Bx = grid.ds["B_x"]
    Bx = grid["B_x"]
    # Assert that the array returned is a u.Quantity
    assert isinstance(Bx, u.Quantity)
    # Assert that the array returned has the right shape
    assert Bx.shape == grid.shape

    # Test require_quantities
    # Test with a key that is there
    req_q = ["B_x"]
    grid.require_quantities(req_q, replace_with_zeros=False)
    req_q = ["B_x", "B_y"]
    # Test with a key that is not there, but can be replaced
    # Do not replace
    with pytest.raises(KeyError):
        grid.require_quantities(req_q, replace_with_zeros=False)
    # Do replace
    with pytest.warns(RuntimeWarning, match="This quantity will be assumed to be zero"):
        grid.require_quantities(req_q, replace_with_zeros=True)
    req_q = ["B_x", "B_y"]
    # Test with a key that is not there, but cannot be replaced because
    # it's not a recognized key
    req_q = ["B_x", "not_a_recognized_key"]
    with pytest.raises(KeyError):
        with pytest.warns(
            RuntimeWarning, match="This quantity will be assumed to be zero"
        ):
            grid.require_quantities(req_q, replace_with_zeros=True)

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
    is_uniform = grid.is_uniform
    shape = grid.shape
    unit = grid.units

    # Grid should be uniform
    assert grid.is_uniform == True

    # Test initializing with a provided grid and a quantity
    q = np.zeros(grid.shape)
    grid2 = grids.CartesianGrid(
        grid.grids[0], grid.grids[1], grid.grids[2], test_quantity=q
    )

    # Test that input with the wrong units will raise an exception
    L0 = [-1 * u.mm, 0 * u.rad, -1 * u.mm]
    L1 = [1 * u.mm, 2 * np.pi * u.rad, 1 * u.mm]
    with pytest.raises(ValueError):
        grid = grids.CartesianGrid(L0, L1, num=10)


def test_grid_methods():
    # ************ UNIFORM CARTESIAN ****************************
    grid = grids.CartesianGrid(
        np.array([-1, -1, -1]) * u.cm, np.array([1, 1, 1]) * u.cm, num=(10, 10, 10)
    )

    # Test on-grid
    pos = np.array([[0.1, -0.3, 0], [3, -0.3, 0]]) * u.cm
    out = grid.on_grid(pos)
    assert np.all(out == np.array([True, False]))

    # Test vector_intersects
    # This vector passes through the grid
    p1, p2 = np.array([0, -5, 0]) * u.cm, np.array([0, 5, 0]) * u.cm
    assert grid.vector_intersects(p1, p2) == True
    # Test going backwards yields the same result
    assert grid.vector_intersects(p2, p1) == True
    # This one doesn't
    p1, p2 = np.array([0, -5, 0]) * u.cm, np.array([0, -5, 10]) * u.cm
    assert grid.vector_intersects(p1, p2) == False
    assert grid.vector_intersects(p2, p1) == False

    # ************ NON-UNIFORM CARTESIAN ****************************

    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=10)

    pos = np.array([[0.1, -0.3, 0], [3, -0.3, 0]]) * u.cm
    out = grid.on_grid(pos)
    assert np.all(out == np.array([True, False]))


def test_interpolate_indices():
    # Create grid
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=25)

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    i = grid.interpolate_indices(pos)[0]
    # Assert that nearest grid cell was found
    pout = grid.grid[int(i[0]), int(i[1]), int(i[2])]
    assert np.allclose(pos, pout, atol=0.1)

    # One position, no units
    pos = np.array([0.1, -0.3, 0])
    i = grid.interpolate_indices(pos)[0]

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
    pout = grid.grid[int(i)]
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

    # Test persistance
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x", "y", persistent=True)
    pout = grid.nearest_neighbor_interpolator(pos, "x", "y", persistent=True)

    # ***********************************************************************

    # Create a non-uniform grid
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=100)

    print(grid.shape)

    grid.add_quantities(x=grid.grids[0], y=grid.grids[1])

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x")
    assert np.allclose(pos[0], pout, atol=0.5)

    # Test persistance
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x", "y", persistent=True)
    pout = grid.nearest_neighbor_interpolator(pos, "x", "y", persistent=True)


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

    # Try running with persistance
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    p1, p2 = grid.volume_averaged_interpolator(pos, "x", "y", persistent=True)
    p1, p2 = grid.volume_averaged_interpolator(pos, "x", "y", persistent=True)
    # Try changing the arg list, make sure it catchs this and auto-reverts
    # to non-persistent interpolation in that case
    p1, p2 = grid.volume_averaged_interpolator(pos, "x", persistent=True)
    assert p1.size == 1


def test_NonUniformCartesianGrid():
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=10)

    pts0, pts1, pts2 = grid.grids

    shape = grid.shape
    units = grid.units

    grid.add_quantities(x=pts0)
    print(grid)

    # Grid should be non-uniform
    assert grid.is_uniform == False

    # Test assigning a quantity
    q1 = np.random.randn(10, 10, 10) * u.kg / u.cm ** 3
    grid.add_quantities(rho=q1)

    # Test grid resolution for non-uniform grids
    assert 0 < grid.grid_resolution < 2

    # Test volume interpolator not implemented yet
    pos = np.array([5, -0.3, 0]) * u.cm
    with pytest.raises(NotImplementedError):
        pout = grid.volume_averaged_interpolator(pos, "x")

    # Test that many properties are unavailable
    with pytest.raises(ValueError):
        grid.ax0
    with pytest.raises(ValueError):
        grid.ax1
    with pytest.raises(ValueError):
        grid.ax2
    with pytest.raises(ValueError):
        grid.dax0
    with pytest.raises(ValueError):
        grid.dax1
    with pytest.raises(ValueError):
        grid.dax2

    # Test that input with the wrong units will raise an exception
    L0 = [-1 * u.mm, 0 * u.rad, -1 * u.mm]
    L1 = [1 * u.mm, 2 * np.pi * u.rad, 1 * u.mm]
    with pytest.raises(ValueError):
        grid = grids.NonUniformCartesianGrid(L0, L1, num=10)


if __name__ == "__main__":
    """
    test_AbstractGrid()
    test_CartesianGrid()
    test_grid_methods()
    test_interpolate_indices()
    test_nearest_neighbor_interpolator()
    test_volume_averaged_interpolator()
    test_NonUniformCartesianGrid()
    """
    pass
