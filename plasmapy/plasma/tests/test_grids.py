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
    ax0, ax1, ax2 = grid.ax0, grid.ax1, grid.ax2
    dax1, dax2, dax3 = grid.dax0, grid.dax1, grid.dax2
    arr0, arr1, arr2 = grid.arr0, grid.arr1, grid.arr2


def test_AbstractGrid_too_many_positionals():
    with pytest.raises(TypeError):
        grid = grids.AbstractGrid(1 * u.cm, 1 * u.eV, 1)

    # ****************************************
    # _load_grid

    # ****************************************
    # _make_grid

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
    grid2 = grids.CartesianGrid(grid.grid)

    # Test initializing using the unit keyword
    grid = grids.CartesianGrid(1, 1, units=u.cm)
    grid2 = grids.CartesianGrid(grid.grid.value, units=grid.units)

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
    pos = np.array([0.1, -0.3, 0]) * u.cm

    # Interpolate indices
    i = grid.interpolate_indices(pos)[0]
    pout = grid.grid[i[0], i[1], i[2], :]
    # Assert that nearest grid cell was found
    assert np.allclose(pos, pout, atol=0.03)

    # Interpolate grid values using nearest-neighbor interpolator
    pout = grid.nearest_neighbor_interpolator(pos, grid.grid)
    assert np.allclose(pos, pout, atol=0.03)

    # Interpolate grid values using volume-weighted interpolator
    pout = grid.volume_averaged_interpolator(pos, grid.grid)
    assert np.allclose(pos, pout, atol=0.03)

    # Test using multiple arguments
    pout, pout2 = grid.nearest_neighbor_interpolator(pos, grid.grid, grid.grid * 2)

    # Test a u-quantity input
    i = grid.interpolate_indices(pos.value)
    pout = grid.nearest_neighbor_interpolator(pos.value, grid.grid)
    # Test multiple grid points
    pos = np.array([[0, 0, 0], [0.5, 0.2, 0]]) * u.cm
    i = grid.interpolate_indices(pos)
    pout = grid.nearest_neighbor_interpolator(pos.value, grid.grid)
    pout = grid.volume_averaged_interpolator(pos, grid.grid)


def test_NonUniformCartesianGrid():
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm)

    # Grid should be non-uniform
    assert grid.is_uniform_grid == False

    pos = np.array([0.1, -0.3, 0]) * u.cm
    # Test interpolator on non-uniform grid
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=100)
    i = grid.interpolate_indices(pos)
    # Get position value corresponding to that index
    pout = grid.grid[i[0][0], i[0][1], i[0][2], :]
    # Assert that grid is sort-of nearby...
    assert np.allclose(pos, pout, atol=0.1)


def test_CylindricalGrid():

    grid = grids.CylindricalGrid(
        [-1, 0, -1], [1, 2 * np.pi, 1], units=[u.cm, u.rad, u.cm]
    )

    array = grid.grid
    rho_arr, theta_arr, z_arr = grid.rho_arr, grid.theta_arr, grid.z_arr
    rho_axis, theta_axis, z_axis = grid.rho_axis, grid.theta_axis, grid.z_axis
    d_rho, d_theta, d_z = grid.d_rho, grid.d_theta, grid.d_z
    is_uniform_grid = grid.is_uniform_grid
    shape = grid.shape
    unit = grid.units

    # Units not all consistent
    with pytest.raises(ValueError):
        grid = grids.CylindricalGrid(-1, 1, units=[u.m, u.m, u.m])
    with pytest.raises(ValueError):
        grid = grids.CylindricalGrid(-1, 1, units=[u.m, u.rad, u.rad])

    # Theta axis larger than 2pi
    with pytest.raises(ValueError):
        grid = grids.CylindricalGrid(
            [-1, 0, -1], [1, 6 * np.pi, 1], units=[u.cm, u.rad, u.cm]
        )

    # Test interpolaotr
    pos = np.array([0.1, 0.1, 0])
    i = grid.interpolate_indices(pos)
    # Get position value corresponding to that index
    pout = grid.grid[i[0][0], i[0][1], i[0][2], :]
    # Assert that grid is sort-of nearby...
    assert np.allclose(pos, pout, atol=0.05)


def test_SphericalGrid():

    grid = grids.SphericalGrid(
        [-1, 0, 0], [1, np.pi, 2 * np.pi], units=[u.cm, u.rad, u.rad]
    )

    array = grid.grid
    r_arr, theta_arr, phi_arr = grid.r_arr, grid.theta_arr, grid.phi_arr
    r_axis, theta_axis, phi_axis = grid.r_axis, grid.theta_axis, grid.phi_axis
    d_r, d_theta, d_phi = grid.d_r, grid.d_theta, grid.d_phi
    is_uniform_grid = grid.is_uniform_grid
    shape = grid.shape
    unit = grid.units

    # Units not all consistent
    with pytest.raises(ValueError):
        grid = grids.SphericalGrid(-1, 1, units=[u.m, u.m, u.m])
    with pytest.raises(ValueError):
        grid = grids.SphericalGrid(-1, 1, units=[u.m, u.rad, u.m])

    # Theta axis larger than pi
    with pytest.raises(ValueError):
        grid = grids.SphericalGrid(
            [-1, 0, 0], [1, 1.5 * np.pi, 2 * np.pi], units=[u.cm, u.rad, u.rad]
        )
    # Phi axis larger than 2pi
    with pytest.raises(ValueError):
        grid = grids.SphericalGrid(
            [-1, 0, 0], [1, np.pi, 6 * np.pi], units=[u.cm, u.rad, u.rad]
        )

    # Test interpolaotr
    pos = np.array([0.1, 0.1, 0.5])
    i = grid.interpolate_indices(pos)
    # Get position value corresponding to that index
    pout = grid.grid[i[0][0], i[0][1], i[0][2], :]
    # Assert that grid is sort-of nearby...
    assert np.allclose(pos, pout, atol=0.05)
