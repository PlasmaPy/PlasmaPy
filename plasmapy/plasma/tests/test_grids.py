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


def test_AbstractGrid_exceptions():
    # ****************************************
    # __init__

    # Too many positional arguments
    with pytest.raises(TypeError):
        grid = grids.AbstractGrid(1 * u.cm, 1 * u.eV, 1)

    # ****************************************
    # _load_grid

    # ****************************************
    # _make_grid

    # Incompatable units
    with pytest.raises(ValueError):
        grid = grids.AbstractGrid(1 * u.cm, 1 * u.eV)

    # Invalid length of start, stop, num, or units
    with pytest.raises(ValueError):
        grid = grids.AbstractGrid(-1, [2, 3], units=[u.m, u.m])


def test_CartesianGrid():

    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm)

    array = grid.grid
    x_arr, y_arr, z_arr = grid.x_arr, grid.y_arr, grid.z_arr
    radius = grid.distance_from_origin
    x_axis, y_axis, z_axis = grid.x_axis, grid.y_axis, grid.z_axis
    d_x, d_y, d_z = grid.d_x, grid.d_y, grid.d_z
    regular_grid = grid.regular_grid
    shape = grid.shape
    unit = grid.units

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


def test_IrregularCartesianGrid():
    grid = grids.IrregularCartesianGrid(-1 * u.cm, 1 * u.cm)

    # Grid should be irregular
    assert grid.regular_grid == False


def test_CylindricalGrid():

    grid = grids.CylindricalGrid(
        [-1, 0, -1], [1, 2 * np.pi, 1], units=[u.cm, u.rad, u.cm]
    )

    array = grid.grid
    rho_arr, theta_arr, z_arr = grid.rho_arr, grid.theta_arr, grid.z_arr
    rho_axis, theta_axis, z_axis = grid.rho_axis, grid.theta_axis, grid.z_axis
    d_rho, d_theta, d_z = grid.d_rho, grid.d_theta, grid.d_z
    regular_grid = grid.regular_grid
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


def test_SphericalGrid():

    grid = grids.SphericalGrid(
        [-1, 0, 0], [1, np.pi, 2 * np.pi], units=[u.cm, u.rad, u.rad]
    )

    array = grid.grid
    r_arr, theta_arr, phi_arr = grid.r_arr, grid.theta_arr, grid.phi_arr
    r_axis, theta_axis, phi_axis = grid.r_axis, grid.theta_axis, grid.phi_axis
    d_r, d_theta, d_phi = grid.d_r, grid.d_theta, grid.d_phi
    regular_grid = grid.regular_grid
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
