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
    xarr, yarr, zarr = grid.xarr, grid.yarr, grid.zarr
    radius = grid.distance_from_origin
    xaxis, yaxis, zaxis = grid.xaxis, grid.yaxis, grid.zaxis
    dx, dy, dz = grid.dx, grid.dy, grid.dz
    regular_grid = grid.regular_grid
    shape = grid.shape
    unit = grid.units

    # Test initializing with a provided grid
    grid2 = grids.CartesianGrid(grid.grid)

    # Test initializing using the unit keyword
    grid = grids.CartesianGrid(1, 1, units=u.cm)
    grid2 = grids.CartesianGrid(grid.grid.value, units=grid.units)


def test_CartesianGrid_exceptions():

    # Units not all consistent
    with pytest.raises(ValueError):
        grid = grids.CartesianGrid(-1, 1, units=[u.m, u.rad, u.rad])


test_CartesianGrid_exceptions()
