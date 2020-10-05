"""
Tests for grids.py
"""

import astropy.units as u
import numpy as np
import pytest
import warnings

from plasmapy.plasma import grids as grids


def test_CartesianGrid():

    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm)

    xarr, yarr, zarr = grid.xarr, grid.yarr, grid.zarr
    radius = grid.distance_from_origin
    xaxis, yaxis, zaxis = grid.xaxis, grid.yaxis, grid.zaxis
    dx, dy, dz = grid.dx, grid.dy, grid.dz
    regular_grid = grid.regular_grid
    shape = grid.shape
    unit = grid.units

    grid2 = grids.CartesianGrid(grid.grid, units=grid.units)



def test_CartesianGrid_exceptions():

    # ****************************************
    # _load_grid

    # Wrong input
    with pytest.raises(ValueError):
        grid = grids.CartesianGrid('str')

    with pytest.raises(ValueError):
        grid = grids.CartesianGrid(1)


    # ****************************************
    # _make_grid

    # Wrong input
    with pytest.raises(ValueError):
        grid = grids.CartesianGrid(1,1)

    # Incompatable units
    with pytest.raises(ValueError):
        grid = grids.CartesianGrid(1*u.cm, 1*u.eV)


    # ****************************************
    # CartesianGrid-specific
    with pytest.raises(ValueError):
        grid = grids.CartesianGrid(1,1, units=(u.m, u.rad))
