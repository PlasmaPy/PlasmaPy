"""
Tests for fields
"""

import astropy.units as u
import numpy as np
import pytest
import warnings

from plasmapy.plasma import fields as fields


def test_posgrid():
    # Create a regular grid
    grid = fields.PosGrid(num=20, length=2 * u.cm, regular_grid=True)

    regular_grid = grid.regular_grid

    xarr, yarr, zarr = grid.xarr, grid.yarr, grid.zarr
    radius = grid.radius
    xaxis, yaxis, zaxis = grid.xaxis, grid.yaxis, grid.zaxis
    dx, dy, dz = grid.dx, grid.dy, grid.dz

    shape = grid.shape
    unit = grid.unit

    # Init a new grid object using the old grid and test for regularity etc.
    grid2 = fields.PosGrid(grid=grid.grid)
    regular_grid = grid2.regular_grid
    nearest_neighbor = grid2.nearest_neighbor

    pt1 = np.array([-10, 0, 0]) * u.cm
    pt2 = np.array([30, 0, 0]) * u.cm
    vector_intersects = grid2.vector_intersects(pt1, pt2)

    # Create an irregular grid
    grid = fields.PosGrid(num=(10, 20, 20), length=1 * u.cm, regular_grid=False)


def test_posgrid_raise_errors():
    # create irregular grid
    grid = fields.PosGrid(num=(10, 20, 20), length=1 * u.cm, regular_grid=False)

    with pytest.raises(ValueError):
        grid.xaxis
    with pytest.raises(ValueError):
        grid.yaxis
    with pytest.raises(ValueError):
        grid.zaxis
    with pytest.raises(ValueError):
        grid.dx
    with pytest.raises(ValueError):
        grid.dy
    with pytest.raises(ValueError):
        grid.dz

    # Throw an error grid is not a u.Quantity
    badgrid = np.zeros([10, 10, 10, 3])
    with pytest.raises(ValueError):
        grid = fields.PosGrid(grid=badgrid)

    # Throw an error if grid is not a valid shape
    badgrid = np.zeros([10, 10, 10])
    with pytest.raises(ValueError):
        grid = fields.PosGrid(grid=badgrid)


def test_example_fields():

    field_models = [
        "no fields",
        "electrostatic gaussian sphere",
        "electrostatic planar shock",
        "axial magnetic field",
    ]

    for model in field_models:
        grid, E, B = fields.example_fields(model=model)

    model = fields.ElectrostaticGaussianSphere(grid, emax=1e9 * u.V / u.m)

    # Test setting and getting emax and bmax
    model.emax = 2e9 * u.V / u.m
    model.bmax = 1 * u.T
    a = model.emax
    a = model.bmax

    grid, E, B = fields.example_fields(
        model="electrostatic planar shock", num=(100, 100, 200)
    )
