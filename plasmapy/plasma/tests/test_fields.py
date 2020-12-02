"""
Tests for fields
"""

import astropy.units as u
import numpy as np
import pytest
import warnings

from plasmapy.plasma import fields as fields
from plasmapy.plasma import grids as grids


@pytest.fixture
def gridobj():
    gridobj = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=(100, 100, 100))
    return gridobj


@pytest.mark.parametrize(
    "model",
    [
        "no fields",
        "electrostatic gaussian sphere",
        "electrostatic planar shock",
        "axial magnetic field",
    ],
)
def test_example_fields(gridobj, model):
    E, B = fields.example_fields(gridobj, model=model)


def test_setting_emax_bmax(gridobj):
    model = fields.ElectrostaticGaussianSphere(gridobj, emax=1e9 * u.V / u.m)
    model.emax = 2e9 * u.V / u.m
    model.bmax = 1 * u.T
    a = model.emax
    a = model.bmax


def test_planar_shock(gridobj):
    E, B = fields.example_fields(gridobj, model="electrostatic planar shock")


if __name__ == "__main__":
    test_example_fields()
