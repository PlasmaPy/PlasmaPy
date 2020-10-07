"""
Tests for fields
"""

import astropy.units as u
import numpy as np
import pytest
import warnings

from plasmapy.plasma import fields as fields
from plasmapy.plasma import grids as grids


def test_example_fields():

    gridobj = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=(100, 100, 100))

    field_models = [
        "no fields",
        "electrostatic gaussian sphere",
        "electrostatic planar shock",
        "axial magnetic field",
    ]

    # Generate all the example fields
    for model in field_models:
        E, B = fields.example_fields(gridobj, model=model)

    # Test setting and getting emax and bmax
    model = fields.ElectrostaticGaussianSphere(gridobj, emax=1e9 * u.V / u.m)
    model.emax = 2e9 * u.V / u.m
    model.bmax = 1 * u.T
    a = model.emax
    a = model.bmax

    E, B = fields.example_fields(gridobj, model="electrostatic planar shock")


test_example_fields()
