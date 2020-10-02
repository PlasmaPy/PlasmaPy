"""
Tests for fields
"""

import astropy.units as u
import numpy as np
import pytest
import warnings

from plasmapy.plasma import fields as fields


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

test_posgrid()