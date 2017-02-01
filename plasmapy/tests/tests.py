import pytest
import numpy as np
import astropy.units as u

from plasmapy import plasma

def test_plasma_defaults():
    test_plasma = plasma.Plasma(domain_x=np.linspace(0, 1, 64)*u.m,
                                domain_y=np.linspace(0, 1, 64)*u.m,
                                domain_z=np.linspace(0, 1, 64)*u.m)

    # Basic grid setup
    assert test_plasma.x.shape == (64,)
    assert test_plasma.y.shape == (64,)
    assert test_plasma.z.shape == (64,)
    assert test_plasma.grid.shape == (3, 64, 64, 64)

    # Core variable units and shapes
    assert test_plasma.density.shape == (64, 64, 64)
    assert test_plasma.density.si.unit == u.kg / u.m**3

    assert test_plasma.momentum.shape == (3, 64, 64, 64)
    assert test_plasma.momentum.si.unit == u.kg / (u.m**2 * u.s)

    assert test_plasma.pressure.shape == (64, 64, 64)
    assert test_plasma.pressure.si.unit == u.Pa

    assert test_plasma.magnetic_field.shape == (3, 64, 64, 64)
    assert test_plasma.magnetic_field.si.unit == u.T
