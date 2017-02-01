import pytest
import numpy as np
import astropy.units as u

from plasmapy import plasma

def test_plasma_defaults():
    test_plasma = plasma.Plasma(domain_x=np.linspace(0, 1, 64)*u.m,
                                domain_y=np.linspace(0, 1, 64)*u.m,
                                domain_z=np.linspace(0, 1, 64)*u.m)

    assert test_plasma.x.shape == (64,)
    assert test_plasma.y.shape == (64,)
    assert test_plasma.z.shape == (64,)
    assert test_plasma.grid.shape == (3, 64, 64, 64)
