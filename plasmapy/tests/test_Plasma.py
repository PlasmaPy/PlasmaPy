import pytest
import numpy as np
import astropy.units as u

from plasmapy import plasma

@pytest.mark.parametrize('grid_dimensions, expected_size', [
    ((100, 1, 1), 100), # Test 1D setup
    ((128, 128, 1), 16384), # 2D
    ((64, 64, 64), 262144), # 3D
])
def test_plasma_setup(grid_dimensions, expected_size):
    x, y, z = grid_dimensions
    test_plasma = plasma.Plasma(domain_x=np.linspace(0, 1, x)*u.m,
                                domain_y=np.linspace(0, 1, y)*u.m,
                                domain_z=np.linspace(0, 1, z)*u.m)

    # Basic grid setup
    assert test_plasma.x.size == x
    assert test_plasma.y.size == y
    assert test_plasma.z.size == z
    assert test_plasma.grid.size == 3 * expected_size

    # Core variable units and shapes
    assert test_plasma.density.size == expected_size
    assert test_plasma.density.si.unit == u.kg / u.m**3

    assert test_plasma.momentum.size == 3 * expected_size
    assert test_plasma.momentum.si.unit == u.kg / (u.m**2 * u.s)

    assert test_plasma.pressure.size == expected_size
    assert test_plasma.pressure.si.unit == u.Pa

    assert test_plasma.magnetic_field.size == 3 * expected_size
    assert test_plasma.magnetic_field.si.unit == u.T


