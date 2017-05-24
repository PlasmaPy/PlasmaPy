import pytest
import numpy as np
import astropy.units as u

from plasmapy.classes import plasma


@pytest.mark.parametrize('grid_dimensions, expected_size', [
    ((100, 1, 1), 100),  # Test 1D setup
    ((128, 128, 1), 16384),  # 2D
    ((64, 64, 64), 262144),  # 3D
])
def test_Plasma_setup(grid_dimensions, expected_size):
    """Function to test basic setup of the Plasma object.

    Tests that a Plasma object initiated with a particular specification
    behaves in the correct way.

    Parameters
    ----------
    grid_dimensions : tuple of ints
        Grid size of the Plasma object to test. Must be a tuple of length 3,
        indicating length of the grid in x, y, and z directions respectively.
        Directions not needed should have a length of 1.

    expected_size : int
        Product of grid dimensions.

    Examples
    --------
    >>> test_Plasma_setup((10, 10, 10), 1000)
    >>> test_Plasma_setup((100, 10, 1), 1000)
    """
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

    assert test_plasma.energy.size == expected_size
    assert test_plasma.energy.si.unit == u.J / u.m**3

    assert test_plasma.magnetic_field.size == 3 * expected_size
    assert test_plasma.magnetic_field.si.unit == u.T


# @pytest.mark.parametrize([()])
def test_Plasma_derived_vars():
    """Function to test derived variables of the Plasma class.

    Tests the shapes, units and values of variables derived from core
    variables.
    The core variables are set with arbitrary uniform values.
    """
    test_plasma = plasma.Plasma(domain_x=np.linspace(0, 1, 8)*u.m,
                                domain_y=np.linspace(0, 1, 8)*u.m,
                                domain_z=np.linspace(0, 1, 1)*u.m)

    # Set an arbitrary uniform values throughout the plasma
    test_plasma.density = 2.0 * np.ones((8, 8)) * u.kg / u.m**3
    test_plasma.momentum = 10.0 * np.ones((3, 8, 8)) * u.kg / (u.m**2 * u.s)
    test_plasma.energy = np.ones((8, 8)) * u.J / u.m**3
    test_plasma.magnetic_field = 0.01 * np.ones((3, 8, 8)) * u.T

    # Test derived variable units and shapes
    assert test_plasma.velocity.shape == test_plasma.momentum.shape
    assert (test_plasma.velocity == 5.0 * u.m / u.s).all()

    assert test_plasma.magnetic_field_strength.shape == test_plasma.magnetic_field.shape[1:]
    assert test_plasma.magnetic_field_strength.si.unit == u.T
    assert np.allclose(test_plasma.magnetic_field_strength.value, 0.017320508)

    assert test_plasma.alfven_speed.shape == test_plasma.density.shape
    assert test_plasma.alfven_speed.unit.si == u.m / u.s
    assert np.allclose(test_plasma.alfven_speed.value, 10.92548431)
