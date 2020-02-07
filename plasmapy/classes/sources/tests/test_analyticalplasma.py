import pytest
from astropy import units as u
import numpy as np
from plasmapy.classes.sources import AnalyticalFields


@pytest.fixture()
def plasma():
    def magnetic_field(r):
        return u.Quantity([0, 0, 1], u.T)

    def electric_field(r):
        return u.Quantity([0, 1, 0], u.V / u.m)

    plasma = AnalyticalFields(magnetic_field, electric_field)
    return plasma


@pytest.fixture()
def particle_locations():
    np.random.seed(0)
    return u.Quantity(np.random.random((20, 3)), u.m)


def test_values(plasma, particle_locations):
    B = plasma.interpolate_B(particle_locations)
    E = plasma.interpolate_E(particle_locations)

    assert u.allclose(B, u.Quantity([0, 0, 1], u.T))
    assert u.allclose(E, u.Quantity([0, 1, 0], u.V / u.m))


def test_shapes(plasma, particle_locations):
    B = plasma.interpolate_B(particle_locations)
    E = plasma.interpolate_E(particle_locations)

    assert B.shape == particle_locations.shape
    assert E.shape == particle_locations.shape
