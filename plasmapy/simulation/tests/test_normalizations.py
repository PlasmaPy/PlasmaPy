from astropy import units as u, constants as const
import pytest

from plasmapy.particles import proton
from plasmapy.simulation.normalizations import AbstractNormalizations, MHDNormalizations


def temporary_test():

    B = 1.1 * u.T
    L = 1.2 * u.m
    n = 1.3 * u.m ** -3
    ion = "p+"

    i = IdealMHDNormalizations(magnetic_field=B, length=L, number_density=n, ion=ion,)

    assert u.isclose(i.magnetic_field, B)
    assert u.isclose(i.magnetic_flux, B * L)
    assert i.ion == proton


class MockNormalizations(AbstractNormalizations):
    """MHD normalizations that were calculated by hand."""
    magnetic_field = 2.0 * u.T
    number_density = 3.0 * u.m ** -3
    length = 5.0 * u.m
    current_density = magnetic_field /


@pytest.mark.fixture
def mock_normalizations():
    return MockNormalizations()


