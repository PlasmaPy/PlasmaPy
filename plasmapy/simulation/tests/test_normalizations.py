import astropy.units as u
import pytest

from plasmapy.particles import proton
from plasmapy.simulation.normalizations import IdealMHDNormalizations, MHDNormalizations


def temporary_test():

    B = 1.1 * u.T
    L = 1.2 * u.m
    n = 1.3 * u.m ** -3
    ion = "p+"

    i = IdealMHDNormalizations(magnetic_field=B, length=L, number_density=n, ion=ion,)

    assert u.isclose(i.magnetic_field, B)
    assert u.isclose(i.magnetic_flux, B * L)
    assert i.ion == proton
