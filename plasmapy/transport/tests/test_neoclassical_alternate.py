import astropy.units as u
import datetime
import hypothesis
import numpy as np
import pytest

from astropy.tests.helper import assert_quantity_allclose
from hypothesis import example, given, settings
from hypothesis import strategies as st

from plasmapy.particles import IonizationStateCollection, Particle
from plasmapy.transport.Houlberg1997 import ExtendedParticleList

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
import numpy as np

from plasmapy.transport.Houlberg1997 import ExtendedParticleList, Particle
all_species = ExtendedParticleList([Particle("p+"), Particle("C +6")],
                     10 * u.eV,
                     u.Quantity([1, 0.08], u.m**-3) * 1e20,
                     u.Quantity(2*[-1], u.eV / u.m),
                     u.Quantity(2*[1e18], u.m**-4),
                    )
N = len(all_species)
M = len(all_species.basic_elements)
thermal_speed = u.Quantity([12676.973, 43769.471], "m/s")
mass = u.Quantity([1.9939269e-26, 1.6727219e-27], "kg")

@pytest.mark.parametrize(
    ["field", "quantity", "tolerances"],
    [
        ("dT", u.Quantity(2*[-11604.518], "K/m"), {}),
        ("mass_density", u.Quantity([1.5951414e-7, 1.6726219e-7], "kg/m3"), {}),
        ("thermal_speed", thermal_speed, {}),
        ("ξ", u.Quantity([1, 1]), {}),
        ("isotopic_thermal_speed", thermal_speed, dict(rtol=1e-4)),
        ("mass", mass, dict(rtol=1e-4)),
        ("isotopic_mass", mass, dict(rtol=1e-4)), 
    ],
)
def test_quantity(field, quantity, tolerances):
    data = getattr(all_species, field)
    assert_quantity_allclose(data, quantity, **tolerances)

@pytest.mark.parametrize(
    ["field", "quantity", "tolerances"],
    [
        ("temperature_ratio", u.Quantity(2*[2*[1]]), {}),
        ("mass_ratio", u.Quantity([[1, 11.920966], [0.083885819, 1]]), dict(rtol=1e-4)),
    ],
)
def test_symmetric_quantity(field, quantity, tolerances):
    data = getattr(all_species, field)
    assert_quantity_allclose(data, quantity, **tolerances)
    assert_quantity_allclose(
        np.diag(np.fliplr(data)),
        np.diag(np.fliplr(data))[::-1]**-1
    )


def is_sorted(arr):
    return (arr.flatten() == np.sort(arr.flatten())).all()

def test_compress():
    N = len(all_species)
    arr = np.arange(N ** 2).reshape(N, N)
    assert is_sorted(arr)
    compressed = all_species.compress(arr, axis=1)
    assert is_sorted(compressed)



