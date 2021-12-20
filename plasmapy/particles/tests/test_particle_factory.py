import astropy.units as u
import pytest

from collections import namedtuple

from plasmapy.particles import deuteron, electron, proton
from plasmapy.particles.particle_class import (
    CustomParticle,
    DimensionlessParticle,
    particle,
    Particle,
)
from plasmapy.particles.particle_collections import ParticleList

test_case = namedtuple("test_case", ["args", "kwargs", "expected"])

mass = 1e-26 * u.kg
charge = 1e-29 * u.C
custom_particle = CustomParticle(mass, charge)
blank_custom_particle = CustomParticle()

test_cases = [
    test_case(args=(proton), kwargs={}, expected=proton),
    test_case(args=("p+"), kwargs=None, expected=proton),
    test_case(args=("H"), kwargs={"Z": 1, "mass_numb": 2}, expected=deuteron),
    test_case(args=("muon"), kwargs={}, expected=Particle("muon")),
    test_case(args=(charge, mass), kwargs={}, expected=custom_particle),
    test_case(args=(mass, charge), kwargs={}, expected=custom_particle),
    test_case(
        args=(), kwargs={"mass": mass, "charge": charge}, expected=custom_particle
    ),
    #    test_case(args=(), kwargs={}, expected=blank_custom_particle),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
]


@pytest.mark.parametrize("args, kwargs, expected", test_cases)
def test_particle_factory(args, kwargs, expected):
    result = particle(*args, **kwargs)
    assert result == expected


test_cases_for_exceptions = [
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
    #    test_case(args=(), kwargs={}, expected=),
]


@pytest.mark.parametrize("args, kwargs, expected", test_cases_for_exceptions)
def test_particle_factory_exceptions(args, kwargs, expected):
    with pytest.raises(expected):
        particle(*args, **kwargs)
