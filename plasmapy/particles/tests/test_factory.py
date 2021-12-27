import astropy.units as u
import pytest

from plasmapy.particles import deuteron, electron, proton
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.particles.factory import physical_particle_factory
from plasmapy.particles.particle_class import (
    CustomParticle,
    DimensionlessParticle,
    Particle,
)
from plasmapy.particles.particle_collections import ParticleList

mass = 1e-26 * u.kg
charge = 1e-29 * u.C
custom_particle = CustomParticle(mass, charge)

test_cases = [
    ([[]], {}, ParticleList()),
    ([proton], {}, proton),
    (["p+"], {}, proton),
    (["H"], {"Z": 1, "mass_numb": 2}, deuteron),
    (["muon"], {}, Particle("muon")),
    ([charge, mass], {}, custom_particle),
    ([mass, charge], {}, custom_particle),
    ([], {"symbol": "ξ"}, CustomParticle(symbol="ξ")),
    ([[proton, electron]], {}, ParticleList([proton, electron])),
    ([], {"mass": mass}, CustomParticle(mass=mass)),
    ([], {"charge": charge}, CustomParticle(charge=charge)),
    (["e-"], {}, electron),
    ([1], {}, Particle("H")),
    (["H"], {"Z": 1, "mass_numb": 1}, Particle("H", Z=1, mass_numb=1)),
    ([custom_particle], {}, custom_particle),
    ([ParticleList(["p+", "e-"])], {}, ParticleList(["p+", "e-"])),
]


@pytest.mark.parametrize("args, kwargs, expected", test_cases)
def test_physical_particle_factory(args, kwargs, expected):
    result = physical_particle_factory(*args, **kwargs)
    assert result == expected
    assert type(result) == type(expected)


test_cases_for_exceptions = [
    ([], {}, TypeError),
    ([DimensionlessParticle()], {}, TypeError),
    ("...", {}, InvalidParticleError),
    (["..."], {}, InvalidParticleError),
]


@pytest.mark.parametrize("args, kwargs, expected", test_cases_for_exceptions)
def test_particle_factory_exceptions(args, kwargs, expected):
    with pytest.raises(expected):
        physical_particle_factory(*args, **kwargs)
