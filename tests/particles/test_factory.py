import astropy.units as u
import pytest

from plasmapy.particles import deuteron, electron, proton
from plasmapy.particles._factory import _physical_particle_factory
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.particles.particle_class import CustomParticle, Particle
from plasmapy.particles.particle_collections import ParticleList

mass = 1e-26 * u.kg
charge = 1e-29 * u.C
custom_particle = CustomParticle(mass=mass, charge=charge)

He4 = Particle("He-4", Z=0)
Z_mean = 1.5

test_cases = [
    ([[]], {}, ParticleList()),
    ([proton], {}, proton),
    (["p+"], {}, proton),
    (["H"], {"Z": 1, "mass_numb": 2}, deuteron),
    (["muon"], {}, Particle("muon")),
    pytest.param([charge, mass], {}, custom_particle),
    ([mass, charge], {"Z": None, "mass_numb": None}, custom_particle),
    ([], {"symbol": "ξ"}, CustomParticle(symbol="ξ")),
    ([[proton, electron]], {}, ParticleList([proton, electron])),
    ([], {"mass": mass}, CustomParticle(mass=mass)),
    ([], {"charge": charge}, CustomParticle(charge=charge)),
    (["e-"], {}, electron),
    ([1], {}, Particle("H")),
    (["H"], {"Z": 1, "mass_numb": 1}, Particle("H", Z=1, mass_numb=1)),
    ([custom_particle], {}, custom_particle),
    ([ParticleList(["p+", "e-"])], {}, ParticleList(["p+", "e-"])),
    ([mass, charge], {}, custom_particle),
    ([charge, mass], {}, custom_particle),
    ([charge], {}, CustomParticle(charge=charge)),
    ([mass], {}, CustomParticle(mass=mass)),
    (
        ["He-4"],
        {"Z": Z_mean},
        CustomParticle(mass=He4.mass - Z_mean * electron.mass, Z=Z_mean),
    ),
    (
        ["He"],
        {"Z": Z_mean, "mass_numb": 4},
        CustomParticle(mass=He4.mass - Z_mean * electron.mass, Z=Z_mean),
    ),
]


@pytest.mark.parametrize(("args", "kwargs", "expected"), test_cases)
def test_physical_particle_factory(args, kwargs, expected) -> None:
    result = _physical_particle_factory(*args, **kwargs)
    assert result == expected
    assert type(result) is type(expected)


test_cases_for_exceptions = [
    ([], {}, TypeError),
    ("not a valid Particle", {}, InvalidParticleError),
    (["not valid for a ParticleList"], {}, InvalidParticleError),
    (["He-4"], {"Z": 2.001}, InvalidParticleError),
    (["He-4"], {"Z": 1 + 1j}, InvalidParticleError),
    (["tau neutrino"], {"Z": 1.3}, InvalidParticleError),
]


@pytest.mark.parametrize(("args", "kwargs", "expected"), test_cases_for_exceptions)
def test_particle_factory_exceptions(args, kwargs, expected) -> None:
    with pytest.raises(expected):
        _physical_particle_factory(*args, **kwargs)


def test_particle_factory_custom_particle_with_none_kwargs() -> None:
    """
    Test that when `_physical_particle_factory` is provided with a
    `CustomParticle` along with ``Z=None`` and ``mass_numb=None`` as
    keyword arguments, then it will return the `CustomParticle`.
    """
    expected = CustomParticle(mass=1.27 * u.kg, charge=1 * u.C)
    actual = _physical_particle_factory(expected, Z=None, mass_numb=None)
    assert expected is actual
