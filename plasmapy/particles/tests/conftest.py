import pytest

from plasmapy.particles import Particle
from plasmapy.particles._special_particles import particle_zoo
from plasmapy.particles.exceptions import InvalidParticleError


@pytest.fixture(params=list(sorted(particle_zoo.everything)))
def particle(request):
    return Particle(request.param)


@pytest.fixture()
def opposite(particle):
    try:
        opposite_particle = ~particle
    except Exception as exc:
        raise InvalidParticleError(
            f"The unary ~ (invert) operator is unable to find the "
            f"antiparticle of {particle}."
        ) from exc
    return opposite_particle
