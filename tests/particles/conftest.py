"""Define fixtures for tests of `plasmapy.particles`."""

import pytest

from plasmapy.particles import Particle
from plasmapy.particles._special_particles import particle_zoo
from plasmapy.particles.exceptions import InvalidParticleError


@pytest.fixture(params=sorted(particle_zoo.everything))
def particle(request):  # noqa: D103
    return Particle(request.param)


@pytest.fixture()
def opposite(particle):  # noqa: D103
    try:
        opposite_particle = ~particle
    except Exception as exc:  # noqa: BLE001
        raise InvalidParticleError(
            f"The unary ~ (invert) operator is unable to find the "
            f"antiparticle of {particle}."
        ) from exc
    return opposite_particle


@pytest.fixture(
    params=sorted(
        [
            ("e-", "e+"),
            ("mu-", "mu+"),
            ("tau-", "tau+"),
            ("p+", "p-"),
            ("n", "antineutron"),
            ("nu_e", "anti_nu_e"),
            ("nu_mu", "anti_nu_mu"),
            ("nu_tau", "anti_nu_tau"),
        ]
    )
)
def particle_antiparticle_pair(request):  # noqa: D103
    return [Particle(p) for p in request.param]
