from __future__ import annotations

from plasmapy.particles import Particle, ParticleLike, particle_input


@particle_input
def function_decorated_with_particle_input(x, particle: ParticleLike):
    return particle


class SomeClass:
    @particle_input
    def __init__(self, particle: ParticleLike):
        self.particle = particle

    @particle_input
    def f(x, particle: ParticleLike):
        return particle


def test_particle_input_from_future_import_annotations_function():
    particle = f(None, "p+")
    assert particle == Particle("p+")


def test_particle_input_from_future_import_annotations_method():
    instance = SomeClass(1, "p+")

    assert isinstance(instance.particle, Particle)
    assert instance.particle == "p+"
