"""Tests for ``@particle_input`` with delayed evaluation of annotations."""

# These tests must be in their own file because the `from __future__`
# import must be at the top.
from __future__ import annotations

from plasmapy.particles import Particle, ParticleLike, particle_input


@particle_input
def function_decorated_with_particle_input(particle: ParticleLike):
    return particle


class DecoratedClass:
    @particle_input
    def __init__(self, particle: ParticleLike):
        self.particle = particle


class UndecoratedClass:
    @particle_input
    def decorated_method(self, particle: ParticleLike):
        return particle


def test_particle_input_from_future_import_annotations_function():
    particle = function_decorated_with_particle_input("p+")
    assert particle == Particle("p+")


def test_particle_input_from_future_import_annotations_instantiation():
    instance = DecoratedClass("p+")

    assert isinstance(instance.particle, Particle)
    assert instance.particle == Particle("p+")


def test_particle_input_from_future_import_annotations_method():
    instance = UndecoratedClass()
    result = instance.decorated_method(particle="p+")
    assert isinstance(result, Particle)
    assert result == Particle("p+")
