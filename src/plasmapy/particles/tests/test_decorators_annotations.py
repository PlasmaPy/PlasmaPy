"""Tests for ``@particle_input`` with postponed evaluation of annotations."""

# These tests must be in their own file because the `from __future__`
# import must be at the top.
from __future__ import annotations

from plasmapy.particles.decorators import particle_input
from plasmapy.particles.particle_class import Particle, ParticleLike


@particle_input
def function_decorated_with_particle_input(particle: ParticleLike) -> Particle:
    return particle  # type: ignore[return-value]


class DecoratedClass:
    @particle_input
    def __init__(self, particle: ParticleLike) -> None:
        self.particle = particle


class UndecoratedClass:
    @particle_input
    def decorated_method(self, particle: ParticleLike) -> Particle:
        return particle  # type: ignore[return-value]


def test_particle_input_postponed_annotations_function() -> None:
    """
    Test that `@particle_input` works with postponed evaluation of
    annotations for a simple function.
    """
    particle = function_decorated_with_particle_input("p+")
    assert particle == Particle("p+")


def test_particle_input_postponed_annotations_instance() -> None:
    """
    Test that `@particle_input` works with postponed evaluation of
    annotations when it decorates the `__init__` method of a class.
    """
    instance = DecoratedClass("p+")
    assert isinstance(instance.particle, Particle)
    assert instance.particle == Particle("p+")


def test_particle_input_postponed_annotations_method() -> None:
    """
    Test that `@particle_input` works with postponed evaluation of
    annotations when it decorates an instance method.
    """
    instance = UndecoratedClass()
    result = instance.decorated_method(particle="p+")
    assert isinstance(result, Particle)
    assert result == Particle("p+")
