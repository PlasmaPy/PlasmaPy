import pytest
import warnings

from ...utils import AtomicError
from ..particle_class import Particle
from ..particle_input import particle_input


@particle_input
def simple_decorated_function(particle: Particle):
    r"""Returns the appropriate instance of the Particle class."""

    if not isinstance(particle, Particle):
        warnings.warn(
            f"In simple_decorated_function, the input particle = {particle} "
            "is not an instance of the Particle class.  Instead")

    return particle


def test_simple_decorated_function():
    symbol = 'p+'
    expected = Particle(symbol)

    result = simple_decorated_function(symbol)

    if not isinstance(result, Particle):
        raise AtomicError(
            f"The result from simple_decorated_function is {repr(result)},"
            f" which should be an instance of the Particle class but is "
            f"actually of type {type(result)}.")

    assert result == expected, \
        (f"The result = {repr(result)} is not equal to expected = "
         f"{expected}.")


def test_class_method():

    class SomeClass:
        @particle_input
        def decorated_method(particle):
            if not isinstance(particle, Particle):
                warnings.warn("particle is not a Particle in "
                              "simpled_decorated_function")
            return particle

    someclass = SomeClass()
    symbol = 'p+'
    assert someclass.decorated_method(symbol) == Particle(symbol)
