import pytest


from ...utils import AtomicError
from ..particle_class import Particle
from ..particle_decorator import particle_input


class Test_particle_input_decorator:

    @particle_input
    def simple_decorated_function(particle, *args, **kwargs):
        r"""Returns the appropriate instance of the Particle class."""
        return particle

    def test_simple_decorated_function(self):
        symbol = 'p+'
        expected = Particle(symbol)

        try:
            result = self.simple_decorated_function(symbol)
        except Exception as e:
            raise AtomicError(
                "simple_decorated_function did not work.") from e

        if not isinstance(result, Particle):
            raise AtomicError(
                f"The result from simple_decorated_function is {repr(result)},"
                f" which should be an instance of the Particle class but is "
                f"actually of type {type(result)}."
            )

        assert result == expected, \
            (f"The result = {repr(result)} is not equal to expected = "
             f"{expected}.  ")
