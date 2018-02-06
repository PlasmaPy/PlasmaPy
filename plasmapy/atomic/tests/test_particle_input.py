import pytest
import warnings

from ...utils import AtomicError, InvalidParticleError
from ..particle_class import Particle
from ..particle_input import particle_input


@particle_input
def func_simple_noparens(
        a, particle: Particle, b=None, Z: int = None, mass_numb: int = None):
    r"""A simple function that, when decorated with `@particle_input`, returns
    the instance of the Particle class corresponding to the inputs."""
    if not isinstance(particle, Particle):
        raise TypeError(
            f"The argument particle in func_simple_noparens is not a Particle")
    return particle


@particle_input()
def func_simple_parens(
        a, particle: Particle, b=None, Z: int = None, mass_numb: int = None):
    r"""A simple function that, when decorated with `@particle_input()`,
    returns the instance of the Particle class corresponding to the inputs."""
    if not isinstance(particle, Particle):
        raise TypeError(
            f"The argument particle in func_simple_parens is not a Particle")
    return particle


particle_input_simple_table = [
    (func_simple_noparens, (1, 'p+'), {'b': 2}, 'p+'),
    (func_simple_parens, (1, 'p+'), {'b': 2}, 'p+'),
    (func_simple_noparens, (1, 'Fe'), {'mass_numb': 56, 'Z': 3}, 'Fe-56 3+'),
    (func_simple_parens, (1, 'Fe'), {'mass_numb': 56, 'Z': 3}, 'Fe-56 3+'),
    (func_simple_parens, (1,), {'particle': 'e-'}, 'e-'),
    (func_simple_noparens, (1,), {'particle': 'e-'}, 'e-'),
    (func_simple_noparens, (1,), {'particle': Particle('e-')}, 'e-'),
    (func_simple_parens, (1,), {'particle': Particle('e-')}, 'e-'),
]


@pytest.mark.parametrize(
    'func, args, kwargs, symbol', particle_input_simple_table)
def test_particle_input_simple(func, args, kwargs, symbol):
    r"""Test that simple functions decorated by particle_input correctly
    return the correct Particle object."""
    try:
        expected = Particle(symbol)
    except Exception as e:
        raise AtomicError(
            f"Cannot create Particle class from symbol {symbol}") from e

    try:
        result = func(*args, **kwargs)
    except Exception as e:
        raise AtomicError(
            f"An exception was raised while trying to execute "
            f"{func} with args = {args} and kwargs = {kwargs}.") from e

    assert result == expected, (
        f"The result {repr(result)} does not equal the expected value of "
        f"{repr(expected)}.\n\n"
        f"func = {func}\n"
        f"args = {args}\n"
        f"kwargs = {kwargs}\nsymbol = {symbol}\n"
        f"{result._attributes}\n"
        f"{expected._attributes}\n"
    )


# function, kwargs, expected_error
particle_input_error_table = [
    (func_simple_noparens, {'a': 1, 'particle': 'asdf'}, InvalidParticleError),
]


@pytest.mark.parametrize(
    'func, kwargs, expected_error', particle_input_error_table)
def test_particle_input_errors(func, kwargs, expected_error):
    r"""Test that functions decorated with particle_input raise the
    expected errors."""
    with pytest.raises(expected_error, message=(
            f"{func} did not raise {expected_error} with kwargs = {kwargs}")):
        func(**kwargs)


class Test_particle_input:

    @particle_input
    def method_noparens(self, particle: Particle):
        return particle

    @particle_input()
    def method_parens(self, particle: Particle):
        return particle


def test_particle_input_classes():
    instance = Test_particle_input()

    symbol = 'muon'
    expected = Particle(symbol)

    try:
        result_noparens = instance.method_noparens(symbol)
    except Exception as e:
        raise AtomicError("Problem with method_noparens") from e

    try:
        result_parens = instance.method_parens(symbol)
    except Exception as e:
        raise AtomicError("Problem with method_parens") from e

    assert result_parens == result_noparens == expected
