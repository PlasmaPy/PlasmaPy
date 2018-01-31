import pytest
import warnings

from ...utils import AtomicError, InvalidParticleError
from ..particle_class import Particle
from ..particle_input import particle_input


@particle_input
def func_simple_noparens(
        a, particle: Particle, b=None, Z: int = None, mass_numb: int = None):
    if not isinstance(particle, Particle):
        raise TypeError(
            f"The argument particle in func_simple_noparens is not a Particle")
    return particle


@particle_input()
def func_simple_parens(
        a, particle: Particle, b=None, Z: int = None, mass_numb: int = None):
    if not isinstance(particle, Particle):
        raise TypeError(
            f"The argument particle in func_simple_parens is not a Particle")
    return particle


particle_input_simple_table = [
    (func_simple_noparens, (1, 'p+'), {'b': 2}, 'p+'),
    (func_simple_parens, (1, 'p+'), {'b': 2}, 'p+'),
    (func_simple_noparens, (1, 'Fe'), {'mass_numb': 56, 'Z': 3}, 'Fe-56 3+'),
    (func_simple_parens, (1, 'Fe'), {'mass_numb': 56, 'Z': 3}, 'Fe-56 3+'),
]


@pytest.mark.parametrize(
    'func, args, kwargs, symbol', particle_input_simple_table)
def test_particle_input_simple(func, args, kwargs, symbol):

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
        f"{repr(expected)}.")
