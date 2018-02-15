import pytest
import warnings

from ...utils import (
    AtomicError,
    InvalidParticleError,
    ChargeError,
    InvalidElementError,
    InvalidIonError,
    InvalidIsotopeError,
)

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
    r"""A sample class with methods to make sure that """

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


# decorator_kwargs, particle, expected_exception
decorator_categories_table = [
    ({'exclude': {'element'}}, 'Fe', AtomicError),
    ({'any_of': {'lepton', 'antilepton'}}, 'tau-', None),
    ({'require': {'isotope', 'ion'}}, 'Fe-56+', None),
    ({'require': {'isotope', 'ion'}}, 'Fe+', AtomicError),
    ({'any_of': {'isotope', 'ion'}}, 'Fe+', None),
    ({'any_of': {'charged', 'uncharged'}}, 'Fe', ChargeError),
    ({'any_of': ['charged', 'uncharged']}, 'Fe', ChargeError),
    ({'any_of': ('charged', 'uncharged')}, 'Fe', ChargeError),

    ({'require': ['fermion', 'charged'],
      'any_of': ['lepton', 'baryon'],
      'exclude': ['antimatter']},
     'p+',
     None),

    ({'require': ['fermion', 'charged'],
      'any_of': ['lepton', 'baryon'],
      'exclude': ['antimatter']},
     'p+',
     None),

    ({'require': ['fermion', 'charged'],
      'any_of': ['lepton', 'baryon'],
      'exclude': ['matter']},
     'p+',
     AtomicError),
]


@pytest.mark.parametrize(
    "decorator_kwargs, particle, expected_exception",
    decorator_categories_table,
)
def test_decorator_categories(decorator_kwargs, particle, expected_exception):
    """Tests the require, any_of, and exclude categories lead to an
    AtomicError being raised when an inputted particle does not meet
    the required criteria, and do not lead to an AtomicError when the
    inputted particle matches the criteria."""

    @particle_input(**decorator_kwargs)
    def decorated_function(argument: Particle):
        return argument

    if expected_exception:
        with pytest.raises(expected_exception):
            decorated_function(particle)
    else:
        decorated_function(particle)


is_element = ['H', 'Fe-56', 'p+', 'alpha', 'Fe', 'D+', 'T 1-']
not_element = ['e-', 'e+', 'n', 'mu-', 'tau+']

is_isotope = ['D', 'T', 'alpha', 'proton', 'Fe-56', 'Be-8']
not_isotope = ['H', 'e-', 'n', 'p-', 'e+', 'Fe', 'Au', 'Og']

is_ion = ['p+', 'D+', 'T+', 'alpha', 'Be-8+', 'Fe 26+']
not_ion = ['D', 'T', 'H-1', 'He-4', 'e-', 'e+', 'n']


@particle_input
def function_with_element_argument(element: Particle):
    return element


@particle_input
def function_with_isotope_argument(isotope: Particle):
    return isotope


@particle_input
def function_with_ion_argument(ion: Particle):
    return ion


@pytest.mark.parametrize('element', is_element)
def test_is_element(element):
    function_with_element_argument(element)


@pytest.mark.parametrize('particle', not_element)
def test_not_element(particle):
    with pytest.raises(InvalidElementError):
        function_with_element_argument(particle)


@pytest.mark.parametrize('isotope', is_isotope)
def test_is_isotope(isotope):
    function_with_isotope_argument(isotope)


@pytest.mark.parametrize('particle', not_isotope)
def test_not_element(particle):
    with pytest.raises(InvalidIsotopeError):
        function_with_isotope_argument(particle)


@pytest.mark.parametrize('ion', is_ion)
def test_is_ion(ion):
    function_with_ion_argument(ion)


@pytest.mark.parametrize('particle', not_ion)
def test_not_ion(particle):
    with pytest.raises(InvalidIonError):
        function_with_ion_argument(particle)
