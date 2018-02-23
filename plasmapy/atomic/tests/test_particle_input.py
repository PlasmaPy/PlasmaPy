import pytest
from typing import Optional

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
        a, particle: Particle, b=None, Z: int = None,
        mass_numb: int = None) -> Particle:
    """
    A simple function that, when decorated with `@particle_input`,
    returns the instance of the Particle class corresponding to the
    inputs.
    """
    if not isinstance(particle, Particle):
        raise TypeError(
            f"The argument particle in func_simple_noparens is not a Particle")
    return particle


@particle_input()
def func_simple_parens(
        a, particle: Particle, b=None, Z: int = None,
        mass_numb: int = None) -> Particle:
    """
    A simple function that, when decorated with `@particle_input()`,
    returns the instance of the Particle class corresponding to the
    inputs.
    """
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
    """
    Test that simple functions decorated by particle_input correctly
    return the correct Particle object.
    """
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
    """
    Test that functions decorated with particle_input raise the
    expected errors.
    """
    with pytest.raises(expected_error, message=(
            f"{func} did not raise {expected_error} with kwargs = {kwargs}")):
        func(**kwargs)


class Test_particle_input:
    """
    A sample class with methods that will be used to check that
    @particle_input works both with and without parentheses.
    """

    @particle_input
    def method_noparens(self, particle: Particle) -> Particle:
        return particle

    @particle_input()
    def method_parens(self, particle: Particle) -> Particle:
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


@particle_input
def function_with_no_annotations():
    """A trivial function that is incorrectly decorated with
    particle_input because no arguments are annotated with Particle."""
    pass


def test_no_annotations_exception():
    """Test that a function decorated with particle_input that has no
    annotated arguments will raise an AtomicError."""
    with pytest.raises(AtomicError):
        function_with_no_annotations()


@particle_input
def ambiguous_keywords(p1: Particle, p2: Particle, Z=None, mass_numb=None):
    """A trivial function with two annotated arguments plus the keyword
    arguments `Z` and `mass_numb`."""
    pass


def test_function_with_ambiguity():
    """Test that a function decorated with particle_input that has two
    annotated arguments"""
    with pytest.raises(AtomicError):
        ambiguous_keywords('H', 'He', Z=1, mass_numb=4)
    with pytest.raises(AtomicError):
        ambiguous_keywords('H', 'He', Z=1)
    with pytest.raises(AtomicError):
        ambiguous_keywords('H', 'He', mass_numb=4)
    # TODO: should particle_input raise an exception when Z and mass_numb
    # are given as keyword arguments but are not explicitly set?


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
    ({'require': 'charged'}, 'Fe 0+', ChargeError),

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
    def decorated_function(argument: Particle) -> Particle:
        return argument

    if expected_exception:
        with pytest.raises(expected_exception):
            decorated_function(particle)
    else:
        decorated_function(particle)


def test_none_shall_pass():
    """Tests the `none_shall_pass` keyword argument in is_particle.  If
    `none_shall_pass=True`, then an annotated argument should allow
    `None` to be passed through to the decorated function.  If
    `none_shall_pass=False`, then particle_input should raise a
    `~plasmapy.utils.AtomicError` if an annotated argument is assigned
    the value of `None`.
    """
    @particle_input(none_shall_pass=True)
    def func_none_shall_pass(particle: Particle) -> Optional[Particle]:
        return particle

    @particle_input(none_shall_pass=False)
    def func_none_shall_not_pass(particle: Particle) -> Particle:
        return particle

    assert func_none_shall_pass(None) is None, \
        ("The none_shall_pass keyword in the particle_input decorator is set "
         "to True, but is not passing through None.")

    with pytest.raises(TypeError, message=(
            "The none_shall_pass keyword in the particle_input decorator is "
            "set to False, but is not raising a TypeError.")):
        func_none_shall_not_pass(None)


# TODO: The following tests might be able to be cleaned up and/or
# further parametrized since there's a fair bit of repetition.

is_element = ['H', 'Fe-56', 'p+', 'alpha', 'Fe', 'D+', 'T 1-']
not_element = ['e-', 'e+', 'n', 'mu-', 'tau+']

is_isotope = ['D', 'T', 'alpha', 'proton', 'Fe-56', 'Be-8']
not_isotope = ['H', 'e-', 'n', 'p-', 'e+', 'Fe', 'Au', 'Og']

is_ion = ['p+', 'D+', 'T+', 'alpha', 'Be-8+', 'Fe 26+']
not_ion = ['D', 'T', 'H-1', 'He-4', 'e-', 'e+', 'n']


@particle_input
def function_with_element_argument(element: Particle) -> Particle:
    """A function decorated with `~plasmapy.atomic.particle_input`
    where the argument annotated with `~plasmapy.atomic.Particle`
    is named `element`.  This function should raise an
    `~plasmapy.utils.InvalidElementError` when the argument is not
    an element, isotope, or ion."""
    return element


@particle_input
def function_with_isotope_argument(isotope: Particle) -> Particle:
    """A function decorated with `~plasmapy.atomic.particle_input`
    where the argument annotated with `~plasmapy.atomic.Particle`
    is named `isotope`.  This function should raise an
    `~plasmapy.utils.InvalidIsotopeError` when the argument is not an
    isotope or an ion of an isotope."""
    return isotope


@particle_input
def function_with_ion_argument(ion: Particle) -> Particle:
    """
    A function decorated with `~plasmapy.atomic.particle_input`
    where the argument annotated with `~plasmapy.atomic.Particle`
    is named `ion`.  This function should raise an
    `~plasmapy.utils.InvalidIonError` when the argument is not an
    ion.
    """
    return ion


@pytest.mark.parametrize('element', is_element)
def test_is_element(element):
    """
    Test that particle_input will not raise an
    `~plasmapy.utils.InvalidElementError` if the annotated argument is
    named 'element' and is assigned values that are elements, isotopes,
    or ions.
    """
    function_with_element_argument(element)


@pytest.mark.parametrize('particle', not_element)
def test_not_element(particle):
    """
    Test that particle_input will raise an
    `~plasmapy.utils.InvalidElementError` if an argument decorated with
    `~plasmapy.atomic.Particle` is named 'element', but the annotated
    argument ends up not being an element, isotope, or ion.
    """
    with pytest.raises(InvalidElementError, message=(
            "@particle_input is not raising an InvalidElementError for "
            f"{repr(particle)} even though the annotated argument is named "
            "'element'.")):
        function_with_element_argument(particle)


@pytest.mark.parametrize('isotope', is_isotope)
def test_is_isotope(isotope):
    """
    Test that particle_input will not raise an
    `~plasmapy.utils.InvalidIsotopeError` if the annotated argument is
    named 'isotope' and is assigned values that are isotopes or
    ions of isotopes."""
    function_with_isotope_argument(isotope)


@pytest.mark.parametrize('particle', not_isotope)
def test_not_isotope(particle):
    """
    Test that particle_input will raise an
    `~plasmapy.utils.InvalidIsotopeError` if an argument decorated with
    `~plasmapy.atomic.Particle` is named 'isotope', but the annotated
    argument ends up not being an isotope or an ion of an isotope.
    """
    with pytest.raises(InvalidIsotopeError, message=(
            "@particle_input is not raising an InvalidIsotopeError for "
            f"{repr(particle)} even though the annotated argument is named "
            "'isotope'.")):
        function_with_isotope_argument(particle)


@pytest.mark.parametrize('ion', is_ion)
def test_is_ion(ion):
    """
    Test that particle_input will not raise an
    `~plasmapy.utils.InvalidIonError` if the annotated argument is
    named 'ion' and is assigned values that are ions.
    """
    function_with_ion_argument(ion)


@pytest.mark.parametrize('particle', not_ion)
def test_not_ion(particle):
    """
    Test that particle_input will raise an
    `~plasmapy.utils.InvalidIonError` if an argument decorated with
    `~plasmapy.atomic.Particle` is named 'ion', but the annotated
    argument ends up not being an ion.
    """
    with pytest.raises(InvalidIonError, message=(
            "@particle_input is not raising an InvalidIonError for "
            f"{repr(particle)} even though the annotated argument is named "
            "'ion'.")):
        function_with_ion_argument(particle)
