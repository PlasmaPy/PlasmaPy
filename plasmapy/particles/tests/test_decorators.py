import astropy.units as u
import inspect
import pytest
import sys

from typing import List, Optional, Tuple, Union

from plasmapy.particles import ParticleList
from plasmapy.particles.decorators import particle_input
from plasmapy.particles.exceptions import (
    ChargeError,
    InvalidElementError,
    InvalidIonError,
    InvalidIsotopeError,
    InvalidParticleError,
    ParticleError,
)
from plasmapy.particles.particle_class import CustomParticle, Particle, ParticleLike
from plasmapy.utils.code_repr import call_string
from plasmapy.utils.decorators.validators import validate_quantities


@particle_input
def func_simple_noparens(
    a, particle: ParticleLike, b=None, Z: int = None, mass_numb: int = None
) -> Particle:
    """
    A simple function that, when decorated with `@particle_input`,
    returns the instance of the Particle class corresponding to the
    inputs.
    """
    if not isinstance(particle, Particle):
        raise TypeError(
            "The argument `particle` in func_simple_noparens should be "
            f"a Particle. Instead, {particle=}."
        )
    return particle


@particle_input()
def func_simple_parens(
    a, particle: ParticleLike, b=None, Z: int = None, mass_numb: int = None
) -> Particle:
    """
    A simple function that, when decorated with `@particle_input()`,
    returns the instance of the Particle class corresponding to the
    inputs.
    """
    if not isinstance(particle, Particle):
        raise TypeError(
            "The argument `particle` in func_simple_parens should be "
            f"a Particle. Instead, {particle=}."
        )
    return particle


particle_input_simple_table = [
    ((1, "p+"), {"b": 2}, "p+"),
    ((1, "Fe"), {"mass_numb": 56, "Z": 3}, "Fe-56 3+"),
    ((1,), {"particle": "e-"}, "e-"),
    ((1,), {"particle": Particle("e-")}, "e-"),
]


@pytest.mark.parametrize("func", [func_simple_parens, func_simple_noparens])
@pytest.mark.parametrize("args, kwargs, symbol", particle_input_simple_table)
def test_particle_input_simple(func, args, kwargs, symbol):
    """
    Test that simple functions decorated by particle_input correctly
    return the correct Particle object.
    """
    expected = Particle(symbol)
    result = func(*args, **kwargs)
    assert result == expected, (
        f"The result {result!r} does not equal the expected value of "
        f"{expected!r} with {func =}, {args =}, and {kwargs =}."
    )


# function, kwargs, expected_error
particle_input_error_table = [
    (func_simple_noparens, {"a": 1, "particle": "asdf"}, InvalidParticleError)
]


@pytest.mark.parametrize("func, kwargs, expected_error", particle_input_error_table)
def test_particle_input_errors(func, kwargs, expected_error):
    """
    Test that functions decorated with `@particle_input` raise the
    expected errors.
    """
    with pytest.raises(expected_error):
        func(**kwargs)
        pytest.fail(f"{func} did not raise {expected_error} with kwargs = {kwargs}")


class ClassWithDecoratedMethods:
    """
    A sample class with methods that will be used to check that
    @particle_input works both with and without parentheses.
    """

    @particle_input
    def method_decorated_with_no_parentheses(self, particle: ParticleLike) -> Particle:
        return particle

    @particle_input()
    def method_decorated_with_parentheses(self, particle: ParticleLike) -> Particle:
        return particle


@pytest.mark.parametrize("symbol", ["muon", "He 2+"])
def test_particle_input_classes(symbol):
    """Test that @particle_input works for instance methods."""
    expected = Particle(symbol)
    instance = ClassWithDecoratedMethods()

    result_noparens = instance.method_decorated_with_no_parentheses(symbol)
    result_parens = instance.method_decorated_with_parentheses(symbol)

    assert isinstance(result_noparens, Particle)
    assert isinstance(result_parens, Particle)

    assert result_parens == result_noparens == expected


@particle_input
def function_with_no_annotations():
    """
    A trivial function that is incorrectly decorated with @particle_input
    because no arguments are annotated with Particle.
    """
    pass


def test_no_annotations_exception():
    """
    Test that a function decorated with @particle_input that has no
    annotated arguments will raise a ParticleError.
    """
    with pytest.raises(ParticleError):
        function_with_no_annotations()


@particle_input
def ambiguous_keywords(p1: ParticleLike, p2: ParticleLike, Z=None, mass_numb=None):
    """
    A trivial function with two annotated arguments plus the keyword
    arguments `Z` and `mass_numb`.
    """
    pass


ambiguous_arguments = [
    [("H", "He"), {"Z": 1, "mass_numb": 4}],
    [("H", "He"), {"Z": 1}],
    [("H", "He"), {"mass_numb": 4}],
]


@pytest.mark.parametrize("args, kwargs", ambiguous_arguments)
def test_function_with_ambiguity(args, kwargs):
    """
    Test that a function decorated with particle_input that has two
    annotated arguments.
    """
    with pytest.raises(ParticleError):
        ambiguous_keywords(*args, **kwargs)


def function_to_test_annotations(particles: Union[Tuple, List], resulting_particles):
    """
    Test that a function with an argument annotated with (Particle,
    Particle, ...) or [Particle] returns a tuple of expected Particle
    instances.

    .. deprecated::

        This capability has been deprecated.

    Arguments
    =========
    particles: tuple or list
        A collection containing many items, each of which may be a valid
        representation of a particle or a |Particle| instance.
    """

    expected = [
        particle if isinstance(particle, Particle) else Particle(particle)
        for particle in particles
    ]

    returned_particle_instances = all(
        isinstance(p, Particle) for p in resulting_particles
    )

    returned_correct_instances = all(
        expected[i] == resulting_particles[i] for i in range(len(particles))
    )

    if not returned_particle_instances:
        raise ParticleError(
            f"A function decorated by particle_input did not return "
            f"a collection of Particle instances for input of "
            f"{particles!r}, and instead returned"
            f"{resulting_particles!r}."
        )

    if not returned_correct_instances:
        raise ParticleError(
            f"A function decorated by particle_input did not return "
            f"{expected!r} as expected, and instead returned "
            f"{resulting_particles!r}."
        )


class TestOptionalArgs:
    def particle_iter(self, particles):
        return [Particle(particle) for particle in particles]

    def test_optional_particle(self):
        particle = "He"

        @particle_input
        def has_default_particle(particle: ParticleLike = particle):
            return particle

        assert has_default_particle() == Particle(particle)
        assert has_default_particle("Ne") == Particle("Ne")


def test_nonhashable_annotation():
    """
    Verify that the annotation for a different parameter can be a
    non-hashable object such as a `list`.

    Notes
    -----
    This problem arose when the collections containing the annotations
    to be processed were sets, and there was an operation like
    ``[] in set()``.
    """

    @particle_input
    def has_nonhashable_annotation(x: [], particle: ParticleLike):
        pass

    has_nonhashable_annotation(1, "e+")


# decorator_kwargs, particle, expected_exception
decorator_categories_table = [
    ({"exclude": {"element"}}, "Fe", ParticleError),
    ({"any_of": {"lepton", "antilepton"}}, "tau-", None),
    ({"require": {"isotope", "ion"}}, "Fe-56+", None),
    ({"require": {"isotope", "ion"}}, "Fe+", ParticleError),
    ({"any_of": {"isotope", "ion"}}, "Fe+", None),
    ({"any_of": {"charged", "uncharged"}}, "Fe", ChargeError),
    ({"any_of": ["charged", "uncharged"]}, "Fe", ChargeError),
    ({"any_of": ("charged", "uncharged")}, "Fe", ChargeError),
    ({"require": "charged"}, "Fe 0+", ChargeError),
    (
        {
            "require": ["fermion", "charged"],
            "any_of": ["lepton", "baryon"],
            "exclude": ["antimatter"],
        },
        "p+",
        None,
    ),
    (
        {
            "require": ["fermion", "charged"],
            "any_of": ["lepton", "baryon"],
            "exclude": ["antimatter"],
        },
        "p+",
        None,
    ),
    (
        {
            "require": ["fermion", "charged"],
            "any_of": ["lepton", "baryon"],
            "exclude": ["matter"],
        },
        "p+",
        ParticleError,
    ),
]


@pytest.mark.parametrize(
    "decorator_kwargs, particle, expected_exception", decorator_categories_table
)
def test_decorator_categories(decorator_kwargs, particle, expected_exception):
    """Tests the require, any_of, and exclude categories lead to an
    ParticleError being raised when an inputted particle does not meet
    the required criteria, and do not lead to an ParticleError when the
    inputted particle matches the criteria."""

    @particle_input(**decorator_kwargs)
    def decorated_function(argument: ParticleLike) -> Particle:
        return argument

    if expected_exception:
        with pytest.raises(expected_exception):
            decorated_function(particle)
            pytest.fail(
                f"{call_string(decorated_function, [], decorator_kwargs)} "
                f"did not raise {expected_exception}"
            )
    else:
        decorated_function(particle)


def test_optional_particle_annotation_argname():
    """Tests the `Optional[Particle]` annotation argument in a function
    decorated by `@particle_input` such that the annotated argument allows
    `None` to be passed through to the decorated function."""

    @particle_input
    def func_optional_particle(particle: Optional[ParticleLike]) -> Optional[Particle]:
        return particle

    assert func_optional_particle(None) is None, (
        "The particle keyword in the particle_input decorator is set "
        "to accept Optional[Particle], but is not passing through None."
    )


def undecorated_function(particle: ParticleLike, distance: u.m):
    return particle, distance


# Both particle_input and validate_quantities can be used with or
# without arguments, so the following list is to test the cases where

decorator_pairs = [
    (particle_input, validate_quantities),
    (particle_input(), validate_quantities),
    (particle_input, validate_quantities()),
    (particle_input(), validate_quantities()),
]


@pytest.mark.parametrize("decorator1, decorator2", decorator_pairs)
def test_stacking_decorators(decorator1, decorator2):
    """
    Test that particle_input and validate_quantities behave as expected.
    """
    decorated_function_1_2 = decorator1(decorator2(undecorated_function))
    decorated_function_2_1 = decorator2(decorator1(undecorated_function))

    particle_1_2, distance_1_2 = decorated_function_1_2("p+", distance=3 * u.cm)
    particle_2_1, distance_2_1 = decorated_function_2_1("p+", distance=3 * u.cm)

    # Test that particle_input is working as expected
    assert isinstance(particle_1_2, Particle)
    assert isinstance(particle_2_1, Particle)
    assert particle_1_2 == particle_2_1 == "p+"

    # Test that validate_quantities is working as expected
    assert distance_1_2.unit == distance_2_1.unit == u.m
    assert distance_1_2 == distance_2_1 == 3 * u.cm


@pytest.mark.parametrize("decorator1, decorator2", decorator_pairs)
def test_preserving_signature_with_stacked_decorators(decorator1, decorator2):
    """
    Test that particle_input & validate_quantities preserve the function
    signature after being stacked.
    """
    decorated_function_1_2 = decorator1(decorator2(undecorated_function))
    decorated_function_2_1 = decorator2(decorator1(undecorated_function))

    undecorated_signature = inspect.signature(undecorated_function)
    decorated_signature_1_2 = inspect.signature(decorated_function_1_2)
    decorated_signature_2_1 = inspect.signature(decorated_function_2_1)

    assert undecorated_signature == decorated_signature_1_2 == decorated_signature_2_1


@pytest.mark.xfail(
    condition=sys.version_info < (3, 9),
    reason="This test fails for Python 3.8 but it is not clear why.",
)
def test_annotated_classmethod():
    """
    Test that `particle_input` behaves as expected for a method that is
    decorated with `classmethod`.
    """

    class HasAnnotatedClassMethod:
        @classmethod
        @particle_input
        def f(cls, particle: ParticleLike):
            return particle

    has_annotated_classmethod = HasAnnotatedClassMethod()
    assert has_annotated_classmethod.f("p+") == Particle("p+")


@pytest.mark.parametrize(
    "inner_decorator, outer_decorator",
    [
        (particle_input, particle_input),
        (particle_input, particle_input()),
        (particle_input(), particle_input),
        (particle_input(), particle_input()),
    ],
)
def test_self_stacked_decorator(outer_decorator, inner_decorator):
    """Test that particle_input can be stacked with itself."""

    @outer_decorator
    @inner_decorator
    def f(x, particle: ParticleLike):
        return particle

    result = f(1, "p+")
    assert result == "p+"
    assert isinstance(result, Particle)


@pytest.mark.parametrize(
    "inner_decorator, outer_decorator",
    [
        (particle_input, particle_input),
        (particle_input, particle_input()),
        (particle_input(), particle_input),
        (particle_input(), particle_input()),
    ],
)
def test_class_stacked_decorator(outer_decorator, inner_decorator):
    """
    Test that particle_input can be stacked with itself for an
    instance method.
    """

    class Sample:
        @outer_decorator
        @inner_decorator
        def __init__(self, particle: ParticleLike):
            self.particle = particle

    result = Sample("p+")
    assert result.particle == "p+"
    assert isinstance(result.particle, Particle)


validate_quantities_ = validate_quantities(
    T_e={"equivalencies": u.temperature_energy()}
)


def test_annotated_init():
    """Test that `particle_input` can decorate an __init__ method."""

    class HasAnnotatedInit:
        @particle_input(require="element")
        def __init__(self, particle: ParticleLike, ionic_fractions=None):
            self.particle = particle

    x = HasAnnotatedInit("H-1", ionic_fractions=32)
    assert x.particle == "H-1"


@pytest.mark.parametrize(
    "outer_decorator, inner_decorator",
    [
        (particle_input, validate_quantities_),
        (particle_input(), validate_quantities_),
        pytest.param(validate_quantities_, particle_input, marks=pytest.mark.xfail),
        pytest.param(validate_quantities_, particle_input(), marks=pytest.mark.xfail),
    ],
)
def test_particle_input_with_validate_quantities(outer_decorator, inner_decorator):
    """Test that particle_input can be stacked with validate_quantities."""

    class C:
        @outer_decorator
        @inner_decorator
        def __init__(
            self,
            particle: ParticleLike,
            T_e: u.K = None,
        ):
            self.particle = particle
            self.T_e = T_e

    instance = C("p+", T_e=3.8 * u.eV)

    assert instance.particle == "p+"
    assert isinstance(instance.particle, Particle)

    assert instance.T_e.unit == u.K


rename_this = [
    ({"allow_custom_particles": False}, CustomParticle()),
    ({"allow_particle_lists": False}, ParticleList(["p+", "Fe"])),
    (
        {"allow_custom_particles": False, "allow_particle_lists": True},
        ParticleList([CustomParticle()]),
    ),
]


@pytest.mark.parametrize("kwargs_to_particle_input, arg", rename_this)
def test_particle_input_verification(kwargs_to_particle_input, arg):
    """Test the allow_custom_particles keyword argument to particle_input."""

    @particle_input(**kwargs_to_particle_input)
    def f(particle: ParticleLike):
        return particle

    with pytest.raises(ParticleError):
        f(arg)


class ParameterNamesCase:
    def __init__(
        self,
        category,
        function,
        particles_in_category,
        particles_not_in_category,
        exception,
    ):
        self.category = category
        self.function = function
        self.particles_in_category = particles_in_category
        self.particles_not_in_category = particles_not_in_category
        self.particles_some_in_category = (
            particles_in_category + particles_not_in_category
        )
        self.exception = exception


@particle_input
def get_element(element: ParticleLike):
    return element


@particle_input
def get_isotope(isotope: ParticleLike):
    return isotope


@particle_input
def get_ion(ion: ParticleLike):
    return ion


cases = [
    ParameterNamesCase(
        category="element",
        function=get_element,
        particles_in_category=["H", "Fe-56", "p+", "alpha", "Fe", "D+", "T 1-"],
        particles_not_in_category=["e-", "e+", "n", "mu-", "tau+"],
        exception=InvalidElementError,
    ),
    ParameterNamesCase(
        category="isotope",
        function=get_isotope,
        particles_in_category=["D", "T", "alpha", "proton", "Fe-56", "Be-8"],
        particles_not_in_category=["H", "e-", "n", "p-", "e+", "Fe", "Au", "Og"],
        exception=InvalidIsotopeError,
    ),
    ParameterNamesCase(
        category="ion",
        function=get_ion,
        particles_in_category=["p+", "D+", "T+", "alpha", "Be-8+", "Fe 26+"],
        particles_not_in_category=["D", "T", "H-1", "He-4", "e-", "e+", "n"],
        exception=InvalidIonError,
    ),
    ParameterNamesCase(
        category="ionic_level",
        function=get_ion,
        particles_in_category=["p+", "D+", "T+", "alpha", "Be-8+", "Fe 26+"],
        particles_not_in_category=["D", "T", "H-1", "He-4", "e-", "e+", "n"],
        exception=InvalidIonError,
    ),
]


@pytest.mark.parametrize("case", cases)
class TestParticleInputParameterNames:
    """
    Test the behavior associated with annotated special parameter names
    such as ``element``, ``isotope``, ``ion``, and ``ionic_level``. In
    particular, make sure that the resulting particle(s) belong to the
    expected categories.
    """

    def test_individual_particles_not_in_category(self, case):
        """
        Test that the appropriate exception is raised when the function
        is provided with individual particles that are not in the
        category.
        """
        for particle in case.particles_not_in_category:
            with pytest.raises(case.exception):
                case.function(particle)

    def test_particle_list_not_in_category(self, case):
        """
        Test that the appropriate exception is raised when the function
        is provided with multiple particles at once that are all not in
        the category.
        """
        with pytest.raises(case.exception):
            case.function(case.particles_not_in_category)

    def test_particle_list_some_in_category(self, case):
        """
        Test that the appropriate exception is raised when the function
        is provided with multiple particles at once, of which some are
        in the category and some are not.
        """
        with pytest.raises(case.exception):
            case.function(case.particles_some_in_category)

    # If "ionic_level" gets added as a particle category, then add
    # assertions to the following test using is_category.

    def test_individual_particles_all_in_category(self, case):
        """
        Test that no exception is raised when the function is provided
        with individual particles that are all in the category.
        """
        for particle in case.particles_in_category:
            case.function(particle)

    def test_particle_list_all_in_category(self, case):
        """
        Test that no exception is raised when the function is provided
        with multiple particles at once which are all in the category.
        """
        case.function(case.particles_in_category)
