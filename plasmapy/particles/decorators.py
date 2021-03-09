"""
Module used to define the framework needed for the `particle_input` decorator.
The decorator takes string and/or integer representations of particles
as arguments and passes through the corresponding instance of the
`~plasmapy.particles.Particle` class.
"""
__all__ = ["particle_input"]

import astropy.constants as const
import astropy.units as u
import collections
import functools
import inspect
import numpy as np

from typing import (
    AbstractSet,
    Any,
    Callable,
    Dict,
    List,
    Mapping,
    NoReturn,
    Optional,
    Set,
    Tuple,
    Union,
)

from plasmapy.particles.exceptions import (
    ChargeError,
    InvalidElementError,
    InvalidIonError,
    InvalidIsotopeError,
    InvalidParticleError,
    ParticleError,
    UnexpectedParticleError,
)
from plasmapy.particles.particle_class import CustomParticle, Particle, ParticleLike
from plasmapy.particles.particle_collections import ParticleList


def _particle_errmsg(
    argname: str,
    argval: str,
    Z: int = None,
    mass_numb: int = None,
    funcname: str = None,
) -> str:
    """
    Return a string with an appropriate error message for an
    `~plasmapy.particles.exceptions.InvalidParticleError`.
    """
    errmsg = f"In {funcname}, {argname} = {repr(argval)} "
    if mass_numb is not None or Z is not None:
        errmsg += "with "
    if mass_numb is not None:
        errmsg += f"mass_numb = {repr(mass_numb)} "
    if mass_numb is not None and Z is not None:
        errmsg += "and "
    if Z is not None:
        errmsg += f"integer charge Z = {repr(Z)} "
    errmsg += "does not correspond to a valid particle."
    return errmsg


def _category_errmsg(particle, require, exclude, any_of, funcname) -> str:
    """
    Return an appropriate error message for when a particle does not
    meet the required categorical specifications.
    """
    category_errmsg = (
        f"The particle {particle} does not meet the required "
        f"classification criteria to be a valid input to {funcname}. "
    )

    errmsg_table = [
        (require, "must belong to all"),
        (any_of, "must belong to any"),
        (exclude, "cannot belong to any"),
    ]

    for condition, phrase in errmsg_table:
        if condition:
            category_errmsg += (
                f"The particle {phrase} of the following categories: {condition}. "
            )

    return category_errmsg


class _ParticleInput:
    @classmethod
    def as_decorator(cls, func=None, **kwargs):
        """
        Convert arguments to methods and functions to
        `~plasmapy.particles.Particle` objects.

        Take positional and keyword arguments that are annotated with
        `~plasmapy.particles.Particle`, and pass through the
        `~plasmapy.particles.Particle` object corresponding to those arguments
        to the decorated function or method.

        Optionally, raise an exception if the particle does not satisfy the
        specified categorical criteria.

        Parameters
        ----------
        wrapped_function : `callable`
            The function or method to be decorated.

        require : `str`, `set`, `list`, or `tuple`, optional
            Categories that a particle must be in.  If a particle is not in
            all of these categories, then an `~plasmapy.particles.exceptions.ParticleError`
            will be raised.

        any_of : `str`, `set`, `list`, or `tuple`, optional
            Categories that a particle may be in.  If a particle is not in
            any of these categories, then an `~plasmapy.particles.exceptions.ParticleError`
            will be raised.

        exclude : `str`, `set`, `list`, or `tuple`, optional
            Categories that a particle cannot be in.  If a particle is in
            any of these categories, then an `~plasmapy.particles.exceptions.ParticleError`
            will be raised.

        allow_particle_lists : `bool`
            If `True`, then this decorator will allow
            `~plasmapy.particles.ParticleList` instances to be created
            and/or passed through.

        allow_custom_particles : `bool`
            If `True`, then this decorator will allow
            `~plasmapy.particles.ParticleList` instances to be created
            and/or passed through.

        Notes
        -----
        If the annotated argument is named ``element``, ``isotope``,
        ``ion``, or ``ionic_level``, then the decorator will check if
        the `plasmapy.particles.ParticleLike` object matches the
        corresponding category.  The argument ``ion`` requires that the
        particle be charged (e.g., ``Particle("He 1+")``, whereas
        ``ionic_level`` allows neutral atoms (e.g., ``Particle("He 0+")``)


        If the annotated argument is named ``element``, ``isotope``, or
        ``ion``, then the decorator will raise an
        `~plasmapy.particles.exceptions.InvalidElementError`,
        `~plasmapy.particles.exceptions.InvalidIsotopeError`, or
        `~plasmapy.particles.exceptions.InvalidIonError` if the particle
        does not correspond to an element, isotope, or ion, respectively.

        If exactly one argument is annotated with `~plasmapy.particles.Particle`,
        then the keywords ``Z`` and ``mass_numb`` may be used to specify the
        integer charge and/or mass number of an ion or isotope.

        Raises
        ------
        `TypeError`
            If the annotated argument is not a `str`, `int`, `tuple`, `list`
            or `~plasmapy.particles.Particle`; or if ``Z`` or ``mass_numb`` is
            not an `int`.

        `ValueError`
            If the number of input elements in a collection does not
            match the number of expected elements.

        `~plasmapy.particles.exceptions.InvalidParticleError`
            If the annotated argument does not correspond to a valid
            particle.

        `~plasmapy.particles.exceptions.InvalidElementError`
            If an annotated argument is named ``element``, and the input
            does not correspond to an element, isotope, or ion.

        `~plasmapy.particles.exceptions.InvalidIsotopeError`
            If an annotated argument is named ``isotope``, and the input
            does not correspond to an isotope or an ion of an isotope.

        `~plasmapy.particles.exceptions.InvalidIonError`
            If an annotated argument is named ``ion`` (or ``ionic_level``),
            and the input does not correspond to an ion (or ionic level).

        `~plasmapy.particles.exceptions.ChargeError`
            If ``'charged'`` is in the ``require`` argument and the particle
            is not explicitly charged, or if ``any_of = {'charged',
            'uncharged'}`` and the particle does not have charge information
            associated with it.

        `~plasmapy.particles.exceptions.ParticleError`
            If an annotated argument does not meet the criteria set by the
            categories in the ``require``, ``any_of``, and ``exclude``
            keywords; if more than one argument is annotated and ``Z`` or
            ``mass_numb`` are used as arguments; or if none of the arguments
            have been annotated with `~plasmapy.particles.Particle`.

        Examples
        --------
        The following simple decorated function returns the
        `~plasmapy.particles.Particle` object created from the function's
        sole argument:

        .. code-block:: python

            from plasmapy.particles import particle_input, Particle
            @particle_input
            def decorated_function(particle: Particle):
                return particle

        This decorator may also be used to accept arguments using tuple
        annotation containing specific number of elements or using list
        annotation which accepts any number of elements in an iterable.
        Returns a tuple of `~plasmapy.particles.Particle`:

        .. code-block:: python

            from plasmapy.particles import particle_input, Particle
            @particle_input
            def decorated_tuple_function(particles: (Particle, Particle)):
                return particles
            sample_particles = decorated_tuple_function(('He', 'Li'))

            @particle_input
            def decorated_list_function(particles: [Particle]):
                return particles
            sample_particles = decorated_list_function(('Al 3+', 'C'))
            sample_particles = decorated_list_function(['He', 'Ne', 'Ar'])

        This decorator may be used for methods in instances of classes, as
        in the following example:

        .. code-block:: python

            from plasmapy.particles import particle_input, Particle
            class SampleClass:
                @particle_input
                def decorated_method(self, particle: Particle):
                    return particle
            sample_instance = SampleClass()
            sample_instance.decorated_method('Fe')

        Some functions may intended to be used with only certain categories
        of particles.  The ``require``, ``any_of``, and ``exclude`` keyword
        arguments enable this functionality.

        .. code-block:: python

            from plasmapy.particles import particle_input, Particle
            @particle_input(
                require={'matter'},
                any_of={'charged', 'uncharged},
                exclude={'neutrino', 'antineutrino'},
            )
            def selective_function(particle: Particle):
                return particle
        """
        # The following code allows the decorator to be used either with
        # or without arguments.  We can invoke the decorator either as
        # `@particle_input` or `@particle_input()`, where the latter
        # form allows the decorator to have arguments.
        self = cls(**kwargs)
        if func is not None and not kwargs:
            return self(func)
        else:
            return self

    def __init__(
        self,
        func=None,
        require=None,
        any_of=None,
        exclude=None,
        allow_particle_lists=True,
        allow_custom_particles=True,
    ):
        self._data = collections.defaultdict(lambda: None)
        self._data["new_kwargs"] = dict()
        self.require = require if require else set()
        self.any_of = any_of if any_of else set()
        self.exclude = exclude if exclude else set()
        self.allow_particle_lists = bool(allow_particle_lists)
        self.allow_custom_particles = bool(allow_custom_particles)

    def __call__(self, wrapped_function: Callable):

        self.wrapped_function = wrapped_function

        assigned = list(functools.WRAPPER_ASSIGNMENTS)
        assigned.append("__signature__")

        @functools.wraps(wrapped_function, assigned=assigned)
        def wrapper(*func_args, **func_kwargs):
            self.func_args = func_args
            self.func_kwargs = func_kwargs
            self.process_arguments()
            return_ = wrapped_function(**self.new_kwargs)
            return return_

        return wrapper

    @property
    def wrapped_function(self) -> Optional[Callable]:
        """The function that is being decorated."""
        return self._data["wrapped_function"]

    @wrapped_function.setter
    def wrapped_function(self, function: Callable):
        self._data["wrapped_function"] = function
        self._data["wrapped_signature"] = inspect.signature(function)
        self._data["annotations"] = getattr(function, "__annotations__", {}).copy()
        if "return" in self._data["annotations"]:
            del self._data["annotations"]["return"]

    @property
    def annotations(self) -> Optional[Dict[str, Any]]:
        """
        Annotations for arguments in the decorated function.

        The keys are the arguments to the decorated function that have
        annotations.  The associated values are the annotations themselves.
        Only arguments with annotations are included.
        """
        return self._data["annotations"]

    @property
    def _wrapped_signature(self) -> Optional[inspect.Signature]:
        return self._data["wrapped_signature"]

    @property
    def _bound_args(self) -> Optional[inspect.BoundArguments]:
        bound_args_ = self._wrapped_signature.bind(*self.func_args, **self.func_kwargs)
        bound_args_.apply_defaults()
        return bound_args_

    @property
    def _default_arguments(self) -> Mapping[str, inspect.Parameter]:
        return self._bound_args.signature.parameters

    @property
    def new_kwargs(self) -> Optional[Dict[str, Any]]:
        """The revised keyword arguments to be supplied to the decorated function."""
        return self._data["new_kwargs"]

    @property
    def func_args(self) -> Optional[tuple]:
        """Positional arguments as originally supplied to the decorated function."""
        return self._data["func_args"]

    @func_args.setter
    def func_args(self, value: tuple):
        self._data["func_args"] = value

    @property
    def func_kwargs(self) -> Optional[Dict[str, Any]]:
        """Keyword arguments as originally supplied to the decorated function."""
        return self._data["func_kwargs"]

    @func_kwargs.setter
    def func_kwargs(self, value: Dict[str, Any]):
        self._data["func_kwargs"] = value

    @property
    def original_values(self) -> Dict[str, Any]:
        """
        The values that were originally passed as positional and keyword
        arguments to the decorated function.
        """
        return self._bound_args.arguments

    @property
    def argument_names(self) -> AbstractSet:
        """The names of the arguments to the original function."""
        return self._bound_args.signature.parameters.keys()

    @property
    def particle_like_annotations(self) -> tuple:
        """Annotations that correspond to `ParticleLike`."""
        return ParticleLike, Optional[ParticleLike]

    @property
    def particle_list_annotations(self) -> tuple:
        """Annotations that correspond to a collection of particles."""
        return ParticleList, Optional[ParticleList]

    @property
    def all_particle_annotations(self) -> tuple:
        """All annotations to be processed"""
        return self.particle_list_annotations + self.particle_like_annotations

    @property
    def optional_particle_annotations(self) -> tuple:
        """Annotations indicating that the value of the argument can be `None`."""
        return Optional[ParticleLike], Optional[ParticleList]

    def has_particle_like_annotation(self, argument) -> bool:
        """
        Return `True` if ``argument`` has an annotation indicating that
        """
        return self.annotations[argument] in self.particle_like_annotations

    def has_particle_list_annotation(self, argument) -> bool:
        annotation = self.annotations[argument]
        if annotation in self.particle_list_annotations:
            return True
        if hasattr(annotation, "__len__"):
            return all(item in self.particle_like_annotations for item in annotation)
        return False

    def particle_like_arguments(self) -> List[Any]:
        """
        The names of arguments that should be processed into particle
        objects.
        """
        return_ = []
        for argument, annotation in self.annotations.items():
            if self.has_particle_like_annotation(argument):
                return_.append(argument)
        return return_

    def particle_list_arguments(self) -> List[Any]:
        """
        The names of arguments that should be processed into
        `ParticleList` objects.
        """
        return_ = []
        for argument, annotation in self.annotations.items():
            if annotation in self.particle_list_annotations:
                return_.append(argument)
        return return_

    @property
    def new_kwargs(self) -> Optional[Dict[str, Any]]:
        """
        The processed keyword arguments that are being provided to the
        decorated function.
        """
        return self._data["new_kwargs"]

    @property
    def arguments_that_should_not_go_in_new_args(self):
        """
        Arguments that are used to create a `Particle` but should not be
        passed to the original function.
        """
        return "mass_numb", "Z"

    def put_original_argument_into_new_args(self, argument: str) -> NoReturn:
        """
        Pass through the value for the original argument to the
        dictionary of new arguments without changing it.
        """
        if argument in self.new_kwargs.keys():
            raise RuntimeError(f"{argument} already in new_kwargs")
        self.new_kwargs[argument] = self.original_values[argument]

    def put_processed_argument_into_new_args(self, argument: str, value) -> NoReturn:
        """
        Pass through a revised value for the original argument to the
        dictionary of new arguments.
        """
        assert argument not in self.new_kwargs.keys()
        self.new_kwargs[argument] = value

    def check_for_errors(self):
        """Check for errors in the arguments that were provided."""
        Z = self.original_values.get("Z", None)
        mass_numb = self.original_values.get("mass_numb", None)

        Z_or_mass_numb_provided = Z is not None or mass_numb is not None
        if len(self.particle_like_arguments()) != 1 and Z_or_mass_numb_provided:
            raise ParticleError

    def process_arguments(self) -> NoReturn:
        """
        Go through the arguments and convert them to the appropriate
        particle or particle collection.
        """
        self.new_kwargs.clear()
        self.check_for_errors()
        for argument in self.original_values.keys():
            if argument in self.arguments_that_should_not_go_in_new_args:
                pass
            elif argument not in self.annotations:
                self.put_original_argument_into_new_args(argument)
            else:
                self.process_annotated_argument(argument)

    def process_annotated_argument(self, argument: str):
        """Process an argument that has an annotation."""
        original_value = self.original_values[argument]
        annotation = self.annotations[argument]

        none_shall_pass = annotation in self.optional_particle_annotations

        if original_value is None and none_shall_pass:
            self.put_original_argument_into_new_args(argument)
        elif self.has_particle_like_annotation(argument):
            self.process_particle_like_argument(argument)
        elif self.has_particle_list_annotation(argument):
            self.process_particle_list_argument(argument)
        else:
            self.put_original_argument_into_new_args(argument)

    def process_particle_like_argument(self, argument: str) -> NoReturn:
        original_value = self.original_values[argument]
        if isinstance(original_value, u.Quantity):
            self.process_quantity_into_custom_particle(argument)
        elif isinstance(original_value, (str, int, Particle)):
            self.process_particle_like_into_particle(argument)

    def process_quantity_into_custom_particle(self, argument: str) -> NoReturn:

        if not self.allow_custom_particles:
            raise ParticleError("Custom particles are not allowed.")

        original_value = self.original_values[argument]

        try:
            n = len(original_value)
        except TypeError:
            pass
        else:
            if n != 1:
                raise NotImplementedError(
                    "Only non-array quantities can be provided at this time."
                )

        is_mass = original_value.unit.physical_type == "mass"
        is_charge = original_value.unit.physical_type == "electrical charge"

        Z = getattr(self.original_values, "Z", None)
        mass_numb = getattr(self.original_values, "Z", None)

        if is_charge and Z is not None:
            raise ParticleError(
                "Redundant or contradictory charge information was provided."
            )

        if mass_numb is not None:
            raise ParticleError(
                "The argument mass_numb cannot be provided for a particle "
                "represented by a mass or a charge."
            )

        mass = original_value if is_mass else np.nan * u.kg

        if Z is not None:
            charge = Z * const.e
        elif is_charge:
            charge = original_value
        else:
            charge = np.nan * u.C

        custom_particle = CustomParticle(mass=mass, charge=charge)
        self.put_processed_argument_into_new_args(argument, custom_particle)

    def process_particle_like_into_particle(self, argument):
        original_value = self.original_values[argument]
        mass_numb = self.original_values.get("mass_numb", None)
        Z = self.original_values.get("Z", None)

        if isinstance(original_value, Particle) and Z is None and mass_numb is None:
            new_particle = original_value
        else:
            new_particle = Particle(original_value, mass_numb=mass_numb, Z=Z)

        if not new_particle.is_category(
            require=self.require, any_of=self.any_of, exclude=self.exclude
        ):
            errmsg = _category_errmsg(
                new_particle,
                require=self.require,
                exclude=self.exclude,
                any_of=self.any_of,
                funcname=self.wrapped_function.__name__,
            )
            exception = ParticleError
            #            if self.any_of.issuperset({"charged", "uncharged"}):
            if {"charged", "uncharged"}.issubset(self.any_of):
                if not new_particle.is_category(any_of={"charged", "uncharged"}):
                    exception = ChargeError
            raise exception(errmsg)

        self.check_special_particle_names(argument, new_particle)
        self.put_processed_argument_into_new_args(argument, new_particle)

    def check_special_particle_names(self, argument, particle):
        funcname = self.wrapped_function.__name__
        if argument == "element" and not particle.is_category("element"):
            errmsg = _category_errmsg(particle, "element", {}, {}, funcname)
            raise InvalidElementError(errmsg)
        elif argument == "isotope" and not particle.is_category("isotope"):
            errmsg = _category_errmsg(particle, "isotope", {}, {}, funcname)
            raise InvalidIsotopeError(errmsg)
        elif argument == "ion" and not particle.is_category("ion"):
            errmsg = _category_errmsg(particle, "ion", {}, {}, funcname)
            raise InvalidIonError(errmsg)
        elif argument == "ionic_level" and not particle.is_category("ion"):
            if not particle.is_category(require="element"):
                raise InvalidElementError
            elif not particle.is_category(any_of={"charged", "uncharged"}):
                raise ChargeError

    def process_particle_list_argument(self, argument):
        original_value = self.original_values[argument]
        annotation = self.annotations[argument]

        if isinstance(original_value, ParticleList):
            particle_list = original_value
        elif not hasattr(original_value, "__len__"):
            particle_list = ParticleList((original_value))
        else:
            particle_list = ParticleList(original_value)

        try:
            expected_number_of_particles = len(annotation)
        except TypeError:
            pass
        else:
            actual_number_of_particles = len(particle_list)
            if actual_number_of_particles != expected_number_of_particles:
                raise ValueError(
                    f"Expecting {expected_number_of_particles} particles, "
                    f"but {actual_number_of_particles} were provided."
                )

        self.put_processed_argument_into_new_args(argument, particle_list)


particle_input = _ParticleInput.as_decorator
