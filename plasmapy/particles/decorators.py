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
import numbers
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


def original_particle_input(
    wrapped_function: Callable = None,
    require: Union[str, Set, List, Tuple] = None,
    any_of: Union[str, Set, List, Tuple] = None,
    exclude: Union[str, Set, List, Tuple] = None,
    none_shall_pass=False,
) -> Any:
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

    none_shall_pass : `bool`, optional
        If set to `True`, then the decorated argument may be set to
        `None` without raising an exception.  In such cases, this
        decorator will pass `None` through to the decorated function or
        method.  If set to `False` and the annotated argument is given
        a value of `None`, then this decorator will raise a `TypeError`.

    Notes
    -----
    If the annotated argument is named `element`, `isotope`, or `ion`,
    then the decorator will raise an
    `~plasmapy.particles.exceptions.InvalidElementError`,
    `~plasmapy.particles.exceptions.InvalidIsotopeError`, or
    `~plasmapy.particles.exceptions.InvalidIonError` if the particle
    does not correspond to an element, isotope, or ion, respectively.

    If exactly one argument is annotated with `~plasmapy.particles.Particle`,
    then the keywords ``Z`` and ``mass_numb`` may be used to specify the
    integer charge and/or mass number of an ion or isotope.  However,
    the decorated function must allow ``Z`` and/or ``mass_numb`` as keywords
    in order to enable this functionality.

    Raises
    ------
    `TypeError`
        If the annotated argument is not a `str`, `int`, `tuple`, `list`
        or `~plasmapy.particles.Particle`; or if ``Z`` or ``mass_numb`` is
        not an `int`.

    `ValueError`
        If the number of input elements in a collection do not match the
        number of expected elements.

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
        If an annotated argument is named ``ion``, and the input does
        not correspond to an ion.

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

    if exclude is None:
        exclude = set()
    if any_of is None:
        any_of = set()
    if require is None:
        require = set()

    def decorator(wrapped_function: Callable):
        wrapped_signature = inspect.signature(wrapped_function)  # DONE

        # add '__signature__' to methods that are copied from
        # wrapped_function onto wrapper
        assigned = list(functools.WRAPPER_ASSIGNMENTS)  # DONE
        assigned.append("__signature__")  # DONE

        @functools.wraps(wrapped_function, assigned=assigned)
        def wrapper(*func_args, **func_kwargs):
            annotations = wrapped_function.__annotations__  # DONE
            bound_args = wrapped_signature.bind(*func_args, **func_kwargs)

            default_arguments = bound_args.signature.parameters
            arguments = bound_args.arguments
            argument_names = bound_args.signature.parameters.keys()

            # Handle optional-only arguments in function declaration
            # Handled in bound_args setter
            for default_arg in default_arguments:
                # The argument is not contained in `arguments` if the
                # user does not explicitly pass an optional argument.
                # In such cases, manually add it to `arguments` with
                # the default value of parameter.
                if default_arg not in arguments:
                    arguments[default_arg] = default_arguments[default_arg].default

            function_name = wrapped_function.__name__

            args_to_become_particles = []
            for argname in annotations.keys():
                if isinstance(annotations[argname], tuple):
                    if argname == "return":
                        continue
                    annotated_argnames = annotations[argname]
                    expected_params = len(annotated_argnames)
                    received_params = len(arguments[argname])
                    if not expected_params == received_params:
                        raise ValueError(
                            f"Number of parameters allowed in the tuple "
                            f"({expected_params} parameters) are "
                            f"not equal to number of parameters passed in "
                            f"the tuple ({received_params} parameters)."
                        )
                elif isinstance(annotations[argname], list):
                    annotated_argnames = annotations[argname]
                    expected_params = len(annotated_argnames)
                    if expected_params > 1:
                        raise TypeError(
                            f"Put in [Particle] as the annotation to "
                            f"accept arbitrary number of Particle arguments."
                        )
                else:
                    annotated_argnames = (annotations[argname],)

                for annotated_argname in annotated_argnames:
                    is_particle = (
                        annotated_argname is Particle
                        or annotated_argname is Optional[Particle]
                    )
                    if is_particle and argname != "return":
                        args_to_become_particles.append(argname)

            if not args_to_become_particles:
                raise ParticleError(
                    f"None of the arguments or keywords to {function_name} "
                    f"have been annotated with Particle, as required "
                    f"by the @particle_input decorator."
                )
            elif len(args_to_become_particles) > 1:
                if "Z" in argument_names or "mass_numb" in argument_names:
                    raise ParticleError(
                        f"The arguments Z and mass_numb in {function_name} are not "
                        f"allowed when more than one argument or keyword is "
                        f"annotated with Particle in functions decorated "
                        f"with @particle_input."
                    )

            for x in args_to_become_particles:
                if (
                    annotations[x] is Particle
                    and isinstance(arguments[x], (tuple, list))
                    and len(arguments[x]) > 1
                ):
                    raise TypeError(
                        f"You cannot pass a tuple or list containing "
                        f"Particles when only single Particle was "
                        f"expected, instead found {arguments[x]}. If you "
                        f"intend to pass more than 1 Particle instance, "
                        f"use a tuple or a list type. "
                        f"That is use (Particle, Particle, ...) or "
                        f"[Particle] in function declaration."
                    )

            # If the number of arguments and keywords annotated with
            # Particle is exactly one, then the Z and mass_numb keywords
            # can be used without potential for ambiguity.

            Z = arguments.get("Z", None)
            mass_numb = arguments.get("mass_numb", None)

            # Go through the argument names and check whether or not they are
            # annotated with Particle.  If they aren't, include the name and
            # value of the argument as an item in the new keyword arguments
            # dictionary unchanged.  If they are annotated with Particle, then
            # either convert the representation of a Particle to a Particle if
            # it is not already a Particle and then do error checks.

            new_kwargs = {}

            for argname in argument_names:
                raw_argval = arguments[argname]
                if isinstance(raw_argval, (tuple, list)):
                    # Input argument value is a tuple or list
                    # of correspoding particles or atomic values.
                    argval_tuple = raw_argval
                    particles = []
                else:
                    # Otherwise convert it to tuple anyway so it can work
                    # with loops too.
                    argval_tuple = (raw_argval,)

                for pos, argval in enumerate(argval_tuple):
                    should_be_particle = argname in args_to_become_particles
                    already_particle = isinstance(argval, Particle)

                    # If the argument is not annotated with Particle, then we just
                    # pass it through to the new keywords without doing anything.

                    if not should_be_particle:
                        new_kwargs[argname] = raw_argval
                        continue

                    # Occasionally there will be functions where it will be
                    # useful to allow None as an argument.

                    # In case annotations[argname] is a collection (which looks
                    # like (Particle, Optional[Particle], ...) or [Particle])
                    if isinstance(annotations[argname], tuple):
                        optional_particle = (
                            annotations[argname][pos] is Optional[Particle]
                        )
                    elif isinstance(annotations[argname], list):
                        optional_particle = annotations[argname] == [Optional[Particle]]
                    else:
                        # Otherwise annotations[argname] must be a Particle itself
                        optional_particle = annotations[argname] is Optional[Particle]

                    if (optional_particle or none_shall_pass) and argval is None:
                        particle = None
                    else:
                        params = (argval, Z, mass_numb)
                        particle = get_particle(
                            argname, params, already_particle, function_name
                        )

                    if isinstance(raw_argval, (tuple, list)):
                        # If passed argument is a tuple or list, keep
                        # appending them.
                        particles.append(particle)
                        # Set appended values if current iteration is the
                        # last iteration.
                        if (pos + 1) == len(argval_tuple):
                            new_kwargs[argname] = tuple(particles)
                            del particles
                    else:
                        # Otherwise directly set values
                        new_kwargs[argname] = particle

            return wrapped_function(**new_kwargs)

        # add '__signature__' if it does not exist
        # - this will preserve parameter hints in IDE's
        if not hasattr(wrapper, "__signature__"):
            wrapper.__signature__ = inspect.signature(wrapped_function)

        return wrapper

    def get_particle(argname, params, already_particle, funcname):
        argval, Z, mass_numb = params
        """
        Convert the argument to a `~plasmapy.particles.Particle` object
        if it is not already one.
        """

        if not already_particle:

            if not isinstance(argval, (numbers.Integral, str, tuple, list)):
                raise TypeError(
                    f"The argument {argname} to {funcname} must be "
                    f"a string, an integer or a tuple or list of them "
                    f"corresponding to an atomic number, or a "
                    f"Particle object."
                )

            try:
                particle = Particle(argval, Z=Z, mass_numb=mass_numb)
            except InvalidParticleError as e:
                raise InvalidParticleError(
                    _particle_errmsg(argname, argval, Z, mass_numb, funcname)
                ) from e

        # We will need to do the same error checks whether or not the
        # argument is already an instance of the Particle class.

        if already_particle:
            particle = argval

        # If the name of the argument annotated with Particle in the
        # decorated function is element, isotope, or ion; then this
        # decorator should raise the appropriate exception when the
        # particle ends up not being an element, isotope, or ion.

        cat_table = [
            ("element", particle.element, InvalidElementError),
            ("isotope", particle.isotope, InvalidIsotopeError),
            ("ion", particle.ionic_symbol, InvalidIonError),
        ]

        for category_name, category_symbol, CategoryError in cat_table:
            if argname == category_name and not category_symbol:
                raise CategoryError(
                    f"The argument {argname} = {repr(argval)} to "
                    f"{funcname} does not correspond to a valid "
                    f"{argname}."
                )

        # Some functions require that particles be charged, or
        # at least that particles have charge information.

        _integer_charge = particle._attributes["integer charge"]

        must_be_charged = "charged" in require
        must_have_charge_info = set(any_of) == {"charged", "uncharged"}

        uncharged = _integer_charge == 0
        lacks_charge_info = _integer_charge is None

        if must_be_charged and (uncharged or must_have_charge_info):
            raise ChargeError(f"A charged particle is required for {funcname}.")

        if must_have_charge_info and lacks_charge_info:
            raise ChargeError(f"Charge information is required for {funcname}.")

        # Some functions require particles that belong to more complex
        # classification schemes.  Again, be sure to provide a
        # maximally useful error message.

        if not particle.is_category(require=require, exclude=exclude, any_of=any_of):
            raise ParticleError(
                _category_errmsg(particle, require, exclude, any_of, funcname)
            )

        return particle

    # The following code allows the decorator to be used either with or
    # without arguments.  This allows us to invoke the decorator either
    # as `@particle_input` or as `@particle_input()`, where the latter
    # call allows the decorator to have keyword arguments.

    if wrapped_function is not None:
        return decorator(wrapped_function)
    else:
        return decorator


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

        none_shall_pass : `bool`, optional
            If set to `True`, then the decorated argument may be set to
            `None` without raising an exception.  In such cases, this
            decorator will pass `None` through to the decorated function or
            method.  If set to `False` and the annotated argument is given
            a value of `None`, then this decorator will raise a `TypeError`.

        Notes
        -----
        If the annotated argument is named `element`, `isotope`, or `ion`,
        then the decorator will raise an
        `~plasmapy.particles.exceptions.InvalidElementError`,
        `~plasmapy.particles.exceptions.InvalidIsotopeError`, or
        `~plasmapy.particles.exceptions.InvalidIonError` if the particle
        does not correspond to an element, isotope, or ion, respectively.

        If exactly one argument is annotated with `~plasmapy.particles.Particle`,
        then the keywords ``Z`` and ``mass_numb`` may be used to specify the
        integer charge and/or mass number of an ion or isotope.  However,
        the decorated function must allow ``Z`` and/or ``mass_numb`` as keywords
        in order to enable this functionality.

        Raises
        ------
        `TypeError`
            If the annotated argument is not a `str`, `int`, `tuple`, `list`
            or `~plasmapy.particles.Particle`; or if ``Z`` or ``mass_numb`` is
            not an `int`.

        `ValueError`
            If the number of input elements in a collection do not match the
            number of expected elements.

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
            If an annotated argument is named ``ion``, and the input does
            not correspond to an ion.

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
