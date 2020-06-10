"""
A decorator to take string and/or integer representations of particles
as arguments and pass through the corresponding instance of the
`~plasmapy.particles.Particle` class.
"""

import functools
import inspect
import numbers
from typing import Any, Callable, List, Optional, Set, Tuple, Union

from plasmapy.particles.exceptions import (
    AtomicError,
    ChargeError,
    InvalidElementError,
    InvalidIonError,
    InvalidIsotopeError,
    InvalidParticleError,
)

from .particle_class import Particle

__all__ = ["particle_input"]


def _particle_errmsg(
    argname: str,
    argval: str,
    Z: int = None,
    mass_numb: int = None,
    funcname: str = None,
) -> str:
    """
    Return a string with an appropriate error message for an
    `~plasmapy.utils.InvalidParticleError`.
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
                f"The particle {phrase} of the following categories: " f"{condition}. "
            )

    return category_errmsg


def particle_input(
    wrapped_function: Callable = None,
    require: Union[str, Set, List, Tuple] = None,
    any_of: Union[str, Set, List, Tuple] = None,
    exclude: Union[str, Set, List, Tuple] = None,
    none_shall_pass: bool = False,
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
        all of these categories, then an `~plasmapy.utils.AtomicError`
        will be raised.

    any_of : `str`, `set`, `list`, or `tuple`, optional
        Categories that a particle may be in.  If a particle is not in
        any of these categories, then an `~plasmapy.utils.AtomicError`
        will be raised.

    exclude : `str`, `set`, `list`, or `tuple`, optional
        Categories that a particle cannot be in.  If a particle is in
        any of these categories, then an `~plasmapy.utils.AtomicError`
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
    then the decorator will raise an `~plasmapy.utils.InvalidElementError`,
    `~plasmapy.utils.InvalidIsotopeError`, or `~plasmapy.utils.InvalidIonError`
    if the particle does not correspond to an element, isotope, or ion,
    respectively.

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

    `~plasmapy/utils/InvalidParticleError`
        If the annotated argument does not correspond to a valid
        particle.

    `~plasmapy/utils/InvalidElementError`
        If an annotated argument is named ``element``, and the input
        does not correspond to an element, isotope, or ion.

    `~plasmapy/utils/InvalidIsotopeError`
        If an annotated argument is named ``isotope``, and the input
        does not correspond to an isotope or an ion of an isotope.

    `~plasmapy/utils/InvalidIonError`
        If an annotated argument is named ``ion``, and the input does
        not correspond to an ion.

    `~plasmapy/utils/ChargeError`
        If ``'charged'`` is in the ``require`` argument and the particle
        is not explicitly charged, or if ``any_of = {'charged',
        'uncharged'}`` and the particle does not have charge information
        associated with it.

    `~plasmapy/utils/AtomicError`
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
        wrapped_signature = inspect.signature(wrapped_function)

        # add '__signature__' to methods that are copied from
        # wrapped_function onto wrapper
        assigned = list(functools.WRAPPER_ASSIGNMENTS)
        assigned.append("__signature__")

        @functools.wraps(wrapped_function, assigned=assigned)
        def wrapper(*args, **kwargs):
            annotations = wrapped_function.__annotations__
            bound_args = wrapped_signature.bind(*args, **kwargs)

            default_arguments = bound_args.signature.parameters
            arguments = bound_args.arguments
            argnames = bound_args.signature.parameters.keys()

            # Handle optional-only arguments in function declaration
            for default_arg in default_arguments:
                # The argument is not contained in `arguments` if the
                # user does not explicitly pass an optional argument.
                # In such cases, manually add it to `arguments` with
                # the default value of parameter.
                if default_arg not in arguments:
                    arguments[default_arg] = default_arguments[default_arg].default

            funcname = wrapped_function.__name__

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
                raise AtomicError(
                    f"None of the arguments or keywords to {funcname} "
                    f"have been annotated with Particle, as required "
                    f"by the @particle_input decorator."
                )
            elif len(args_to_become_particles) > 1:
                if "Z" in argnames or "mass_numb" in argnames:
                    raise AtomicError(
                        f"The arguments Z and mass_numb in {funcname} are not "
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

            for argname in argnames:
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
                            argname, params, already_particle, funcname
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

        # Convert the argument to a Particle object if it is not
        # already one.

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
            raise AtomicError(
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
