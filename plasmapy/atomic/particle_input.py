import functools
import inspect

from .particle_class import Particle
from ..utils import (AtomicError,
                     InvalidParticleError,
                     InvalidElementError,
                     InvalidIonError,
                     InvalidIsotopeError,
                     ChargeError)


from typing import Callable, Union, Any, Set, List, Tuple

# TODO: make sure particle_input works with classes and class methods





def particle_input(wrapped_function: Callable = None,
                   must_be: Union[str, Set, List, Tuple] = set(),
                   cannot_be: Union[str, Set, List, Tuple] = set(),
                   any: bool = False,
                   **kwargs) -> Any:
    r"""A decorator to take arguments and keywords related to a
    particle, and pass through the Particle class to the callable instead.

    Parameters
    ----------

    wrapped_function : callable
        The function to be decorated.

    must_be : str, set, list, or tuple; optional
        The categories that the Particle must be in; otherwise the


    Notes
    -----

    This version of particle_input only works with functions, but in the future
    should be extended to work with classes and methods.

    """

    def _errmsg(argname: str, argval: str, Z: int = None,
                mass_numb: int = None, funcname: str = None) -> str:
        r"""Returns a string with an appropriate error message for an
        InvalidParticleError."""

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

    def decorator(wrapped_function: Callable):
        wrapped_signature = inspect.signature(wrapped_function)

        @functools.wraps(wrapped_function)
        def wrapper(*args, **kwargs):

            # The following couple of statements will likely need to be
            # modified in order to work with methods.

            bound_args = wrapped_signature.bind(*args, **kwargs)
            arguments = bound_args.arguments
            argnames = bound_args.arguments.keys()
            funcname = wrapped_function.__name__
            annotations = wrapped_function.__annotations__

            args_to_become_particles = [
                argname for argname in annotations.keys()
                if annotations[argname] is Particle
            ]

            if not args_to_become_particles:
                raise AtomicError(
                    f"None of the arguments or keywords to {funcname} "
                    f"have been annotated with Particle as required "
                    f"by the particle_input decorator.")
            elif len(args_to_become_particles) > 1:
                if 'Z' in argnames or 'mass_numb' in argnames:
                    raise AtomicError(
                        f"The arguments Z and mass_numb in {funcname} are not "
                        f"allowed when more than one argument or keyword is "
                        f"annotated with Particle in functions decorated "
                        f"with @particle_input.")

            # If the number of arguments and keywords annotated with
            # Particle is exactly one, then the Z and mass_numb keywords
            # can be used without potential for ambiguity.

            Z = arguments.get('Z', None)
            mass_numb = arguments.get('mass_numb', None)

            if Z is not None and not isinstance(Z, int):
                raise TypeError(
                    f"The argument Z = {repr(Z)} in {funcname} is not an "
                    f"integer.")

            if mass_numb is not None and not isinstance(mass_numb, int):
                raise TypeError(
                    f"The argument mass_numb = {repr(mass_numb)} in "
                    f"{funcname} is not an integer.")

            # Go through the argument names and check whether or not they are
            # annotated with Particle.  If they aren't, include the name and
            # value of the argument as an item in the new keyword arguments
            # dictionary unchanged.  If they are annotated with Particle, then
            # either convert the representation of a Particle to a Particle if
            # it is not already a Particle.

            new_kwargs = {}

            for argname in argnames:
                argval = arguments[argname]

                should_be_particle = argname in args_to_become_particles
                not_already_particle = not isinstance(argval, Particle)

                if should_be_particle and not_already_particle:

                    if not isinstance(argval, (int, str)):
                        raise TypeError

                    try:
                        particle = Particle(argval, Z=Z, mass_numb=mass_numb)
                    except InvalidParticleError as e:
                        raise InvalidParticleError(
                            _errmsg(argname, argval, Z, mass_numb, funcname)
                        ) from e

                    if argname == 'element' and not particle.element:
                        raise InvalidElementError
                    if argname == 'isotope' and not particle.isotope:
                        raise InvalidIsotopeError
                    if argname == 'ion' and not particle.ion:
                        raise InvalidIonError
                    if 'charged' in must_be and not particle._integer_charge:
                        raise ChargeError
                    if not particle.is_category(must_be, exclude=cannot_be):
                        raise AtomicError


                    new_kwargs[argname] = particle
                else:
                    new_kwargs[argname] = argval

            return wrapped_function(**new_kwargs)

        return wrapper

    # The following code helps allow the decorator to be used either with
    # or without arguments.  In particular, this helps allow us to invoke
    # the decorator either as `@particle_input` or as `@particle_input()`,
    # where the latter call allows the decorator to have keyword arguments.

    if wrapped_function is not None:
        return decorator(wrapped_function)
    else:
        return decorator
