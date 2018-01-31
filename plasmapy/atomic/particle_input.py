import functools
import inspect

from .particle_class import Particle
from ..utils import AtomicError

# TODO: make sure particle_input works with classes and class methods


def particle_input(wrapped_function=None, **kwargs):
    r"""A decorator to take arguments and keywords related to a
    particle, and pass through the Particle class to the callable instead.

    Notes
    -----

    This version of particle_input only works with functions, but in the future
    should be extended to work with classes and methods.

    """

    def decorator(wrapped_function):
        wrapped_signature = inspect.signature(wrapped_function)

        @functools.wraps(wrapped_function)
        def wrapper(*args, **kwargs):

            # The following couple of statements will likely need to be
            # modified in order to work with classes.

            bound_args = wrapped_signature.bind(*args, **kwargs)
            arguments = bound_args.arguments
            argnames = bound_args.arguments.keys()

            annotations = wrapped_function.__annotations__

            args_to_become_particles = [
                argname for argname in annotations.keys()
                if annotations[argname] is Particle
            ]

            if not args_to_become_particles:
                raise AtomicError(
                    f"None of the arguments or keywords to {func.__name__} "
                    f"have been annotated with the Particle class as required "
                    f"by the particle_input decorator. For example,\n\n"
                    f"    def func(arg, particle: Particle, kwarg=None):\n\n"
                    f"will convert the particle argument to ")

            # If the number of arguments and keywords annotated with
            # Particle is exactly one, then the Z and mass_numb keywords
            # can be used without potential for ambiguity.

            if len(args_to_become_particles) == 1:
                Z = arguments.get('Z', None)
                mass_numb = arguments.get('mass_numb', None)
            else:
                if 'Z' in argnames:
                    raise AtomicError
                elif 'mass_numb' in argnames:
                    raise AtomicError
                else:
                    Z = None
                    mass_numb = None

            # Go through the argument names and check whether or not they are
            # annotated with Particle.  If they aren't, include the name and
            # value of the argument as an item in the new keyword arguments
            # dictionary unchanged.  If they are annotated with Particle, then
            # either convert the representation of a Particle to a Particle if
            # it is not already a Particle.

            new_kwargs = {}

            for argname in argnames:
                argval = arguments[argname]

                should_be_particle = argname in annotations.keys()
                not_already_particle = not isinstance(argval, Particle)

                if should_be_particle and not_already_particle:
                    new_kwargs[argname] = \
                        Particle(argval, Z=Z, mass_numb=mass_numb)
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
