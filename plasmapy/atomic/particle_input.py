import functools

from .particle_class import Particle
from ..utils import AtomicError


def particle_input(*particle_input_args, **particle_input_kwargs):
    r"""A decorator that takes inputs related to particles and passes
    through the corresponding instance of the Particle class."""
    def decorator(func, *decorator_input_args, **decorator_input_kwargs):
        @functools.wraps(func)
        def wrapper(*wrapper_args, **wrapper_kwargs):
            particle = Particle(*args, **kwargs)
            return func(particle)
        return wrapper
    return decorator
