import functools

from .particle_class import Particle
from ..utils import AtomicError


class particle_input:

    def __init__(self, func):
        self.func = func

    def __call__(self, *args, **kwargs):
        return self.func(Particle(*args, **kwargs))
