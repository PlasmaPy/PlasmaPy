"""Typing constructs for `plasmapy`."""

from __future__ import annotations

import astropy.units as u

from numbers import Integral
from typing import Iterable, Union

from plasmapy.particles.particle_class import CustomParticle, Particle
from plasmapy.particles.particle_collections import ParticleList

RealParticleLike = Union[str, Integral, Particle]
CustomParticleLike = Union[u.Quantity, CustomParticle]
SingleParticleLike = Union[RealParticleLike, CustomParticleLike]
ParticleListLike = Union[Iterable[SingleParticleLike], ParticleList]
ParticleLike = Union[SingleParticleLike, ParticleListLike]


def _get_docstring(types: str):
    """
    Produce a standardized docstring for typing constructs for
    `plasmapy.particles`.
    """
    return f"""
A typing construct for the creation of {types} objects.

For more information, please refer to |particle_like|.
"""


RealParticleLike.__doc__ = _get_docstring("|Particle|")
CustomParticleLike.__doc__ = _get_docstring("|CustomParticle|")
SingleParticleLike.__doc__ = _get_docstring("|Particle| or |CustomParticle|")
ParticleListLike.__doc__ = _get_docstring("ParticleList")
ParticleLike.__doc__ = _get_docstring("|Particle|, |CustomParticle|, or |ParticleList|")
