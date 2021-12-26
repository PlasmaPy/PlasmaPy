"""A factory function for creating particle objects."""

import astropy.units as u

from numbers import Integral, Real
from typing import Optional, Union

from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.particles.particle_class import CustomParticle, Particle
from plasmapy.particles.particle_collections import ParticleList


def particle(
    *args,
    Z: Optional[Real] = None,
    mass_numb: Optional[Integral] = None,
    mass: Optional[Union[u.Quantity, Real]] = None,
    charge: Optional[Union[u.Quantity, Real]] = None,
):
    """
    Return a representation of a physical particle.

    Parameters
    ----------
    Z : `~numbers.Real`, keyword-only, optional

    mass_numb : `~numbers.Real`, keyword-only, optional

    mass : `~astropy.units.Quantity` or `~numbers.Real`, keyword-only, optional

    charge : `~astropy.units.Quantity` or `~numbers.Real`, keyword-only, optional

    Raises
    ------
    `InvalidParticleError`

    """
    particle_types = [Particle, CustomParticle, ParticleList]

    # Should DimensionlessParticle be included?

    atomic_info_in_kwargs = Z is not None or mass_numb is not None
    custom_info_in_kwargs = mass is not None or charge is not None

    if atomic_info_in_kwargs and custom_info_in_kwargs:
        raise InvalidParticleError("Cannot provide ")

    no_kwargs = not atomic_info_in_kwargs and not custom_info_in_kwargs

    # Return the original object if it is already an instance of one of
    # the appropriate particle classes.
    if all([no_kwargs, len(args) == 1, isinstance(args[0], particle_types)]):
        return args[0]

    if atomic_info_in_kwargs and custom_info_in_kwargs:
        raise InvalidParticleError("...")

    if not args and no_kwargs:
        return CustomParticle(mass=np.nan * u.kg, charge=np.nan * u.C)
