"""
A module containing an interface function that accepts inputs intended
for |Particle|, |CustomParticle|, or |ParticleList| and returns the
appropriate instance of one of those three classes.
"""

__all__ = []

import astropy.units as u
import contextlib

from astropy.units.physical import electrical_charge, mass
from numbers import Integral
from typing import Union

from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.particles.particle_class import CustomParticle, Particle
from plasmapy.particles.particle_collections import ParticleList


def _physical_particle_factory(
    *args, **kwargs
) -> Union[Particle, CustomParticle, ParticleList]:
    """
    Return a representation of one or more physical particles.

    This function will select the appropriate type among |Particle|,
    |CustomParticle|, and |ParticleList|.

    .. caution::

       If |Quantity| instances are provided to this function as
       positional arguments, then they must presently be in the order
       expected by |CustomParticle|.

    Parameters
    ----------
    *args
        Positional arguments to be supplied to |Particle|,
        |CustomParticle|, or |ParticleList|.

    **kwargs
        Keyword arguments to be supplied to |Particle|,
        |CustomParticle|, or |ParticleList|.

    Returns
    -------
    |Particle|, |CustomParticle|, or |ParticleList|

    Raises
    ------
    `InvalidParticleError`
        If an appropriate particle could not be constructed.

    `TypeError`
        If no positional arguments and no keyword arguments were
        provided.

    See Also
    --------
    ~plasmapy.particles.particle_class.Particle
    ~plasmapy.particles.particle_class.CustomParticle
    ~plasmapy.particles.particle_class.ParticleList

    Notes
    -=---
    If ``Z`` or ``mass_numb`` is provided as keyword arguments but are
    equal to `None`, then they will not be provided to any of the calls.
    This is to allow |CustomParticle| instances to be created.

    Examples
    --------
    >>> from plasmapy.particles._factory import _physical_particle_factory
    >>> import astropy.units as u
    >>> _physical_particle_factory("p+")
    Particle("p+")
    >>> _physical_particle_factory(mass = 9e-26 * u.kg, charge = 8e20 * u.C)
    CustomParticle(mass=9e-26 kg, charge=8e+20 C)
    >>> _physical_particle_factory(["p+", "e-"])
    ParticleList(['p+', 'e-'])
    """

    # We need to remove Z and mass_numb from kwargs when they are `None`
    # because they are not allowed as arguments to `CustomParticle`, and
    # are not needed in kwargs if they are their default values. Note
    # that this affects `not kwargs` below.
    for parameter in ["Z", "mass_numb"]:
        if parameter in kwargs and kwargs[parameter] is None:
            kwargs.pop(parameter)

    if (
        len(args) == 1
        and not kwargs
        and isinstance(args[0], (Particle, CustomParticle, ParticleList))
    ):
        return args[0]

    if not args and not kwargs:
        raise TypeError("Particle information has not been provided.")

    for particle_type in (Particle, CustomParticle, ParticleList):
        with contextlib.suppress(TypeError, InvalidParticleError):
            return particle_type(*args, **kwargs)

    if isinstance(args[0], u.Quantity):
        physical_type = u.get_physical_type(args[0])
        if physical_type not in (electrical_charge, mass):
            raise u.UnitConversionError(
                "Cannot create a particle object with a Quantity with a "
                f"physical type of {physical_type}."
            )

    if not isinstance(args[0], (str, Integral, CustomParticle, Particle, ParticleList)):
        raise TypeError("Invalid type for particle.")

    raise InvalidParticleError(
        f"Unable to create an appropriate particle object with "
        f"args={args} and kwargs={kwargs}."
    )
