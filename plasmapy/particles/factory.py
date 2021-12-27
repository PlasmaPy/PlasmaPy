"""A factory function for creating particle objects."""

__all__ = ["physical_particle_factory"]

from typing import Union

from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.particles.particle_class import (
    CustomParticle,
    DimensionlessParticle,
    Particle,
)
from plasmapy.particles.particle_collections import ParticleList


def physical_particle_factory(
    *args, **kwargs
) -> Union[Particle, CustomParticle, ParticleList]:
    """
    Return a representation of one or more physical particles.

    This function will select the appropriate type among
    `~plasmapy.particles.particle_class.Particle`,
    `~plasmapy.particles.particle_class.CustomParticle`, and
    `~plasmapy.particle_collections.particle_class.ParticleList`.

    Parameters
    ----------
    *args
        Positional arguments to be supplied to the appropriate particle
        class.

    Z : array-like, keyword-only, optional
        The charge number of the particle.

    mass_numb : integer or array of integers, keyword-only, optional
        The mass number of an isotope.

    mass : `~astropy.units.Quantity`, keyword-only, optional
        The mass of the particle.

    charge : array-like, keyword-only, optional
        The electrical charge of the particle.

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

    Examples
    --------
    >>> from plasmapy.particles.factory import physical_particle_factory
    >>> import astropy.units as u
    >>> physical_particle_factory("p+")
    Particle("p+")
    >>> physical_particle_factory(mass = 9e-26 * u.kg, charge = 8e20 * u.C)
    CustomParticle(mass=9e-26 kg, charge=8e+20 C)
    >>> physical_particle_factory(["p+", "e-"])
    ParticleList(['p+', 'e-'])

    """
    if (
        len(args) == 1
        and not kwargs
        and isinstance(args[0], (Particle, CustomParticle, ParticleList))
    ):
        return args[0]

    if args and isinstance(args[0], DimensionlessParticle):
        raise TypeError("Unable to create a ")

    if not args and not kwargs:
        raise TypeError("No particle information has been provided.")

    for particle_type in (Particle, CustomParticle, ParticleList):
        try:
            return particle_type(*args, **kwargs)
        except (TypeError, InvalidParticleError):
            pass

    raise InvalidParticleError(
        f"Unable to create an appropriate particle object with "
        f"args={args} and kwargs={kwargs}."
    )
