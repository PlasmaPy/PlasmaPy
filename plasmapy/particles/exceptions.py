"""
Collection of `Exceptions` and `Warnings` for PlasmaPy particles.
"""
__all__ = [
    "ParticleError",
    "ParticleTypeError",
    "ParticleValueError",
    "ParticleWarning",
    "ChargeError",
    "InvalidElementError",
    "InvalidIonError",
    "InvalidIsotopeError",
    "InvalidParticleError",
    "MissingParticleDataError",
    "MissingParticleDataWarning",
    "UnexpectedParticleError",
]

from plasmapy.utils import PlasmaPyError, PlasmaPyWarning


class ParticleError(PlasmaPyError):
    """Base exception for errors in the `~plasmapy.particles` subpackage."""

    pass


class ParticleTypeError(ParticleError, TypeError):
    """
    An exception for errors in the `~plasmapy.particles` subpackage
    caused because an incorrect type was provided.
    """

    pass


class ParticleValueError(ParticleError, ValueError):
    """
    An exception for errors in the `~plasmapy.particles` subpackage
    caused because an incorrect value was provided.
    """

    pass


class MissingParticleDataError(ParticleError):
    """An exception for missing atomic or particle data."""

    pass


class ChargeError(ParticleError):
    """An exception for incorrect or missing charge information."""

    pass


class UnexpectedParticleError(ParticleError):
    """An exception for when a particle is not of the expected category."""

    pass


class InvalidIonError(UnexpectedParticleError):
    """
    An exception for when an argument is a valid particle but not a
    valid ion.
    """

    pass


class InvalidIsotopeError(UnexpectedParticleError):
    """
    An exception for when an argument is a valid particle but not a
    valid isotope.
    """

    pass


class InvalidElementError(UnexpectedParticleError):
    """
    An exception for when an argument is a valid particle is not a
    valid element.
    """

    pass


class InvalidParticleError(ParticleError):
    """An exception for when a particle is invalid."""

    pass


class ParticleWarning(PlasmaPyWarning):
    """The base warning for the `~plasmapy.particles` subpackage."""

    pass


class MissingParticleDataWarning(ParticleWarning):
    """Warning for use when atomic or particle data is missing."""

    pass
