"""Collection of exceptions and warnings for `plasmapy.particles`."""

__all__ = [
    "ParticleError",
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


class MissingParticleDataError(ParticleError):
    """
    An exception for missing atomic or particle data in the
    `~plasmapy.particles` subpackage.
    """

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
