from plasmapy.utils import PlasmaPyError, PlasmaPyWarning


class AtomicError(PlasmaPyError):
    """An exception for errors in the `~plasmapy.particles` subpackage."""

    pass


class MissingAtomicDataError(AtomicError):
    """An exception for missing atomic or particle data."""

    pass


class ChargeError(AtomicError):
    """An exception for incorrect or missing charge information."""

    pass


class UnexpectedParticleError(AtomicError):
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


class InvalidParticleError(AtomicError):
    """An exception for when a particle is invalid."""

    pass


class AtomicWarning(PlasmaPyWarning):
    """The base warning for the `~plasmapy.particles` subpackage."""

    pass


class MissingAtomicDataWarning(AtomicWarning):
    """Warning for use when atomic or particle data is missing."""

    pass
