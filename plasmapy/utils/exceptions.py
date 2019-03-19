"""Exceptions and warnings specific to PlasmaPy."""


# ----------
# Exceptions
# ----------

class PhysicsError(Exception, ValueError):
    """
    The base exception for physics-related errors.
    """
    pass


class RelativityError(PhysicsError):
    """
    An exception for speeds greater than the speed of light.
    """
    pass


class AtomicError(Exception):
    """An exception for errors in the `~plasmapy.atomic` subpackage."""
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


class DataStandardError(Exception):
    """An exception for when HDF5 is not defined by OpenPMD standard."""
    pass


# ----------
# Warnings:
# ----------

class PhysicsWarning(Warning):
    """The base warning for `~plasmapy.physics` related warnings."""
    pass


class RelativityWarning(PhysicsWarning):
    """
    A warning for when relativistic velocities are being used in or are
    returned by non-relativistic functionality.
    """
    pass


class CouplingWarning(PhysicsWarning):
    """
    A warning for functions that rely on a particular coupling regime to be valid.
    """


class AtomicWarning(Warning):
    """The base warning for the `~plasmapy.atomic` subpackage."""
    pass


class MissingAtomicDataWarning(AtomicWarning):
    """Warning for use when atomic or particle data is missing."""
    pass
