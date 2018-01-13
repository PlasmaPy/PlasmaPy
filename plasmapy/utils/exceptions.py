"""
Custom Error and Warning names to improve readability
"""


# ----------
# Exceptions:
# ----------

class PlasmaPyError(Exception):
    r"""
    Base class of PlasmaPy custom errors.

    All custom exceptions raised by PlasmaPy should inherit from this
    class and be defined in this module.

    Custom exceptions can inherit from other exception types too.
    Thus, if code already knows how to handle a ValueError, it won't
    need any specific modification.
    """
    pass


class PhysicsError(PlasmaPyError, ValueError):
    r"""Error for use of a physics value outside PlasmaPy theoretical
    bounds"""
    pass


class RelativityError(PhysicsError):
    r"""Error for use of a speed greater than or equal to the speed of
    light"""
    pass


class AtomicError(PlasmaPyError):
    r"""An exception for errors occurring in the atomic subpackage."""
    pass


class MissingAtomicDataError(AtomicError):
    r"""An exception for when atomic data is missing."""
    pass


class ChargeError(AtomicError):
    r"""An exception for when charge information is incorrect or
    missing."""
    pass


class UnexpectedParticleError(AtomicError):
    r"""An exception for when a particle is not of the expected
    category."""
    pass


class InvalidIonError(UnexpectedParticleError):
    r"""An exception for when an argument is not an ion."""
    pass


class InvalidIsotopeError(UnexpectedParticleError):
    r"""An exception for when an argument is not an isotope."""
    pass


class InvalidElementError(UnexpectedParticleError):
    r"""An exception for when an argument is not an element."""
    pass


class InvalidParticleError(AtomicError):
    r"""An exception for when a particle is invalid."""
    pass


# ----------
# Warnings:
# ----------

class PlasmaPyWarning(Warning):
    r"""Base class of PlasmaPy custom warnings.

    All PlasmaPy custom warnings should inherit from this class and be defined
    in this module.

    Warnings should be issued using warnings.warn, which will not break
    execution if unhandled.
    """
    pass


class PhysicsWarning(PlasmaPyWarning):
    r"""Warning for using a mildly worrisome physics value"""
    pass


class RelativityWarning(PhysicsWarning):
    r"""Warning for use of a speed quantity approaching the speed of light"""
    pass


class AtomicWarning(PlasmaPyWarning):
    r"""Warnings for use in the atomic subpackage."""
    pass


class MissingAtomicDataWarning(AtomicWarning):
    r"""Warning for use when atomic data is missing."""
    pass
