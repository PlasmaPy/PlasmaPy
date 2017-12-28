"""
Custom Error and Warning names to improve readability
"""


# ----------
# Exceptions:
# ----------

class PlasmaPyError(Exception):
    """
    Base class of PlasmaPy custom errors.

    All custom exceptions raised by PlasmaPy should inherit from this class
    and be defined in this module.

    Custom exceptions can inherit from other exception types too. Thus, if code
    already knows how to handle a ValueError, it won't need any specific
    modification.
    """
    pass


class PhysicsError(PlasmaPyError, ValueError):
    """Error for use of a physics value outside PlasmaPy theoretical bounds"""
    pass


class RelativityError(PhysicsError):
    """Error for use of a speed greater than or equal to the speed of light"""
    pass


class AtomicError(PlasmaPyError):
    """Error for use by an atomic subpackage"""
    pass


class MissingAtomicDataError(AtomicError):
    """Error for use when atomic data is missing."""
    pass


class NoChargeInfoError(AtomicError):
    """Error for use when charge information is needed but missing."""


class IonError(NoChargeInfoError):
    """Error for use when an ion is invalid."""
    pass


class IsotopeError(AtomicError):
    """Error for use when an isotope is invalid."""
    pass


class ElementError(IsotopeError, IonError):
    """Error for use when an element is invalid."""
    pass


class ParticleError(ElementError):
    """Error for use when a particle is invalid."""
    pass


# ----------
# Warnings:
# ----------

class PlasmaPyWarning(Warning):
    """Base class of PlasmaPy custom warnings.

    All PlasmaPy custom warnings should inherit from this class and be defined
    in this module.

    Warnings should be issued using warnings.warn, which will not break
    execution if unhandled.
    """
    pass


class PhysicsWarning(PlasmaPyWarning):
    """Warning for using a mildly worrisome physics value"""
    pass


class RelativityWarning(PhysicsWarning):
    """Warning for use of a speed quantity approaching the speed of light"""
    pass


class AtomicWarning(PlasmaPyWarning):
    """Warnings for use in the atomic subpackage."""
    pass


class MissingAtomicDataWarning(AtomicWarning):
    """Warning for use when atomic data is missing."""
    pass
