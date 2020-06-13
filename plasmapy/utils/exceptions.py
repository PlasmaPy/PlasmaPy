"""Exceptions and warnings specific to PlasmaPy."""
__all__ = [
    "PlasmaPyError",
    "PhysicsError",
    "RelativityError",
    "PlasmaPyWarning",
    "CouplingWarning",
    "ImplicitUnitConversionWarning",
    "PhysicsWarning",
    "RelativityWarning",
]

from astropy.units import UnitsWarning

# ----------
# Exceptions
# ----------


class PlasmaPyError(Exception):
    """
    Base class of PlasmaPy custom errors.

    All custom exceptions raised by PlasmaPy should inherit from this
    class and be defined in this module.
    """

    pass


class PhysicsError(PlasmaPyError, ValueError):
    """
    The base exception for physics-related errors.
    """

    pass


class RelativityError(PhysicsError):
    """
    An exception for speeds greater than the speed of light.
    """

    pass


# ----------
# Warnings:
# ----------


class PlasmaPyWarning(Warning):
    """
    Base class of PlasmaPy custom warnings.

    All PlasmaPy custom warnings should inherit from this class and be
    defined in this module.

    Warnings should be issued using `~warnings.warn`, which will not break
    execution if unhandled.

    """

    pass


class PhysicsWarning(PlasmaPyWarning):
    """The base warning for warnings related to non-physical situations."""

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


class ImplicitUnitConversionWarning(PlasmaPyWarning, UnitsWarning):
    """
    A warning for an implicit conversion between equivalent :mod:`astropy` units.
    """

    pass
