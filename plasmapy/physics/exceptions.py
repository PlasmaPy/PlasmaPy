"""Exceptions and warnings related to plasma parameter calculations."""


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


