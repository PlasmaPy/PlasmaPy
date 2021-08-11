"""Exceptions and warnings for functionality defined in `plasmapy.plasma`."""
__all__ = ["DataStandardError"]

from plasmapy.utils import PlasmaPyError


class DataStandardError(PlasmaPyError):
    """An exception for when HDF5 is not defined by OpenPMD standard."""

    pass
