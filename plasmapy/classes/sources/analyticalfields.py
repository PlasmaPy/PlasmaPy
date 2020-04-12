"""
Defines the AnalyticalFields class used by PlasmaPy to represent electric and magnetic fields
varying as known functions of position.
"""

__all__ = ["AnalyticalFields"]

import astropy.units as u

from typing import Callable
from plasmapy.classes import GenericPlasma
import numpy as np

_volt_over_meter = u.V / u.m  # for performance reasons


class AnalyticalFields(GenericPlasma):
    """
    Allows passing analytical functions as fields.

    This is primarily helpful for `plasmapy.simulation.ParticleTracker`.

    Parameters
    ----------
    magnetic_field : Callable
        Function of position returning the magnetic field.
    electric_field : Callable
        Function of position returning the electric field.

    Notes
    -----
    Units are assumed SI - they will be overwritten with Tesla and V/m, respectively.
    """

    def __init__(self, magnetic_field: Callable, electric_field: Callable):
        test_array = np.zeros((4, 3))
        if magnetic_field(test_array).shape != test_array.shape:
            self._interpolate_B = lambda r: np.array([magnetic_field(ri) for ri in r])
        else:
            self._interpolate_B = magnetic_field
        if electric_field(test_array).shape != test_array.shape:
            self._interpolate_E = lambda r: np.array([electric_field(ri) for ri in r])
        else:
            self._interpolate_E = electric_field

    def interpolate_B(self, r: u.m) -> u.T:
        """
        Returns the magnetic field, with units.

        Parameters
        ----------
        r : u.m
            Position; can be a (3,) array or a (N, 3) array.

        Returns
        -------
        Magnetic field array of the same shape as the input.
        """
        return u.Quantity(self._interpolate_B(r), u.T, copy=False)

    def interpolate_E(self, r: u.m) -> u.V / u.m:
        """
        Returns the electric field, with units.

        Parameters
        ----------
        r : u.m
            Position; can be a (3,) array or a (N, 3) array.

        Returns
        -------
        electric field array of the same shape as the input.
        """
        return u.Quantity(self._interpolate_E(r), _volt_over_meter, copy=False)
