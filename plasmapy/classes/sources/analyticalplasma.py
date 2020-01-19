"""
Defines the core Plasma class used by PlasmaPy to represent plasma properties.
"""
import astropy.units as u

from plasmapy.classes import GenericPlasma
import numpy as np

__all__ = [
    "AnalyticalPlasma"
]


class AnalyticalPlasma(GenericPlasma):
    """
    Work-in-progress class for passing analytical functions as fields 

    This is primarily helpful for `plasmapy.simulation.ParticleTracker`.
    """

    def __init__(self, magnetic_field, electric_field):
        """
        Initialize plasma paramters.
        The most basic description is composition (ion), temperature,
        density, and ionization.
        """
        test_array = np.zeros((4, 3))
        if magnetic_field(test_array).shape != test_array.shape:
            self._interpolate_B = lambda r: np.array([magnetic_field(ri) for ri in r])
        else: 
            self._interpolate_B = magnetic_field
        if electric_field(test_array).shape != test_array.shape:
            self._interpolate_E = lambda r: np.array([electric_field(ri) for ri in r])
        else: 
            self._interpolate_E = electric_field

    def interpolate_B(self, r):
        return u.Quantity(self._interpolate_B(r), u.T, copy = False)

    def interpolate_E(self, r):
        return u.Quantity(self._interpolate_E(r), u.V / u.m, copy = False)

    @classmethod
    def is_datasource_for(cls, **kwargs):
        match = 'interpolate_E' in kwargs.keys() and 'interpolate_B' in kwargs.keys()
        return match
