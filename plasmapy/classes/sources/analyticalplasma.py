"""
Defines the core Plasma class used by PlasmaPy to represent plasma properties.
"""
import astropy.units as u

from plasmapy.classes import GenericPlasma

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
        self.interpolate_B = lambda r: u.Quantity(len(r) * [magnetic_field(r)])  # TODO I think this may be off...
        self.interpolate_E = lambda r: u.Quantity(len(r) * [electric_field(r)])

    @classmethod
    def is_datasource_for(cls, **kwargs):
        match = 'interpolate_E' in kwargs.keys() and 'interpolate_B' in kwargs.keys()
        return match
