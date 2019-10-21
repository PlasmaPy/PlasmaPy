"""
Defines the core Plasma class used by PlasmaPy to represent plasma properties.
"""
import warnings

import astropy.units as u

from plasmapy.formulary.parameters import _grab_charge
from plasmapy.formulary.dimensionless import quantum_theta
from plasmapy.formulary import coupling_parameter
from plasmapy.atomic import particle_mass

from plasmapy.utils import CouplingWarning
from plasmapy.utils.pytest_helpers import call_string

from plasmapy.classes import GenericPlasma

__all__ = [
    "PlasmaBlob"
]


class AnalyticalPlasma(GenericPlasma):
    """
    Class for describing and calculating plasma parameters without
    spatial/temporal description.
    """

    def __init__(self, magnetic_field, electric_field):
        """
        Initialize plasma paramters.
        The most basic description is composition (ion), temperature,
        density, and ionization.
        """
        self.interpolate_B = magnetic_field
        self.interpolate_E = electric_field

    # @classmethod
    # def is_datasource_for(cls, **kwargs):
    #     match = 'E' in kwargs.keys() and 'n_e' in kwargs.keys()
    #     return match
