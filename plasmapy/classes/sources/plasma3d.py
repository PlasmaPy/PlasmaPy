"""
Defines the core Plasma class used by PlasmaPy to represent plasma properties.
"""
import warnings

import numpy as np
import astropy.units as u

from plasmapy.constants import (m_p,
                                m_e,
                                c,
                                mu0,
                                k_B,
                                e,
                                eps0,
                                pi,
                                )

from plasmapy.classes import GenericPlasma

__all__ = [
    "Plasma3D"
]


class Plasma3D(GenericPlasma):
    """
    Core class for describing and calculating plasma parameters with
    spatial dimensions.

    Attributes
    ----------
    x : `astropy.units.Quantity`
        x-coordinates within the plasma domain. Equal to the
        `domain_x` input parameter.
    y : `astropy.units.Quantity`
        y-coordinates within the plasma domain. Equal to the
        `domain_y` input parameter.
    z : `astropy.units.Quantity`
        z-coordinates within the plasma domain. Equal to the
        `domain_z` input parameter.
    grid : `astropy.units.Quantity`
        (3, x, y, z) array containing the values of each coordinate at
        every point in the domain.
    domain_shape : tuple
        Shape of the plasma domain.
    density : `astropy.units.Quantity`
        (x, y, z) array of mass density at every point in the domain.
    momentum : `astropy.units.Quantity`
        (3, x, y, z) array of the momentum vector at every point in
        the domain.
    pressure : `astropy.units.Quantity`
        (x, y, z) array of pressure at every point in the domain.
    magnetic_field : `astropy.units.Quantity`
        (3, x, y, z) array of the magnetic field vector at every point
        in the domain.

    Parameters
    ----------
    domain_x : `astropy.units.Quantity`
        1D array of x-coordinates for the plasma domain. Must have
        units convertable to length.
    domain_y : `astropy.units.Quantity`
        1D array of y-coordinates for the plasma domain. Must have
        units convertable to length.
    domain_z : `astropy.units.Quantity`
        1D array of z-coordinates for the plasma domain. Must have
        units convertable to length.

    """
    @u.quantity_input(domain_x=u.m, domain_y=u.m, domain_z=u.m)
    def __init__(self, domain_x, domain_y, domain_z):
        #GenericPlasma.__init__(self, data, header, **kwargs)

        # Define domain sizes
        self.x = domain_x
        self.y = domain_y
        self.z = domain_z

        self.grid = np.array(np.meshgrid(self.x, self.y, self.z,
                                         indexing='ij'))
        self.domain_shape = (len(self.x), len(self.y), len(self.z))

        # Initiate core plasma variables
        self.density = np.zeros(self.domain_shape) * u.kg / u.m**3
        self.momentum = np.zeros((3, *self.domain_shape)) * u.kg / (u.m ** 2 * u.s)
        self.pressure = np.zeros(self.domain_shape) * u.Pa
        self.magnetic_field = np.zeros((3, *self.domain_shape)) * u.T
        self.electric_field = np.zeros((3, *self.domain_shape)) * u.V / u.m

    @property
    def velocity(self):
        return self.momentum / self.density

    @property
    def magnetic_field_strength(self):
        B = self.magnetic_field
        return np.sqrt(np.sum(B * B, axis=0))

    @property
    def electric_field_strength(self):
        E = self.electric_field
        return np.sqrt(np.sum(E * E, axis=0))

    @property
    def alfven_speed(self):
        B = self.magnetic_field
        rho = self.density
        return np.sqrt(np.sum(B * B, axis=0) / (mu0 * rho))

    @classmethod
    def is_datasource_for(cls, **kwargs):
        if len(kwargs) == 3:
            match = kwargs.get('domain_x') and \
                    kwargs.get('domain_y') and \
                    kwargs.get('domain_z')
        else:
            match = False
        return match
