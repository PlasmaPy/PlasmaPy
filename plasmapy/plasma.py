"""
plasmapy.plasma
===============

Defines the core Plasma class used by PlasmaPy to represent plasma properties.
"""

import numpy as np
import astropy.units as u

class Plasma():
    @u.quantity_input(domain_x=u.m, domain_y=u.m, domain_z=u.m)
    def __init__(self, domain_x, domain_y, domain_z):
        self.x = domain_x
        self.y = domain_y
        self.z = domain_z

        x, y, z = self.x.si.value, self.y.si.value, self.z.si.value
        self.grid = np.meshgrid(x, y, z, indexing='ij') * u.m
