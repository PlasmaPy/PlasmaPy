"""
Module of dimensionless plasma parameters.

These are especially important for determining what regime a plasma is
in. (e.g., turbulent, quantum, collisional, etc.).

"""

import numpy as np
from astropy import units

from plasmapy.constants import k_B
from plasmapy.physics import quantum

def quantum_theta(T, n_e):
    """
    Compares Fermi energy to thermal kinetic energy to check if quantum
    effects are important.
    """
    fermi_energy = quantum.Fermi_energy(n_e)
    thermal_energy = k_B * T
    theta = thermal_energy / fermi_energy
    return theta
