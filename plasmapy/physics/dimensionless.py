"""
Module of dimensionless plasma parameters.

These are especially important for determining what regime a plasma is
in. (e.g., turbulent, quantum, collisional, etc.).

"""

import numpy as np
from astropy import units

from plasmapy.constants import k_B
from plasmapy.physics import quantum, parameters

def quantum_theta(T, n_e):
    """
    Compares Fermi energy to thermal kinetic energy to check if quantum
    effects are important.
    """
    fermi_energy = quantum.Fermi_energy(n_e)
    thermal_energy = k_B * T
    theta = thermal_energy / fermi_energy
    return theta


def beta(T, n, B):
    """
    The ratio of thermal pressure to magnetic pressure.
    """
    thermal_pressure = parameters.thermal_pressure(T, n)
    magnetic_pressure = parameters.magnetic_pressure(B)
    return thermal_pressure / magnetic_pressure
