""" This module gathers functions relating to ionization states and the properties thereof.
"""
__all__ = ["Z_bal"]

from astropy import units as u
from plasmapy.utils.decorators import (
    angular_freq_to_hz,
    check_relativistic,
    validate_quantities,
)
from astropy.constants import c, a0, k_B
from numpy import pi, exp, sqrt, log


def Z_bal(n: u.m ** -3, T: u.K):
    r"""
    Z_bal is derived from the Saha equation with the assumptions that
    the atoms are of a single species,
    are either hydrogenic or completely ionized,
    and that there is a balance between ionization and recombination,
    meaning that the number of atoms in either state are equal.
    The Saha equation and therefore Z_bal are more accurate when the plasma
    is at a high density and temperature.

    Parameters
    ----------

    T : `~astropy.units.Quantity`
        The electron temperature.

    n : ~astropy.units.Quantity`
        The electron number density of the plasma.
    Returns
    -------
    Z: `~astropy.units.Quantity`
        The average ionization state of the ions in the plasma
        assuming that the number of ions in each state are equal.

    """
    E_H = 13.6 * u.eV

    A = sqrt(k_B * T / E_H)
    B = log(1 / (4 * n * a0 ** 3) * (k_B * T / (pi * E_H)) ** (3 / 2))

    Z = A * sqrt(B) - 1 / 2
    return Z
