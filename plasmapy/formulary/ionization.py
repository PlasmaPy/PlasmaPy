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


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n={"can_be_negative": False}
)
def Z_bal(n: u.m ** -3, T_e: u.K) -> u.dimensionless_unscaled:
    r"""
    Z_bal is the estimate average ionization level of a plasma in thermal equilibrium.
    Z_bal is derived from the Saha equation with the assumptions that
    the atoms are of a single species,
    are either hydrogenic or completely ionized,
    and that there is a balance between ionization and recombination,
    meaning that the number of atoms in either state are equal.
    The Saha equation and therefore Z_bal are more accurate when the plasma
    is at a high density and temperature.

    .. math::

        Z\_bal = \sqrt{\frac{k_B T_e}{E_H}} \sqrt{\ln{\frac{1}{4 n a_{0}^3}
        (\frac{k_B T_e}{\pi E_H})^{3/2}}} - \frac{1}{2}

    Where :math:`k_B` is the Boltzmann constant,
    :math:`a_0` is the Bohr radius, and
    :math:`E_H` is the ionization energy of Hydrogen

    Parameters
    ----------

    T_e : `~astropy.units.Quantity`
        The electron temperature.

    n : ~astropy.units.Quantity`
        The electron number density of the plasma.
    Returns
    -------
    Z: `~astropy.units.Quantity`
        The average ionization state of the ions in the plasma
        assuming that the number of ions in each state are equal.

    """
    E_H = 1 * u.Ry

    A = sqrt(k_B * T_e / E_H)
    B = log(1 / (4 * n * a0 ** 3) * (k_B * T_e / (pi * E_H)) ** (3 / 2))

    Z = A * sqrt(B) - 1 / 2
    return Z
