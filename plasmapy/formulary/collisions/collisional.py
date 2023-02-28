"""
Module containing the Collisional Analysis formulation.
"""
__all__ = [
    "coal"
]

import astropy.units as u
import math
import numpy as np

from astropy.constants.si import e
from numbers import Real

from plasmapy import particles, utils
from plasmapy.formulary.collisions import frequencies


def coal(
    r_0: u.m,
    r_n: u.m,
    n_p: 1/u.cm**-3,
    eta,
    v_p: u.m/u.s,
    T_p: u.K,
    theta,
    n_step: int = 1000
):
    r"""
    Compute the

    Parameters
    ----------
    r_0 : `~astropy.units.Quantity`
        Starting position in units of meters or astronomical units.

    r_1 : `~astropy.units.Quantity`
        Final position in units of meters or astronomical units.

    n_p : `~astropy.units.Quantity`
        The proton number density in units convertible to m\ :sup:`-3`.

    eta :
        The alpha number density divided by the proton number density,
        densities should be in units convertible to m\ :sup:`-3`.

    v_p : `~astropy.units.Quantity`, optional
        The velocity of the protons in units convertible to ms\ :sup:`-1`.

    T_p : `~astropy.units.Quantity`
        The proton temperature in units of :math:`K` or :math:`eV`.

    n_step : `int`, optional
        The step number for the operation

    Returns
    -------
    theta : `float`
        The dimensionless alpha-proton temperature ratio prediction
        for distance provided

    Raises
    ------

    Notes
    -----
    Using th

    Examples
    --------


    """

    # Initialize the alpha-proton charge and mass ratios.
    z_a = 2.
    mu_a = 4.

    # Initialise.
    d_r = (r_0 - r_n) / (1. * n_step)

    # Loop.
    for i in range(n_step):

        r = r_n + ((i + 1) * d_r)

        n_p = n_p * (r / r_n) ** -1.8
        v_p = v_p * (r / r_n) ** -0.2
        t_p = t_p * (r / r_n) ** -0.77

        alpha = (theta_ap + mu_a)
        beta = theta_ap

        if alpha == 0:
            alpha = float('Nan')
        if beta == 0:
            beta = float('Nan')

        charlie = 1 + (z_a ** 2 * eta / beta)
        if charlie < 0:
            charlie = 0

        arg_ = ((n_p ** 0.5 / t_p ** 1.5) * (z_a * (mu_a + 1) / alpha) *
                (charlie) ** 0.5)

        if arg_ == 0:
            arg_ = math.exp(9)
        elif arg_ < 0:
            arg_ = math.exp(9)
        else:
            pass

        lambda_ap = 9 - np.log(arg_)

        x = (v_p * t_p ** 1.5)
        y = (eta + 1) ** 2.5

        if theta_ap == float('Nan'):
            z = float('Nan')
        else:
            z = (theta_ap + mu_a) ** 1.5

        if x == 0:
            x = float('Nan')
        if y == 0:
            y = float('Nan')
        if z == 0:
            z = float('Nan')
        elif z < 0:
            z = float('Nan')

        d_theta_ap = ((-2.60e7) * ((n_p / x)) * (mu_a ** 0.5 * z_a ** 2 / y) * (
                    (theta_ap - 1.) * (eta * theta_ap + 1.) ** 2.5 / z) * (lambda_ap) * (d_r))

        theta_ap = theta_ap + d_theta_ap

    return theta_ap


