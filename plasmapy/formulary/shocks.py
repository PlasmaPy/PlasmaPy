""" This module gathers functions relating to shocks and the properties thereof.
"""
__all__ = [""]


from astropy import units as u
from plasmapy.utils.decorators import (
    angular_freq_to_hz,
    check_relativistic,
    validate_quantities,
)
from astropy.constants import c, a0, k_B
from numpy import pi, exp, sqrt, log


validate_quantities(
    rho_1={"can_be_negative": False}, rho_2={"can_be_negative": False})
def entropy_across_shock_polytorpic(c_v: u.J / u.K, p_1: u.Bar, p_2: u.Bar, rho_1: u.kg / u.m ** 3, rho_2: u.kg / u.m ** 3, gamma) -> u.J / u.K:
    r"""


     .. math::

        s_2 - s_1 =
        c_v ln\left[ \frac{p_2}{p_1} \left( \frac{\rho_1}{\rho_2} \right)^{\gamma} \right]
    Here variables indexed by 1 and 2 are upstream (pre shock)
    and downstream (post shock) respectively.
    """
