"""Functions to calculate transport coefficients."""

from astropy import units
import numpy as np
from ..utils import _check_quantity
from ..constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, h, hbar,
                         ion_mass, charge_state)
from .parameters import Debye_length
from .quantum import deBroglie_wavelength


def Coulomb_logarithm(T=1e6*units.K, n_e=1e19*units.m**-3,
                      particles=('e', 'p')):
    r"""Estimates the Coulomb logarithm.

    Parameters
    ----------
    T : Quantity
        Temperature in units of temperature or energy per particle.

    n_e : Quantity
        The electron density in units convertible to per cubic meter,
        defaulting to

    particles : tuple containing two objects
        A tuple containing

    Returns
    -------
    lnLambda : float or numpy.ndarray
        An estimate of the Coulomb logarithm

    Raises
    ------



    Notes
    -----
    The Coulomb logarithm is given by

    .. math::
    \ln{\Lambda} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    where :math:`b_{min}` and :math:`b_{max}` are the inner and outer
    impact parameters.

    The outer impact parameter is generally accepted to be the Debye
    length: `b_{min} = \lambda_D`.  At distances greater than the
    Debye length,

    At this distance, charges from particles are screened out by other
    particles.

    There is some disagreement on what the inner impact parameter
    should be.

    Examples
    --------
    >>> from astropy import units as u
    >>> Coulomb_logarithm(T=1e6*units.K, n_e=1e19*units.m**-3)

    See also
    --------

    References
    ----------

    [1] Bittencourt

    [2] Mulser, Alber, and Murakami (2014)

    """

    # The temperatures of the two species are assumed to be the same.
    # The temperature comes up in the calculation of the Debye length
    # and the relative velocity between particles.  If the collisions
    # involve electrons, it would probably make more sense to choose
    # the electron temperature.

    _check_quantity(T, 'T', 'Coulomb_logarithm', units.K)
    _check_quantity(n_e, 'n_e', 'Coulomb_logarithm', units.m**-3)

    if len(particles) != 2:
        raise ValueError("Incorrect number of particles in Coulomb_logarithm")

    # The outer impact parameter is the Debye length.  At distances
    # greater than the Debye length, the electrostatic potential of a
    # single particle is screened out by the electrostatic potentials
    # of other particles.  Past this distance, the electric fields of
    # individual particles do not affect each other.

    bmax = Debye_length(T, n_e)

    # The choice of inner impact parameter is somewhat more
    # controversial.


    # Heuristic arguments suggest that the minimum impact paramter is
    # the maximum of either the de Broglie wavelength (below which
    # quantum effects dominate) or the angle of perpendicular
    # deflection.

    q1 = np.abs(e*charge_state(particles[0]))
    q2 = np.abs(e*charge_state(particles[1]))

    bperp = q1*q2/(12*pi*eps0*k_B*T) # Bittencourt eq



    # The second possibility is the electron de Broglie wavelength.
    # The wave nature of ions can usually be neglected (Spitzer 1962).
    # 

    m1 = ion_mass(particles[0])
    m2 = ion_mass(particles[1])

    reduced_mass = m1*m2/(m1+m2)

    # What thermal velocity?  How to define this between arbitrary particles?

#    b_deBroglie = hbar/

    bmin = bperp  # temporary

    # The minimum

    lnLambda = np.log(bmax/bmin)

    return lnLambda


