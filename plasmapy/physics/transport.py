"""Functions to calculate transport coefficients."""

from astropy import units
import numpy as np
from ..utils import _check_quantity
from ..constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, h, hbar,
                         ion_mass, charge_state)
from .parameters import Debye_length
from .quantum import deBroglie_wavelength


def Coulomb_logarithm(n_e, T, particles, V=None):
    r"""Estimates the Coulomb logarithm with an accuracy of order 10%.

    Parameters
    ----------
    n_e : Quantity
        The electron density in units convertible to per cubic meter.

    T : Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    particles : tuple containing two objects
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

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
    length: :math:`b_{min} = \lambda_D` which is a function of electron
    temperature and electron density.  At distances greater than the
    Debye length, charges from other particles will be screened out
    due to electrons rearranging themselves.

    The inner impact parameter is generally thought to be the maximum
    of :math:`b_\perp`, the angle at which a test particle is
    deflected by 90 degrees, and the test particle de Broglie wavelength,
    `\lambda_B`



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
    # and the relative velocity between particles.

    _check_quantity(T, 'T', 'Coulomb_logarithm', units.K)
    _check_quantity(n_e, 'n_e', 'Coulomb_logarithm', units.m**-3)

    if V is not None:
        _check_quantity(V, 'V', 'Coulomb_logarithm', units.m/units.s)

    if len(particles) != 2:
        raise ValueError("Incorrect number of particles in Coulomb_logarithm")

    m_test = ion_mass(particles[0])
    m_target = ion_mass(particles[1])

    reduced_mass = (m_test * m_target) / (m_test + m_target)

    q_test = np.abs(e*charge_state(particles[0]))
    q_target = np.abs(e*charge_state(particles[1]))

    # The outer impact parameter is the Debye length.  At distances
    # greater than the Debye length, the electrostatic potential of a
    # single particle is screened out by the electrostatic potentials
    # of other particles.  Past this distance, the electric fields of
    # individual particles do not affect each other much.  This
    # expression neglects screening by heavier ions.

    b_max = Debye_length(T, n_e)

    # The choice of inner impact parameter is more controversial.
    # There are two broad possibilities, and the expressions in the
    # literature often differ by factors of order unity or by
    # interchanging the reduced mass with the test particle mass.

    # The relative velocity is a source of uncertainty.  It is
    # reasonable to make an assumption relating the thermal energy to
    # the kinetic energy: reduced_mass*velocity**2 is approximately
    # equal to 3*k_B*T.  If a velocity was inputted, then we use that
    # instead.

    if V is None:
        velocity = np.sqrt(3*k_B*T/reduced_mass)
    else:
        velocity = V

    # The first possibility is that the inner impact parameter
    # corresponding to a deflection of 90 degrees, which is important
    # when classical effects dominate.

    b_perp = q_test*q_target/(4*pi*eps0*reduced_mass*velocity**2)

    # The second possibility is that the inner impact parameter is a
    # de Broglie wavelength.  There remains some ambiguity as to which
    # mass to choose to go into the de Broglie wavelength calculation.
    # Here we use the reduced mass, which will be comparable to the

    b_deBroglie = hbar/(2*reduced_mass*velocity)

    # Coulomb-style collisions will not happen for impact parameters
    # shorter than either of these two impact parameters, so we choose
    # the larger one.

    if b_perp > b_deBroglie:
        b_min = b_perp
    else:
        b_min = b_deBroglie

    # Now that we know how many approximations have to go into plasma
    # transport theory, we shall celebrate by returning the Coulomb
    # logarithm.

    lnLambda = np.log(b_max/b_min)

    return lnLambda
