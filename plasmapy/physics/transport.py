"""Functions to calculate transport coefficients."""

from astropy import units
import numpy as np

import plasmapy.atomic as atomic
import plasmapy.utils as utils

from ..constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, h, hbar)
from .parameters import Debye_length


@utils.check_quantity({"T": {"units": units.K, "can_be_negative": False},
                       "n_e": {"units": units.m**-3}
                       })
def Coulomb_logarithm(T, n_e, particles, V=None):
    r"""Estimates the Coulomb logarithm.

    Parameters
    ----------

    T : Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    n_e : Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    V : Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 3 k_B T`
        where `mu` is the reduced mass.

    Returns
    -------
    lnLambda : float or numpy.ndarray
        An estimate of the Coulomb logarithm that is accurate to
        roughly its reciprocal.

    Raises
    ------
    ValueError
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    UnitConversionError
        If the units on any of the inputs are incorrect

    UserWarning
        If the inputted velocity is greater than 80% of the speed of
        light.

    TypeError
        If the n_e, T, or V are not Quantities.

    Notes
    -----
    The Coulomb logarithm is given by

    .. math::
        \ln{\Lambda} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    where :math:`b_{min}` and :math:`b_{max}` are the inner and outer
    impact parameters for Coulomb collisions [1]_.

    The outer impact parameter is given by the Debye length:
    :math:`b_{min} = \lambda_D` which is a function of electron
    temperature and electron density.  At distances greater than the
    Debye length, electric fields from other particles will be
    screened out due to electrons rearranging themselves.

    The choice of inner impact parameter is more controversial. There
    are two main possibilities.  The first possibility is that the
    inner impact parameter corresponds to a deflection angle of 90
    degrees.  The second possibility is that the inner impact
    parameter is a de Broglie wavelength, :math:`\lambda_B`
    corresponding to the reduced mass of the two particles and the
    relative velocity between collisions.  This function uses the
    standard practice of choosing the inner impact parameter to be the
    maximum of these two possibilities.  Some inconsistencies exist in
    the literature on how to define the inner impact parameter [2]_.

    Errors associated with the Coulomb logarithm are of order its
    inverse If the Coulomb logarithm is of order unity, then the
    assumptions made in the standard analysis of Coulomb collisions
    are invalid.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*units.m**-3
    >>> T = 1e6*units.K
    >>> particles = ('e', 'p')
    >>> Coulomb_logarithm(T, n, particles)
    14.748259780491056
    >>> Coulomb_logarithm(T, n, particles, 1e6*u.m/u.s)
    11.363478214139432

    References
    ----------
    .. [1] Physics of Fully Ionized Gases, L. Spitzer (1962)

    .. [2] Comparison of Coulomb Collision Rates in the Plasma Physics
       and Magnetically Confined Fusion Literature, W. Fundamenski and
       O.E. Garcia, EFDA–JET–R(07)01
       (http://www.euro-fusionscipub.org/wp-content/uploads/2014/11/EFDR07001.pdf)

    """

    if not isinstance(particles, (list, tuple)) or len(particles) != 2:
        raise ValueError("The third input of Coulomb_logarithm must be a list "
                         "or tuple containing representations of two charged "
                         "particles.")

    masses = np.zeros(2)*units.kg
    charges = np.zeros(2)*units.C

    for particle, i in zip(particles, range(2)):

        try:
            masses[i] = atomic.ion_mass(particles[i])
        except Exception:
            raise ValueError("Unable to find mass of particle: " +
                             str(particles[i]) + " in Coulomb_logarithm.")

        try:
            charges[i] = np.abs(e*atomic.charge_state(particles[i]))
            if charges[i] is None:
                raise
        except Exception:
            raise ValueError("Unable to find charge of particle: " +
                             str(particles[i]) + " in Coulomb_logarithm.")

    reduced_mass = masses[0] * masses[1] / (masses[0] + masses[1])

    # The outer impact parameter is the Debye length.  At distances
    # greater than the Debye length, the electrostatic potential of a
    # single particle is screened out by the electrostatic potentials
    # of other particles.  Past this distance, the electric fields of
    # individual particles do not affect each other much.  This
    # expression neglects screening by heavier ions.

    T = T.to(units.K, equivalencies=units.temperature_energy())

    b_max = Debye_length(T, n_e)

    # The choice of inner impact parameter is more controversial.
    # There are two broad possibilities, and the expressions in the
    # literature often differ by factors of order unity or by
    # interchanging the reduced mass with the test particle mass.

    # The relative velocity is a source of uncertainty.  It is
    # reasonable to make an assumption relating the thermal energy to
    # the kinetic energy: reduced_mass*velocity**2 is approximately
    # equal to 3*k_B*T.

    # If no relative velocity is inputted, then we make an assumption
    # that relates the thermal energy to the kinetic energy:
    # reduced_mass*velocity**2 is approximately equal to 3*k_B*T.

    if V is None:
        V = np.sqrt(3 * k_B * T / reduced_mass)
    else:
        utils._check_relativistic(V, 'Coulomb_logarithm', betafrac=0.8)

    # The first possibility is that the inner impact parameter
    # corresponds to a deflection of 90 degrees, which is valid when
    # classical effects dominate.

    b_perp = charges[0] * charges[1] / (4 * pi * eps0 * reduced_mass * V**2)

    # The second possibility is that the inner impact parameter is a
    # de Broglie wavelength.  There remains some ambiguity as to which
    # mass to choose to go into the de Broglie wavelength calculation.
    # Here we use the reduced mass, which will be of the same order as
    # mass of the smaller particle and thus the longer de Broglie
    # wavelength.

    b_deBroglie = hbar / (2 * reduced_mass * V)

    # Coulomb-style collisions will not happen for impact parameters
    # shorter than either of these two impact parameters, so we choose
    # the larger of these two possibilities.

    b_min = np.zeros_like(b_perp)

    for i in range(b_min.size):

        if b_perp.flat[i] > b_deBroglie.flat[i]:
            b_min.flat[i] = b_perp.flat[i]
        else:
            b_min.flat[i] = b_deBroglie.flat[i]

    # Now that we know how many approximations have to go into plasma
    # transport theory, we shall celebrate by returning the Coulomb
    # logarithm.

    ln_Lambda = np.log(b_max/b_min)
    ln_Lambda = ln_Lambda.to(units.dimensionless_unscaled).value

    return ln_Lambda
