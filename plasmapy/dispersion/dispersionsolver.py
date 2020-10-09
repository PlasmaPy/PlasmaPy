__all__ = ["two_fluid_dispersion_solution"]

import numpy as np
import astropy.units as u

def two_fluid_dispersion_solution(n, B, T_i, theta, k, m_e=None, m_i=None, T_e=None, gamma_e=3., gamma_i=3):

    r"""
    Equation 34 of Bellan2012 (doi:10.1029/2012JA017856)

    .. math::
    \zeta^3 - A\zeta^2 + B\zeta - C = 0

    """

    if ((n is None) or (n <= 0)):
        raise ValueError("Number density must be a positive real number")

    if ((T_e is None) or (T_e <= 0)):
        raise ValueError("Electron temperature must be a positive real number")

    if ((T_i is None) or (T_i <= 0)):
        raise ValueError("Ion temperature must be a positive real number")

    if (theta is None):
        raise ValueError("Propagation direction can't be None")

    if ((k is None) or (k <= 0)):
        raise ValueError("Wave number must be a positive real number")

    # The required physcial constants
    gamma_sigma = 1
    m_i  = 1.672622E-24
    m_e  = 9.109384E-28
    mu_0 = 1
    e    = 4.8032E-10
    c    = 2.99792458E10

    # Required derived parameters
    c_s = np.sqrt(gamma_sigma * T * m_i)
    v_A = np.sqrt(B**2/(mu_0 * n * m_i))
    omega_ci = mu_0 * B/(c * m_i)
    omega_pe = np.sqrt((mu_0 * n * e**2)/m_e)

    alpha = np.cos(theta)**2
    beta = c_s**2/v_A**2
    Lambda = k**2 * v_A**2/omega_ci**2

    Q = 1 + k**2 * c**2/omega_pe**2

    A = (Q + Q**2 * beta + Q * alpha + alpha * Lambda)/Q**2
    B = alpha * (1 + 2 * Q * beta + Lambda * beta)/Q**2
    C = alpha**2 * beta/Q**2

    p = (3 * B - A**2)/3
    q = (9 * A * B - 2 * A**3 - 27 * C)/27


    keys = ['fast_mode', 'alfven_mode', 'acoustic_mode']
    zeta_sol = dict.fromkeys(keys)
    omega    = dict.fromkeys(keys)

    for (j,key) in zip(range(3), keys):

        # The solution corresponding to equation 37
        zeta_sol[key] = 2 * np.sqrt(-p/3) * np.cos(1/3 * np.arccos(3 * q/p
        * np.sqrt(-3/ p)) - 2 * np.pi/3 * j) + A/3

        # The solution corresponding to equation 38
        omega[key] = omega_ci * np.sqrt( 2 * Lambda * np.sqrt(-p/3)
        * np.cos(1/3 * np.arccos(3 * q/p * np.  sqrt(-3/ p)) - 2 * np.pi/3 * j) + A/3 )

    #omega = k**2 * v_A**2 * zeta_sol

    return zeta_sol, omega
