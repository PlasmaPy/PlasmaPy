import numpy as np


def dispersion_solver(self, n, B, beta, k):

    r"""
    Equation 34 of Bellan2012 (doi:10.1029/2012JA017856)

    .. math::
    \zeta^3 - A\zeta^2 + B\zeta - C = 0

    """
    zeta = omega**2/(k**2 * v_A**2)
    beta = c_s**2/v_A**2
    Lambda = k**2 * v_A**2/omega_ci**2
    alpha = cos(theta)**2
    Q = 1 + k**2 * c**2/omega_pe**2

    A = (Q + Q**2 * beta + Q * alpha + alpha * Lambda)/Q**2
    B = alpha * (1 + 2 * Q * beta + Lambda * beta)/Q**2
    C = alpha**2 * beta/Q**2

    p = (3 * B - A**2)/3
    q = (9 * A * B - 2 * A**3 - 27 * C)/27

    zeta_sol = np.full((3), np.nan)

    for i in range(3):
        zeta_sol[i] = 2 * np.sqrt(-p/3) * np.cos(1/3 * np.arccos(3 * q/p * np.sqrt(-3/ p)) - 2 * np.pi/3 * i) + A/3

    return zeta_sol

