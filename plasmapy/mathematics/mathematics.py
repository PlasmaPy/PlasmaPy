"""Functions related to the plasma dispersion function"""

import numpy as np
from scipy import special
from astropy import units as u
from mpmath import polylog
from scipy.special import wofz as Faddeeva_function


def plasma_dispersion_func(zeta):
    r"""
    Calculate the plasma dispersion function

    Parameters
    ----------
    zeta : complex, int, float, ndarray, or Quantity
        Argument of plasma dispersion function.

    Returns
    -------
    Z : complex, float, or ndarray
        Value of plasma dispersion function.

    Raises
    ------
    TypeError
        If the argument is invalid.
    UnitsError
        If the argument is a Quantity but is not dimensionless
    ValueError
        If the argument is not entirely finite

    See also
    --------
    plasma_dispersion_func_deriv

    Notes
    -----
    The plasma dispersion function is defined as:

    .. math::
        Z(\zeta) = \pi^{-0.5} \int_{-\infty}^{+\infty}
        \frac{e^{-x^2}}{x-\zeta} dx

    where the argument is a complex number [fried.conte-1961]_.

    In plasma wave theory, the plasma dispersion function appears
    frequently when the background medium has a Maxwellian
    distribution function.  The argument of this function then refers
    to the ratio of a wave's phase velocity to a thermal velocity.

    References
    ----------
    .. [fried.conte-1961] Fried, Burton D. and Samuel D. Conte. 1961.
       The Plasma Dispersion Function: The Hilbert Transformation of the
       Gaussian. Academic Press (New York and London).

    Examples
    --------
    >>> plasma_dispersion_func(0)
    1.7724538509055159j
    >>> plasma_dispersion_func(1j)
    0.757872156141312j
    >>> plasma_dispersion_func(-1.52+0.47j)
    (0.6088888957234254+0.33494583882874024j)

    """

    if not isinstance(zeta, (int, float, complex, np.ndarray, u.Quantity)):
        raise TypeError("The argument to plasma_dispersion_function "
                        "must be one of the following types: complex, float, "
                        "int, ndarray, or Quantity.")

    if isinstance(zeta, u.Quantity):
        if zeta.unit == u.dimensionless_unscaled:
            zeta = zeta.value
        else:
            raise u.UnitsError("The argument to plasma_dispersion_function "
                               "must be dimensionless if it is a Quantity")

    if not np.all(np.isfinite(zeta)):
        raise ValueError("The argument to plasma_dispersion_function is "
                         "not finite.")

    Z = 1j * np.sqrt(np.pi) * Faddeeva_function(zeta)

    return Z


def plasma_dispersion_func_deriv(zeta):
    r"""Calculate the derivative of the plasma dispersion function

    Parameters
    ----------
    zeta : complex, int, float, ndarray, or Quantity
        Argument of plasma dispersion function.

    Returns
    -------
    Zprime : complex, int, float, or ndarray
        First derivative of plasma dispersion function.

    Raises
    ------
    TypeError
        If the argument is invalid.
    UnitsError
        If the argument is a Quantity but is not dimensionless
    ValueError
        If the argument is not entirely finite

    See also
    --------
    plasma_dispersion_func

    Notes
    -----
    The derivative of the plasma dispersion function is defined as:

    .. math::
        Z'(\zeta) = \pi^{-0.5} \int_{-\infty}^{+\infty}
        \frac{e^{-x^2}}{(x-\zeta)^2} dx

    where the argument is a complex number [fried.conte-1961]_.

    Examples
    --------
    >>> plasma_dispersion_func_deriv(0)
    (-2+0j)
    >>> plasma_dispersion_func_deriv(1j)
    (-0.48425568771737604+0j)
    >>> plasma_dispersion_func_deriv(-1.52+0.47j)
    (0.16587133149822897+0.44587978805935047j)

    """

    if not isinstance(zeta, (int, float, complex, np.ndarray, u.Quantity)):
        raise TypeError("The argument to plasma_dispersion_function_deriv "
                        "must be one of the following types: complex, float, "
                        "int, ndarray, or Quantity.")

    if isinstance(zeta, u.Quantity):
        if zeta.unit == u.dimensionless_unscaled:
            zeta = zeta.value
        else:
            raise u.UnitsError("The argument to "
                               "plasma_dispersion_function_deriv "
                               "must be dimensionless if it is a Quantity")

    if not np.all(np.isfinite(zeta)):
        raise ValueError("The argument to plasma_dispersion_function_deriv is "
                         "not finite.")

    Zprime = -2 * (1 + zeta * plasma_dispersion_func(zeta))

    return Zprime


def Fermi_integral(x, j):
    r"""Calculate the complete Fermi-Dirac integral.

    Parameters
    ----------
    x : Quantity
        Argument of the Fermi-Dirac integral function.

    j : Quantity
        Order/index of the Fermi-Dirac integral function.

    Returns
    -------
    integral : Quantity, complex
        Complete Fermi-Dirac integral for given argument and order.

    Raises
    ------
    TypeError
        If the argument is invalid.
    UnitsError
        If the argument is a Quantity but is not dimensionless
    ValueError
        If the argument is not entirely finite

    See also
    --------
    plasma_dispersion_func

    Notes
    -----
    The complete Fermi-Dirac integral is defined as: [1]_

    .. math::
        F_j (x) = \frac{1}{\Gamma (j+1)} \int_0^{\inf} \frac{t^j}{\exp{(t-x)} + 1} dt

    for j > 0.

    This is equivalent to the following polylogarithm [2]_ function:

    .. math::
        F_j (x) = -Li_{j+1}\left(-e^{x}\right)

    Warning: at present this function is limited to relatively small
    arguments due to limitations in the mpmath package's implementation
    of polylog.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Complete_Fermi-Dirac_integral
    .. [2] https://en.wikipedia.org/wiki/Polylogarithm

    Examples
    --------
    >>> Fermi_integral(0, 0)
    (0.6931471805599453-0j)
    >>> Fermi_integral(1, 0)
    (1.3132616875182228-0j)
    >>> Fermi_integral(1, 1)
    (1.8062860704447743-0j)
    """
    if isinstance(x, (int, float, complex)):
        arg = -np.exp(x)
        integral = -1 * complex(polylog(j + 1, arg))
        return integral
    elif isinstance(x, np.ndarray):
        integral_arr = np.zeros_like(x, dtype='complex')
        for idx, val in enumerate(x):
            integral_arr[idx] = -1 * complex(polylog(j + 1, -np.exp(val)))
        return integral_arr
    else:
        raise ValueError(f"Improper type {type(x)} given for argument x.")
