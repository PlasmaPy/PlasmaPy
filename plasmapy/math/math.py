"""Functions related to the plasma dispersion function"""

import numpy as np
from scipy import special
from astropy import units as u


def plasma_dispersion_func(zeta):
    r"""Calculate the plasma dispersion function

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
    Z(\zeta) = \sqrt{\pi}
    \int_{-\infty}^{+\infty} dx \frac{e^{-x^2}}{x-\zeta}

    where the argument is a complex number [fried.conte-1961].

    In plasma wave theory, the plasma dispersion function appears
    frequently when the background medium has a Maxwellian
    distribution function.  The argument of this function then refers
    to the ratio of a wave's phase velocity to a thermal velocity.

    References
    ----------
    .. [fried.conte-1961]
    Fried, Burton D. and Samuel D. Conte. 1961. The Plasma Dispersion
    Function: The Hilbert Transformation of the Gaussian. Academic
    Press (New York and London).

    Examples
    --------
    >>> plasma_dispersion_func(0)
    1.7724538509055159j
    >>> plasma_dispersion_func(1j)
    0.75787215614131187j
    >>> plasma_dispersion_func(-1.52+0.47j)
    (0.60888889572342553+0.33494583882874029j)

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

    Z = 1j * np.sqrt(np.pi) * np.exp(-zeta**2) * (1.0 + special.erf(1j * zeta))

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
    plasma_dispersion_func_deriv

    Notes
    -----
    The plasma dispersion function is defined as:

    .. math:: Z(\zeta) \equiv \sqrt{\pi}
    \int_{-\infty}^{+\infty} dx \frac{e^{-x^2}}{x-\zeta}

    where the argument is a complex number [fried.conte-1961].

    References
    ----------
    .. [fried.conte-1961]
    Fried, Burton D. and Samuel D. Conte. 1961. The Plasma Dispersion
    Function: The Hilbert Transformation of the Gaussian. Academic
    Press (New York and London).

    Examples
    --------
    >>> plasma_dispersion_func_deriv(0)
    (-2+0j)
    >>> plasma_dispersion_func_deriv(1j)
    (-0.48425568771737626+0j)
    >>> plasma_dispersion_func_deriv(-1.52+0.47j)
    (0.16587133149822941+0.44587978805935069j)

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
