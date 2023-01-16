"""
Module containing functionality focused on the plasma dispersion function
:math:`Z(ζ)`.
"""
__all__ = ["plasma_dispersion_func", "plasma_dispersion_func_deriv"]

__lite_funcs__ = ["plasma_dispersion_func_lite", "plasma_dispersion_func_deriv_lite"]

import astropy.units as u
import numbers
import numpy as np

from scipy.special import wofz as Faddeeva_function
from typing import Union

from plasmapy.utils.decorators import bind_lite_func, preserve_signature

__all__ += __lite_funcs__


# TODO: Use cython to speed up the Faddeeva_function execution in
#       plasma_dispersion_func_lite


@preserve_signature
def plasma_dispersion_func_lite(zeta):
    r"""
    The :term:`lite-function` version of
    `~plasmapy.dispersion.dispersionfunction.plasma_dispersion_func`.
    Performs the same calculation as
    `~plasmapy.dispersion.dispersionfunction.plasma_dispersion_func`,
    but is intended for computational use and, thus, has all data
    conditioning safeguards removed.

    Parameters
    ----------
    zeta : |array_like| of real or complex values
        Argument of the plasma dispersion function. ``zeta`` is
        dimensionless.

    Returns
    -------
    Z : |array_like| of complex values
        Value of plasma dispersion function.

    See Also
    --------
    ~plasmapy.dispersion.dispersionfunction.plasma_dispersion_func

    """

    return 1j * np.sqrt(np.pi) * Faddeeva_function(zeta)


@bind_lite_func(plasma_dispersion_func_lite)
def plasma_dispersion_func(
    zeta: Union[complex, int, float, np.ndarray, u.Quantity]
) -> Union[complex, float, np.ndarray, u.Quantity]:
    r"""
    Calculate the plasma dispersion function.

    Parameters
    ----------
    zeta : complex, int, float, ~numpy.ndarray, or ~astropy.units.Quantity
        Argument of plasma dispersion function.

    Returns
    -------
    Z : complex, float, or ~numpy.ndarray
        Value of plasma dispersion function.

    Raises
    ------
    TypeError
        If the argument is of an invalid type.

    ~astropy.units.UnitsError
        If the argument is a `~astropy.units.Quantity` but is not
        dimensionless.

    ValueError
        If the argument is not entirely finite.

    See Also
    --------
    plasma_dispersion_func_deriv

    Notes
    -----
    The plasma dispersion function is defined as:

    .. math::
        Z(\zeta) = \pi^{-0.5} \int_{-\infty}^{+\infty}
        \frac{e^{-x^2}}{x-\zeta} dx

    where the argument is a complex number :cite:p:`fried:1961`.

    In plasma wave theory, the plasma dispersion function appears
    frequently when the background medium has a Maxwellian
    distribution function.  The argument of this function then refers
    to the ratio of a wave's phase velocity to a thermal velocity.

    Examples
    --------
    >>> plasma_dispersion_func(0)
    1.7724538509055159j
    >>> plasma_dispersion_func(1j)
    0.757872156141312j
    >>> plasma_dispersion_func(-1.52+0.47j)
    (0.6088888957234254+0.33494583882874024j)

    For user convenience
    `~plasmapy.dispersion.dispersionfunction.plasma_dispersion_func_lite`
    is bound to this function and can be used as follows:

    >>> plasma_dispersion_func.lite(0)
    1.7724538509055159j
    >>> plasma_dispersion_func.lite(1j)
    0.757872156141312j
    >>> plasma_dispersion_func.lite(-1.52+0.47j)
    (0.6088888957234254+0.33494583882874024j)

    """
    if not isinstance(
        zeta, (numbers.Integral, numbers.Real, numbers.Complex, np.ndarray, u.Quantity)
    ):
        raise TypeError(
            "The argument to plasma_dispersion_function "
            "must be one of the following types: complex, float, "
            "int, ndarray, or Quantity."
        )

    if isinstance(zeta, u.Quantity):
        if zeta.unit == u.dimensionless_unscaled:
            zeta = zeta.value
        else:
            raise u.UnitsError(
                "The argument to plasma_dispersion_function "
                "must be dimensionless if it is a Quantity"
            )

    if not np.all(np.isfinite(zeta)):
        raise ValueError("The argument to plasma_dispersion_function is not finite.")

    return plasma_dispersion_func_lite(zeta)


@preserve_signature
def plasma_dispersion_func_deriv_lite(zeta):
    r"""
    The :term:`lite-function` version of
    `~plasmapy.dispersion.dispersionfunction.plasma_dispersion_func_deriv`.
    Performs the same calculation as
    `~plasmapy.dispersion.dispersionfunction.plasma_dispersion_func_deriv`,
    but is intended for computational use and, thus, has all data
    conditioning safeguards removed.

    Parameters
    ----------
    zeta : |array_like| of real or complex values
        Argument of the plasma dispersion function. ``zeta`` is
        dimensionless.

    Returns
    -------
    Zprime : |array_like| of complex values
        First derivative of plasma dispersion function.

    See Also
    --------
    ~plasmapy.dispersion.dispersionfunction.plasma_dispersion_func_deriv

    """

    return -2 * (1 + zeta * plasma_dispersion_func_lite(zeta))


@bind_lite_func(plasma_dispersion_func_deriv_lite)
def plasma_dispersion_func_deriv(
    zeta: Union[complex, int, float, np.ndarray, u.Quantity]
) -> Union[complex, float, np.ndarray, u.Quantity]:
    r"""
    Calculate the derivative of the plasma dispersion function.

    Parameters
    ----------
    zeta : complex, int, float, ~numpy.ndarray, or ~astropy.units.Quantity
        Argument of plasma dispersion function.

    Returns
    -------
    Zprime : complex, float, or ~numpy.ndarray
        First derivative of plasma dispersion function.

    Raises
    ------
    TypeError
        If the argument is invalid.

    ~astropy.units.UnitsError
        If the argument is a `~astropy.units.Quantity` but is not
        dimensionless.

    ValueError
        If the argument is not entirely finite.

    See Also
    --------
    plasma_dispersion_func

    Notes
    -----
    The derivative of the plasma dispersion function is defined as:

    .. math::
        Z'(\zeta) = \pi^{-1/2} \int_{-\infty}^{+\infty}
        \frac{e^{-x^2}}{(x-\zeta)^2} dx

    where the argument is a complex number :cite:p:`fried:1961`.

    Examples
    --------
    >>> plasma_dispersion_func_deriv(0)
    (-2+0j)
    >>> plasma_dispersion_func_deriv(1j)
    (-0.484255687717376...+0j)
    >>> plasma_dispersion_func_deriv(-1.52+0.47j)
    (0.165871331498228...+0.445879788059350...j)

    For user convenience
    `~plasmapy.dispersion.dispersionfunction.plasma_dispersion_func_deriv_lite`
    is bound to this function and can be used as follows:

    >>> plasma_dispersion_func_deriv.lite(0)
    (-2+0j)
    >>> plasma_dispersion_func_deriv.lite(1j)
    (-0.484255687717376...+0j)
    >>> plasma_dispersion_func_deriv.lite(-1.52+0.47j)
    (0.165871331498228...+0.445879788059350...j)

    """

    if not isinstance(
        zeta, (numbers.Integral, numbers.Real, numbers.Complex, np.ndarray, u.Quantity)
    ):
        raise TypeError(
            "The argument to plasma_dispersion_function_deriv "
            "must be one of the following types: complex, float, "
            "int, ndarray, or Quantity."
        )

    if isinstance(zeta, u.Quantity):
        if zeta.unit == u.dimensionless_unscaled:
            zeta = zeta.value
        else:
            raise u.UnitsError(
                "The argument to plasma_dispersion_function_deriv "
                "must be dimensionless if it is a Quantity."
            )

    if not np.all(np.isfinite(zeta)):
        raise ValueError(
            "The argument to plasma_dispersion_function_deriv is not finite."
        )

    return plasma_dispersion_func_deriv_lite(zeta)
