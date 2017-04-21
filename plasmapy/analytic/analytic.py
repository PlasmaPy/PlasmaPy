"""Functions related to the plasma dispersion function"""

def plasma_dispersion_func(zeta):
    """Calculate the plasma dispersion function

    Parameters
    ----------
    zeta : complex number or array
        Argument of plasma dispersion function.

    Returns
    -------
    Z : complex number or array
        Value of plasma dispersion function.

    Raises
    ------
    ValueError
        Invalid argument.

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
    Here we calculate Z(0) and Z(1j).

    >>> from plasmapy.analytic_functions import plasma_dispersion_func
    >>> Z_of_0 = plasma_dispersion_func(0)
    >>> print(Z_of_0)
    1.77245385091j
    >>> Z_of_1j = plasma_dispersion_func(1j)
    >>> print(Z_of_1j)
    >>> Z_of_complex_number = plasma_dispersion_func(-1.52+0.47j)
    print(Z_of_complex_number)
    (0.608888895723+0.334945838829j)

    """

    import numpy as np
    from scipy import special

    if not np.all(np.isfinite(zeta)):
        raise ValueError("Argument of plasma_dispersion_function " + \
                             "must be a complex number or array")

    Z = 1j * np.sqrt(np.pi) * np.exp(-zeta**2) * \
        (1.0 + special.erf(1j*zeta))

    return Z

def plasma_dispersion_func_deriv(zeta):
    """Calculate the derivative of the plasma dispersion function
    
    Parameters
    ----------
    zeta : complex number or array
        Argument of plasma dispersion function.
 
    Returns
    -------
    Zprime : complex number or array
        Derivative of plasma dispersion function.

    Raises
    ------
    ValueError
        Invalid argument.

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

    """
    import numpy as np

    if not np.all(np.isfinite(zeta)):
        raise ValueError("Argument of plasma_dispersion_function " + \
                             "must be a complex number or array")
    
    Zprime = -2*(1+zeta*plasma_dispersion_func(zeta))

    return Zprime
