"""Functions to calculated fundamental plasma length parameters."""
__all__ = ["Debye_length"]
__aliases__ = ["lambdaD_"]

import astropy.units as u
import numpy as np

from astropy.constants.si import e, eps0, k_B

from plasmapy.utils.decorators import validate_quantities

__all__ += __aliases__


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def Debye_length(T_e: u.K, n_e: u.m ** -3) -> u.m:
    r"""Calculate the characteristic decay length for electric fields,
     due to charge screening.

    **Aliases:** `lambdaD_`

    Parameters
    ----------
    T_e : `~astropy.units.Quantity`
        Electron temperature.

    n_e : `~astropy.units.Quantity`
        Electron number density.

    Returns
    -------
    lambda_D : `~astropy.units.Quantity`
        The Debye length in meters.

    Raises
    ------
    `TypeError`
        If either argument is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If either argument is in incorrect units.

    `ValueError`
        If either argument contains invalid values.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    The Debye length is the exponential scale length for charge
    screening and is given by

    .. math::
        λ_D = \sqrt{\frac{ε_0 k_b T_e}{n_e e^2}}

    for an electron plasma with nearly stationary ions.

    The electrical potential will drop by a factor of :math:`1/e` every Debye
    length.

    Plasmas will generally be quasineutral on length scales significantly
    larger than the Debye length.

    See Also
    --------
    Debye_number

    Examples
    --------
    >>> from astropy import units as u
    >>> Debye_length(5e6*u.K, 5e15*u.m**-3)
    <Quantity 0.002182... m>

    """
    lambda_D = np.sqrt(eps0 * k_B * T_e / (n_e * e ** 2))
    return lambda_D


lambdaD_ = Debye_length
"""Alias to `~plasmapy.formulary.parameters.Debye_length`."""
