"""Functions to calculate plasma density parameters."""
__all__ = [
    "plasma_critical_density",
]

import astropy.units as u

from astropy.constants.si import e, eps0, m_e

from plasmapy.utils.decorators import validate_quantities


@validate_quantities(
    omega_p={"can_be_negative": False},
    validations_on_return={
        "units": [u.m**-3],
        "equivalencies": [
            (u.F * u.kg * u.rad**2 / (u.C**2 * u.m * u.s**2), u.m**-3)
        ],
    },
)
def plasma_critical_density(omega_p: u.rad / u.s) -> u.m**-3:
    r"""Calculate the plasma critical density.

    Parameters
    ----------
    omega_p : `~astropy.units.Quantity`
        The plasma frequency in units of angular frequency.

    Returns
    -------
    n_c : `~astropy.units.Quantity`
        The plasma critical density.

    Raises
    ------
    `TypeError`
        If ``omega_p`` is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitTypeError`
        If ``omega_p`` is not in correct units.

    `ValueError`
        If ``omega_p`` is an invalid value.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----


    Examples
    --------
    >>> from astropy import units as u
    >>> plasma_critical_density(1 * u.rad/u.s)
    <Quantity 0.00031421 1 / m3>

    """

    return m_e * eps0 * omega_p**2 / (e**2)
