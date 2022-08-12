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

    Notes
    -----
    The plasma critical density is given by the formula

    .. math::
        n_{c}=\frac{m_{e}\varepsilon_0\omega_{p}^{2}}{e^{2}}

    where :math:`m_{e}` is the mass of an electron,
    :math:`\varepsilon_0` is the permittivity of free space, :math:`\omega_{p}`
    is the plasma frequency, and :math:`e` is the elementary charge.

    Examples
    --------
    >>> from astropy import units as u
    >>> plasma_critical_density(1 * u.rad/u.s)
    <Quantity 0.00031421 1 / m3>

    """

    n_c = m_e * eps0 * omega_p**2 / (e**2)

    return (n_c/u.rad**2).to(u.m**-3)
