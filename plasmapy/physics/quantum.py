import numpy as np
from astropy import units
from ..constants import c, h, ion_mass
from ..utils import _check_quantity, _check_relativistic
from .relativity import Lorentz_factor


def deBroglie_wavelength(V, particle):
    """Calculates the de Broglie wavelength.

    Parameters
    ----------
    V : Quantity
        Particle velocity in units convertible to meters per second.

    particle : string or Quantity
        Representation of the particle species (e.g., 'e', 'p', 'D+',
        or 'He-4 1+', or the particle mass in units convertible to
        kilograms.

    Returns
    -------
    lambda_dB : Quantity
        The de Broglie wavelength in units of meters.

    Raises
    ------
    TypeError
        The velocity is not a Quantity and cannot be converted into a
        Quantity.

    UnitConversionError
        If the velocity is not in appropriate units.

    ValueError
        If the magnitude of V is faster than the speed of light.

    UserWarning
        If V is not a Quantity, then a UserWarning will be raised and
        units of meters per second will be assumed.

    Notes
    -----
    The de Broglie wavelength is given by

    .. math::
    \lambda_{dB} = \frac{h}{p} = \frac{h}{\gamma m V}.

    where :math:`h` is the Planck constant, :math:`p` is the
    relativistic momentum of the particle, :math:`gamma` is the
    Lorentz factor, `m` is the particle's mass, and :math:`V` is the
    particle's velocity.

    Examples
    --------
    >>> from astropy import units as u
    >>> velocity = 1.4e7*u.m/u.s
    >>> deBroglie_wavelength(velocity, 'e')
    <Quantity 5.1899709519786425e-11 m>
    >>> deBroglie_wavelength(V = 0*u.m/u.s, particle = 'D+')
    <Quantity inf m>

    """

    _check_quantity(V, 'V', 'deBroglie_wavelength', units.m/units.s)

    V = np.abs(V)

    if np.any(V >= c):
        raise ValueError("Velocity input in deBroglie_wavelength cannot be "
                         "greater than or equal to the speed of light.")

    try:
        m = ion_mass(particle)  # Replace with more general routine!
    except Exception:
        raise ValueError("Unable to find particle mass.")

    if V.size > 1:

        is_zero = np.where(V == 0*units.m/units.s)
        is_nonzero = np.invert(is_zero)

        lambda_dB = np.zeros_like(V.value)
        lambda_dB[is_zero] = np.inf*units.m
        lambda_dB[nonzero] = h / (m * V[nonzero] * Lorentz_factor(V[nonzero]))

    else:

        if V == 0*units.m/units.s:
            lambda_dB = np.inf*units.m
        else:
            lambda_dB = h / (Lorentz_factor(V) * m * V)

    return lambda_dB.to(units.m)
