"""
Functions for quantum parameters, including electron degenerate
gases and warm dense matter.
"""
import numpy as np
from astropy import units
from ..constants import c, h, hbar, m_e, eps0, e, k_B
from ..atomic import ion_mass
from ..utils import _check_quantity, _check_relativistic, check_quantity
from .relativity import Lorentz_factor


def deBroglie_wavelength(V, particle):
    r"""Calculates the de Broglie wavelength.

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

    if not isinstance(particle, units.Quantity):
        try:
            m = ion_mass(particle)  # TODO: Replace with more general routine!
        except Exception:
            raise ValueError("Unable to find particle mass.")
    else:
        try:
            m = particle.to(units.kg)
        except Exception:
            raise units.UnitConversionError("The second argument for deBroglie"
                                            " wavelength must be either a "
                                            "representation of a particle or a"
                                            " Quantity with units of mass.")

    if V.size > 1:

        lambda_dBr = np.ones(V.shape) * np.inf * units.m
        indices = V.value != 0
        lambda_dBr[indices] = h / (m*V[indices]*Lorentz_factor(V[indices]))

    else:

        if V == 0*units.m/units.s:
            lambda_dBr = np.inf*units.m
        else:
            lambda_dBr = h / (Lorentz_factor(V) * m * V)

    return lambda_dBr.to(units.m)


@check_quantity({
    'T_e': {'units': units.K, 'can_be_negative': False}
})
def thermal_deBroglie_wavelength(T_e):
    r"""Calculate the thermal deBroglie wavelength for electrons.

    Parameters
    ----------
    T_e: Quantity
        Electron temperature

    Returns
    -------
    lambda_dbTh: Quantity
        The thermal deBroglie wavelength for electrons in meters

    Raises
    ------
    TypeError
        If argument is not a Quantity

    UnitConversionError
        If argument is in incorrect units

    ValueError
        If argument contains invalid values

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----
    The thermal deBroglie wavelength is approximately the average deBroglie
    wavelength for electrons in an ideal gas and is given by

    .. math::
    \lambda_dbTh = \frac{h}{\sqrt{2 \pi m_e k_B T_e}}

    See also
    --------
    

    Example
    -------
    >>> from astropy import units as u
    >>> thermal_deBroglie_wavelength(1 * u.eV)
    <Quantity 6.919367518364532e-10 m>

    """
    T_e = T_e.to(units.K, equivalencies=units.temperature_energy())
    try:
        lambda_dbTh = h / np.sqrt(2 * np.pi * m_e * k_B * T_e)
    except Exception:
        raise ValueError("Unable to find thermal deBroglie wavelength.")
    return lambda_dbTh.to(units.m)

@check_quantity({
    'n_e': {'units': units.m**-3, 'can_be_negative': False}
})
def Fermi_energy(n_e):
    r"""Calculate the Fermi energy.

    Parameters
    ----------
    n_e: Quantity
        Electron number density

    Returns
    -------
    energy_F: Quantity
        The Fermi energy in Joules

    Raises
    ------
    TypeError
        If argument is not a Quantity

    UnitConversionError
        If argument is in incorrect units

    ValueError
        If argument contains invalid values

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----
    The Fermi energy is the kinetic energy in a degenerate electron gas
    and is given by

    .. math::
    E_F = \frac{\pi^2 \hbar^2}{2 m_e}\left(\frac{3 n_e}{\pi}\right)^{2/3}

    This quantity is often used in place of thermal energy for analysis
    of cold, dense plasmas (e.g. warm dense matter, condensed matter).

    See also
    --------
    Thomas_Fermi_length

    Example
    -------
    >>> from astropy import units as u
    >>> Fermi_energy(1e23 * u.cm**-3)
    <Quantity 1.2586761116196002e-18 J>

    """
    try:
        coeff = (np.pi * hbar) ** 2 / (2 * m_e)
        energy_F = coeff * (3 * n_e / np.pi) ** (2/3)
    except Exception:
        raise ValueError("Unable to find Fermi energy.")
    return energy_F.to(units.Joule)

@check_quantity({
    'n_e': {'units': units.m**-3, 'can_be_negative': False}
})
def Thomas_Fermi_length(n_e):
    r"""Calculate the Thomas-Fermi screening length.

    Parameters
    ----------
    n_e: Quantity
        Electron number density

    Returns
    -------
    lambda_TF: Quantity
        The Thomas-Fermi screening length in meters

    Raises
    ------
    TypeError
        If argument is not a Quantity

    UnitConversionError
        If argument is in incorrect units

    ValueError
        If argument contains invalid values

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----
    The Thomas-Fermi screening length is the exponential scale length for 
    charge screening and is given by

    .. math::
    \lambda_TF = \sqrt{\frac{2 \epsilon_0 E_F}{3 n_e e^2}}

    for an electron degenerate gas.
    
    This quantity is often used in place of the Debye length for analysis
    of cold, dense plasmas (e.g. warm dense matter, condensed matter).

    The electrical potential will drop by a factor of 1/e every Thomas-Fermi
    screening length.

    Plasmas will generally be quasineutral on length scales significantly
    larger than the Thomas-Fermi screening length.

    See also
    --------
    Fermi_energy

    Example
    -------
    >>> from astropy import units as u
    >>> Thomas_Fermi_length(1e23 * u.cm**-3)
    <Quantity 5.379914085596706e-11 m>

    """
    try:
        energy_F = Fermi_energy(n_e)
        lambda_TF = np.sqrt(2 * eps0 * energy_F / (3 * n_e * e ** 2))
    except Exception:
        raise ValueError("Unable to find Thomas Fermi screening length.")
    
    return lambda_TF.to(units.m)