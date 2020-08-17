"""
Module of dimensionless plasma parameters.

These are especially important for determining what regime a plasma is
in. (e.g., turbulent, quantum, collisional, etc.).

For example, plasmas at high (much larger than 1) Reynolds numbers are
highly turbulent, while turbulence is negligible at low Reynolds
numbers.
"""
__all__ = [
    "beta",
    "Mag_Reynolds",
    "quantum_theta",
    "Re_",
    "Reynolds_number",
    "Rm_",
]

from astropy import constants
from astropy import units as u
from astropy.constants import c
from astropy.constants.codata2010 import mu0

from plasmapy.formulary import electron_viscosity, parameters, quantum
from plasmapy.utils.decorators import validate_quantities


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def quantum_theta(T: u.K, n_e: u.m ** -3) -> u.dimensionless_unscaled:
    """
    Compares Fermi energy to thermal kinetic energy to check if quantum
    effects are important.

    Parameters
    ----------
    T : ~astropy.units.Quantity
        The temperature of the plasma.
    n_e : ~astropy.units.Quantity
          The electron number density of the plasma.

    Examples
    --------
    >>> import astropy.units as u
    >>> quantum_theta(1*u.eV, 1e20*u.m**-3)
    <Quantity 127290.619...>
    >>> quantum_theta(1*u.eV, 1e16*u.m**-3)
    <Quantity 59083071...>
    >>> quantum_theta(1*u.eV, 1e26*u.m**-3)
    <Quantity 12.72906...>
    >>> quantum_theta(1*u.K, 1e26*u.m**-3)
    <Quantity 0.00109...>

    Returns
    -------
    theta: ~astropy.units.Quantity

    """
    fermi_energy = quantum.Fermi_energy(n_e)
    thermal_energy = constants.k_B * T
    theta = thermal_energy / fermi_energy
    return theta


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n={"can_be_negative": False},
)
def beta(T: u.K, n: u.m ** -3, B: u.T) -> u.dimensionless_unscaled:
    """
    The ratio of thermal pressure to magnetic pressure.

    Parameters
    ----------
    T : ~astropy.units.Quantity
        The temperature of the plasma.
    n : ~astropy.units.Quantity
        The particle density of the plasma.
    B : ~astropy.units.Quantity
        The magnetic field in the plasma.

    Examples
    --------
    >>> import astropy.units as u
    >>> beta(1*u.eV, 1e20*u.m**-3, 1*u.T)
    <Quantity 4.0267...e-05>
    >>> beta(8.8e3*u.eV, 1e20*u.m**-3, 5.3*u.T)
    <Quantity 0.01261...>

    Returns
    -------
    beta: ~astropy.units.Quantity
        Dimensionless quantity.

    """
    thermal_pressure = parameters.thermal_pressure(T, n)
    magnetic_pressure = parameters.magnetic_pressure(B)
    return thermal_pressure / magnetic_pressure


@validate_quantities(U={"can_be_negative": True})
def Reynolds_number(
    rho: u.kg / u.m ** 3, U: u.m / u.s, L: u.m, mu: u.kg / (u.m * u.s)
) -> u.dimensionless_unscaled:
    r"""
    The Reynolds Number is a dimensionless quantity
    that is used to predict flow patterns in fluids.
    The Reynolds Number is defined as the ratio of inertial forces to viscous forces.
    A low Reynolds Number describes smooth, laminar flow
    while a high Reynolds Number describes rough, turbulent flow.

    .. math::

        Re = \frac{\rho U L }{\mu}

    **Aliases:** `Re_`

    Parameters
    ----------
    rho : `~astropy.units.Quantity`
        The density of the plasma.
    U : `~astropy.units.Quantity`
        The flow velocity of the plasma.
    L : `~astropy.units.Quantity`
        The characteristic length scale.
    mu : `~astropy.units.Quantity`
        The dynamic viscosity of the plasma.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed.

    Raises
    ------
    TypeError
        The `U` is not a `~astropy.units.Quantity` and cannot be
        converted into a ~astropy.units.Quantity.

    ~astropy.units.UnitConversionError
        If the `U` is not in appropriate units.

    :exc:`~plasmapy.utils.exceptions.RelativityError`
        If the velocity `U` is greater than the speed of light.

    Examples
    --------
    >>> import astropy.units as u
    >>> rho = 1000 * u.kg / u.m ** 3
    >>> U = 10 * u.m / u.s
    >>> L = 1 * u.m
    >>> mu = 8.9e-4 * u.kg / (u.m * u.s)
    >>> Reynolds_number(rho, U, L, mu)
    <Quantity 11235955.05617978>
    >>> rho = 1490 * u.kg / u.m ** 3
    >>> U = 0.1 * u.m / u.s
    >>> L = 0.05 * u.m
    >>> mu = 10 * u.kg / (u.m * u.s)
    >>> Reynolds_number(rho, U, L, mu)
    <Quantity 0.745>

    Returns
    -------
    Re: `~astropy.Quantity`
        Dimensionless quantity.

    """
    Re = abs(rho * U * L / mu)
    return Re


Re_ = Reynolds_number
""" Alias to :func:`Reynolds_number`. """


@validate_quantities(U={"can_be_negative": True})
def Mag_Reynolds(U: u.m / u.s, L: u.m, sigma: u.S / u.m) -> u.dimensionless_unscaled:
    r"""
    The Magnetic Reynolds number is a dimensionless quantity that
    estimates the relative contributions of advection and induction
    to magnetic diffusion in a conducting medium.

    .. math::

        Rm = \frac{U L}{\eta}

    where :math:`\eta = \frac{1}{\mu_0 \sigma}`
    and :math:`\mu_0` is the permeability of free space.

    **Aliases:** `Rm_`

    Parameters
    ----------
    U : `~astropy.units.Quantity`
        The velocity scale of the plasma.
    L : `~astropy.units.Quantity`
        The length scale of the plasma.
    sigma : `~astropy.units.Quantity`
        The conductivity of the plasma.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed.

    Raises
    ------
    TypeError
        The `U` is not a `~astropy.units.Quantity` and cannot be
        converted into a ~astropy.units.Quantity.

    ~astropy.units.UnitConversionError
        If the `U` is not in appropriate units.



    Examples
    --------
    >>> import astropy.units as u
    >>> sigma = 5.96e7 * u.S / u.m
    >>> U = 10 * u.m / u.s
    >>> L = 1 * u.cm
    >>> Mag_Reynolds(U, L, sigma)
    <Quantity 7.48955689>
    >>> rho = 1e-8 * u.S / u.m
    >>> U = 0.1 * u.m / u.s
    >>> L = 0.05 * u.m
    >>> Mag_Reynolds(U, L, sigma)
    <Quantity 0.37447784>

    Returns
    -------
    Rm: `~astropy.Quantity`
        Dimensionless quantity.

    """
    eta = 1 / (mu0 * sigma)
    Rm = abs(U * L / eta)
    return Rm


Rm_ = Mag_Reynolds
""" Alias to :func:`Mag_Reynolds`. """
