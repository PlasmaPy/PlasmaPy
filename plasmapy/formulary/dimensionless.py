"""
Module of dimensionless plasma parameters.

These are especially important for determining what regime a plasma is
in. (e.g., turbulent, quantum, collisional, etc.).

For example, plasmas at high (much larger than 1) Reynolds numbers are
highly turbulent, while turbulence is negligible at low Reynolds
numbers.
"""
__all__ = ["beta", "quantum_theta", "magnetic_prandtl_number"]

from astropy import constants
from astropy import units as u
from plasmapy.formulary import quantum, parameters
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


@validate_quantities(
    v = {"can_be_negative": False},
    alpha = {"can_be_negative": False},
)
def prandtl_number(v: u.m ** 2 / u.s,
                   alpha: u.m ** 2 / u.s) -> u.dimensionless_unscaled:
    r"""
    The ratio of momentum diffusivity to thermal diffusivity..

    Parameters
    ----------
    v : ~astropy.units.Quantity
        The momentum diffusivity (kinematic viscosity).
    alpha : ~astropy.units.Quantity
        The thermal diffusivity.

    Returns
    -------
    prandtl_number: ~astropy.units.Quantity
        Dimensionless quantity.

    Notes
    -----
    The Prandtl number is given as [1]_:

    .. math::
        \mathrm{Pr} =  \frac{\nu}{\alpha}

    Examples
    --------
    >>> import astropy.units as u
    >>> prandtl_number(8e-2*u.m**2/u.s, 3e2*u.m**2/u.s)
    <Quantity 0.00026667>
    >>> prandtl_number(4e2*u.m**2/u.s,2e-10*u.m**2/u.s)
    <Quantity 2.e+12>


    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Prandtl_number

    """

    return v / alpha


@validate_quantities(
    v = {"can_be_negative": False},
    n = {"can_be_negative": False},
)
def magnetic_prandtl_number(v: u.m ** 2 / u.s,
                            n: u.m ** 2 / u.s) -> u.dimensionless_unscaled:
    r"""
    The ratio of momentum diffusivity (viscosity) and magnetic diffusivity.

    Parameters
    ----------
    v : ~astropy.units.Quantity
        The momentum diffusivity (kinematic viscosity).
    n : ~astropy.units.Quantity
        The magnetic diffusivity.

    Returns
    -------
    magnetic_prandtl_number: ~astropy.units.Quantity
        Dimensionless quantity.

    Notes
    -----
    The magnetic Prandtl number is given as [1]_:

    .. math::
        \mathrm{Pr}_\mathrm{m} =  \frac{\nu}{\eta}

    Examples
    --------
    >>> import astropy.units as u
    >>> magnetic_prandtl_number(8e-2*u.m**2/u.s, 3e2*u.m**2/u.s)
    <Quantity 0.00026667>
    >>> magnetic_prandtl_number(4e2*u.m**2/u.s,2e-10*u.m**2/u.s)
    <Quantity 2.e+12>


    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Magnetic_Prandtl_number

    """

    return v / n
