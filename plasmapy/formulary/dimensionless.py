"""
Module of dimensionless plasma parameters.

These are especially important for determining what regime a plasma is
in. (e.g., turbulent, quantum, collisional, etc.).

For example, plasmas at high (much larger than 1) Reynolds numbers are
highly turbulent, while turbulence is negligible at low Reynolds
numbers.
"""
__all__ = ['beta', 'quantum_theta']

from astropy import constants
from astropy import units as u
from plasmapy.formulary import (quantum, parameters)
from plasmapy.utils.decorators import validate_quantities


@validate_quantities(T={'can_be_negative': False,
                        'equivalencies': u.temperature_energy()},
                     n_e={'can_be_negative': False})
def quantum_theta(T: u.K, n_e: u.m**-3) -> u.dimensionless_unscaled:
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
    <Quantity 59083071.8...>
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


@validate_quantities(T={'can_be_negative': False,
                        'equivalencies': u.temperature_energy()},
                     n={'can_be_negative': False})
def beta(T: u.K, n: u.m**-3, B: u.T) -> u.dimensionless_unscaled:
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
