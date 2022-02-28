"""Functions to for miscellaneous plasma parameter calculations."""

__all__ = ["mass_density", "thermal_pressure"]
__aliases__ = ["rho_", "pth_"]

import astropy.units as u
import numbers
import numpy as np

from astropy.constants.si import k_B
from typing import Optional, Union

from plasmapy import particles
from plasmapy.particles import Particle
from plasmapy.utils.decorators import validate_quantities


def _grab_charge(ion: Particle, z_mean=None):
    """
    Merge two possible inputs for particle charge.

    Parameters
    ----------
    ion : `~plasmapy.particles.particle_class.Particle`
        a string representing a charged particle, or a Particle object.

    z_mean : `float`
        An optional float describing the average ionization of a particle
        species.

    Returns
    -------
    float
        if ``z_mean`` was passed, ``z_mean``, otherwise, the charge number
        of ``ion``.

    """
    if z_mean is None:
        # warnings.warn("No z_mean given, defaulting to atomic charge",
        #               PhysicsWarning)
        Z = particles.charge_number(ion)
    else:
        # using average ionization provided by user
        Z = z_mean
    return Z


@validate_quantities(
    density={"can_be_negative": False}, validations_on_return={"can_be_negative": False}
)
def mass_density(
    density: (u.m ** -3, u.kg / (u.m ** 3)),
    particle: Union[Particle, str],
    z_ratio: Optional[numbers.Real] = 1,
) -> u.kg / u.m ** 3:
    r"""
    Calculate the mass density from a number density.

    .. math::

        \rho = \left| \frac{Z_{s}}{Z_{particle}} \right| n_{s} m_{particle}
              = | Z_{ratio} | n_{s} m_{particle}

    where :math:`m_{particle}` is the particle mass, :math:`n_{s}` is a number
    density for plasma species :math:`s`, :math:`Z_{s}` is the charge number of
    species :math:`s`, and :math:`Z_{particle}` is the charge number of
    ``particle``.  For example, if the electron density is given for :math:`n_s`
    and ``particle`` is a doubly ionized atom, then :math:`Z_{ratio} = -1 / 2`\ .

    **Aliases:** `rho_`

    Parameters
    ----------
    density : `~astropy.units.Quantity`
        Either a particle number density (in units of m\ :sup:`-3` or
        equivalent) or a mass density (in units of kg / m\ :sup:`3` or
        equivalent).  If ``density`` is a mass density, then it will be passed
        through and returned without modification.

    particle : `~plasmapy.particles.particle_class.Particle`
        The particle for which the mass density is being calculated for.  Must
        be a `~plasmapy.particles.particle_class.Particle` or a value convertible to
        a `~plasmapy.particles.particle_class.Particle` (e.g., ``'p'`` for protons,
        ``'D+'`` for deuterium, or ``'He-4 +1'`` for singly ionized helium-4).

    z_ratio : `int`, `float`, optional
        The ratio of the charge numbers corresponding to the plasma species
        represented by ``density`` and the ``particle``.  For example, if the
        given ``density`` is and electron density and ``particle`` is doubly
        ionized ``He``, then ``z_ratio = -0.5``.  Default is ``1``.

    Raises
    ------
    `~astropy.units.UnitTypeError`
        If the ``density`` does not have units equivalent to a number density
        or mass density.

    `TypeError`
        If ``density`` is not of type `~astropy.units.Quantity`, or convertible.

    `TypeError`
        If ``particle`` is not of type or convertible to
        `~plasmapy.particles.particle_class.Particle`.

    `TypeError`
        If ``z_ratio`` is not of type `int` or `float`.

    `ValueError`
        If ``density`` is negative.

    Returns
    -------
    `~astropy.units.Quantity`
        The mass density for the plasma species represented by ``particle``.

    Examples
    --------
    >>> import astropy.units as u
    >>> mass_density(1 * u.m ** -3, 'p')
    <Quantity 1.67262...e-27 kg / m3>
    >>> mass_density(4 * u.m ** -3, 'D+')
    <Quantity 1.33743...e-26 kg / m3>
    >>> mass_density(2.e12 * u.cm ** -3, 'He')
    <Quantity 1.32929...e-08 kg / m3>
    >>> mass_density(2.e12 * u.cm ** -3, 'He', z_ratio=0.5)
    <Quantity 6.64647...e-09 kg / m3>
    >>> mass_density(1.0 * u.g * u.m ** -3, "")
    <Quantity 0.001 kg / m3>
    """
    if density.unit.is_equivalent(u.kg / u.m ** 3):
        return density

    if not isinstance(particle, Particle):
        try:
            particle = Particle(particle)
        except TypeError:
            raise TypeError(
                f"If passing a number density, you must pass a plasmapy Particle "
                f"(not type {type(particle)}) to calculate the mass density!"
            )

    if not isinstance(z_ratio, (float, np.floating, int, np.integer)):
        raise TypeError(
            f"Expected type int or float for keyword z_ratio, got type {type(z_ratio)}."
        )

    return abs(z_ratio) * density * particle.mass


rho_ = mass_density
"""Alias to `~plasmapy.formulary.misc.mass_density`."""


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n={"can_be_negative": False},
)
def thermal_pressure(T: u.K, n: u.m ** -3) -> u.Pa:
    r"""
    Return the thermal pressure for a Maxwellian distribution.

    **Aliases:** `pth_`

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        The particle temperature in either kelvin or energy per particle.

    n : `~astropy.units.Quantity`
        The particle number density in units convertible to m\ :sup:`-3`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> thermal_pressure(1*u.eV, 1e20/u.m**3)
    <Quantity 16.021... Pa>
    >>> thermal_pressure(10*u.eV, 1e20/u.m**3)
    <Quantity 160.21... Pa>

    Returns
    -------
    p_th : `~astropy.units.Quantity`
        Thermal pressure.

    Raises
    ------
    `TypeError`
        The temperature or number density is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the particle temperature is not in units of temperature or
        energy per particle.

    Notes
    -----
    The thermal pressure is given by:

    .. math::
        T_{th} = n k_B T
    """
    return n * k_B * T


pth_ = thermal_pressure
"""Alias to `~plasmapy.formulary.misc.thermal_pressure`."""
