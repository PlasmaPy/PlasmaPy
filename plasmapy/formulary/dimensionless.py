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
    "Debye_number",
    "Hall_parameter",
    "Mag_Reynolds",
    "quantum_theta",
    "Reynolds_number",
]
__aliases__ = ["betaH_", "nD_", "Re_", "Rm_"]

import astropy.units as u
import numpy as np

from astropy.constants.si import k_B, mu0

from plasmapy.formulary import frequencies, lengths, misc, quantum
from plasmapy.particles import Particle
from plasmapy.utils.decorators import validate_quantities

__all__ += __aliases__


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def Debye_number(T_e: u.K, n_e: u.m**-3) -> u.dimensionless_unscaled:
    r"""Return the number of electrons within a sphere with a radius
    of the Debye length.

    **Aliases:** `nD_`

    Parameters
    ----------
    T_e : `~astropy.units.Quantity`
        Electron temperature.

    n_e : `~astropy.units.Quantity`
        Electron number density.

    Raises
    ------
    `TypeError`
        If either argument is not a `~astropy.units.Quantity`.

    `astropy.units.UnitConversionError`
        If either argument is in incorrect units.

    `ValueError`
        If either argument contains invalid values.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Returns
    -------
    N_D : `~astropy.units.Quantity`
        Number of electrons within a sphere with a radius of the Debye
        length.

    Notes
    -----
    The Debye number is the number of electrons contained within a
    sphere with a radius of a Debye length and is given by

    .. math::
        N_D = \frac{4π}{3} n_e λ_D^3

    The Debye number is also known as the plasma parameter.

    Collective behavior requires :math:`N_D ≫ 1`\ .

    See Also
    --------
    ~plasmapy.formulary.lengths.Debye_length

    Examples
    --------
    >>> from astropy import units as u
    >>> Debye_number(5e6*u.K, 5e9*u.cm**-3)
    <Quantity 2.17658...e+08>

    """

    lambda_D = lengths.Debye_length(T_e, n_e)
    return (4 / 3) * np.pi * n_e * lambda_D**3


nD_ = Debye_number
"""Alias to `~plasmapy.formulary.dimensionless.Debye_number`."""


@validate_quantities(
    n={"can_be_negative": False},
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
def Hall_parameter(
    n: u.m**-3,
    T: u.K,
    B: u.T,
    ion: Particle,
    particle: Particle,
    coulomb_log=None,
    V=None,
    coulomb_log_method="classical",
):
    r"""
    Calculate the ``particle`` Hall parameter for a plasma.

    The Hall parameter for plasma species :math:`s` (``particle``) is given by:

    .. math::

        β_{s} = \frac{Ω_{c s}}{ν_{s s^{\prime}}}

    where :math:`Ω_{c s}` is the gyrofrequncy for plasma species :math:`s`
    (``particle``) and :math:`ν_{s s^{\prime}}` is the collision frequency
    between plasma species :math:`s` (``particle``) and species
    :math:`s^{\prime}` (``ion``).

    **Aliases:** `betaH_`

    Parameters
    ----------
    n : `~astropy.units.quantity.Quantity`
        The number density associated with ``particle``.

    T : `~astropy.units.quantity.Quantity`
        The temperature of associated with ``particle``.

    B : `~astropy.units.quantity.Quantity`
        The magnetic field.

    ion : `~plasmapy.particles.particle_class.Particle`
        The type of ion ``particle`` is colliding with.

    particle : `~plasmapy.particles.particle_class.Particle`
        The particle species for which the Hall parameter is calculated for.
        Representation of the particle species (e.g., ``'p'`` for protons,
        ``'D+'`` for deuterium, or ``'He-4 +1'`` for singly ionized helium-4).
        If no charge state information is provided, then the particles are
        assumed to be singly charged.

    coulomb_log : `float`, optional
        Preset value for the Coulomb logarithm. Used mostly for testing purposes.

    V : `~astropy.units.quantity.Quantity`
        The relative velocity between ``particle`` and ``ion``.  If not provided,
        then the ``particle`` thermal velocity is assumed
        (`~plasmapy.formulary.speeds.thermal_speed`).

    coulomb_log_method : `str`, optional
        The method by which to compute the Coulomb logarithm.
        The default method is the classical straight-line Landau-Spitzer
        method (``"classical"`` or ``"ls"``). The other 6 supported methods
        are ``"ls_min_interp"``, ``"ls_full_interp"``, ``"ls_clamp_mininterp"``,
        ``"hls_min_interp"``, ``"hls_max_interp"``, and ``"hls_full_interp"``.
        Please refer to the docstring of
        `~plasmapy.formulary.collisions.Coulomb_logarithm` for more
        information about these methods.

    See Also
    --------
    ~plasmapy.formulary.frequencies.gyrofrequency
    ~plasmapy.formulary.collisions.fundamental_electron_collision_freq
    ~plasmapy.formulary.collisions.fundamental_ion_collision_freq
    ~plasmapy.formulary.collisions.Coulomb_logarithm

    Returns
    -------
    `~astropy.units.quantity.Quantity`
        Hall parameter for ``particle``.

    Notes
    -----
    * For calculating the collision frequency
      `~plasmapy.formulary.collisions.fundamental_electron_collision_freq` is used
      when ``particle`` is an electron and
      `~plasmapy.formulary.collisions.fundamental_ion_collision_freq` when
      ``particle`` is an ion.
    * The collision frequencies are calculated assuming a slowly moving
      Maxwellian distribution.

    Examples
    --------
    >>> import astropy.units as u
    >>> import pytest
    >>> from plasmapy.utils.exceptions import RelativityWarning

    >>> Hall_parameter(1e10 * u.m**-3, 2.8e2 * u.eV, 2.3 * u.T, 'He-4 +1', 'e-')
    <Quantity 2.500...e+15>
    >>> with pytest.warns(RelativityWarning):
    ...     Hall_parameter(1e10 * u.m**-3, 5.8e3 * u.eV, 2.3 * u.T, 'He-4 +1', 'e-')
    <Quantity 2.11158...e+17>
    """
    from plasmapy.formulary.collisions import (
        fundamental_electron_collision_freq,
        fundamental_ion_collision_freq,
    )

    gyro_frequency = frequencies.gyrofrequency(B, particle)
    gyro_frequency = gyro_frequency / u.radian
    if Particle(particle).symbol == "e-":
        coll_rate = fundamental_electron_collision_freq(
            T, n, ion, coulomb_log, V, coulomb_log_method=coulomb_log_method
        )
    else:
        coll_rate = fundamental_ion_collision_freq(T, n, ion, coulomb_log, V)
    return gyro_frequency / coll_rate


betaH_ = Hall_parameter
"""Alias to `~plasmapy.formulary.dimensionless.Hall_parameter`."""


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def quantum_theta(T: u.K, n_e: u.m**-3) -> u.dimensionless_unscaled:
    r"""
    Compare Fermi energy to thermal kinetic energy to check if quantum
    effects are important.

    The quantum theta (:math:`θ`) of a plasma is defined by

    .. math::
        θ = \frac{E_T}{E_F}

    where :math:`E_T` is the thermal energy of the plasma
    and :math:`E_F` is the Fermi energy of the plasma.

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        The temperature of the plasma.

    n_e : `~astropy.units.Quantity`
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
    theta : `~astropy.units.Quantity`

    Notes
    -----
    The thermal energy of the plasma (:math:`E_T`) is defined by

    .. math::
        E_T = k_B T

    where :math:`k_B` is the Boltzmann constant
    and :math:`T` is the temperature of the plasma.

    See Also
    --------
    ~plasmapy.formulary.quantum.Fermi_energy
    """
    fermi_energy = quantum.Fermi_energy(n_e)
    thermal_energy = k_B * T
    return thermal_energy / fermi_energy


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n={"can_be_negative": False},
)
def beta(T: u.K, n: u.m**-3, B: u.T) -> u.dimensionless_unscaled:
    r"""
    Compute the ratio of thermal pressure to magnetic pressure.

    The beta (:math:`β`) of a plasma is defined by

    .. math::
        β = \frac{p_{th}}{p_{mag}}

    where :math:`p_{th}` is the thermal pressure of the plasma
    and :math:`p_{mag}` is the magnetic pressure of the plasma.

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        The temperature of the plasma.

    n : `~astropy.units.Quantity`
        The particle density of the plasma.

    B : `~astropy.units.Quantity`
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
    beta: `~astropy.units.Quantity`
        Dimensionless quantity.

    See Also
    --------
    ~plasmapy.formulary.misc.thermal_pressure
    ~plasmapy.formulary.misc.magnetic_pressure
    """
    thermal_pressure = misc.thermal_pressure(T, n)
    magnetic_pressure = misc.magnetic_pressure(B)
    return thermal_pressure / magnetic_pressure


@validate_quantities(U={"can_be_negative": True})
def Reynolds_number(
    rho: u.kg / u.m**3, U: u.m / u.s, L: u.m, mu: u.kg / (u.m * u.s)
) -> u.dimensionless_unscaled:
    r"""
    Compute the Reynolds number.

    The Reynolds number is a dimensionless quantity that is used to
    predict flow patterns in fluids. The Reynolds number is defined as
    the ratio of inertial forces to viscous forces. A low Reynolds
    number describes smooth, laminar flow while a high Reynolds number
    describes rough, turbulent flow.

    .. math::

        Re = \frac{ρ U L}{μ}

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
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Raises
    ------
    `TypeError`
        If ``U`` is not a `~astropy.units.Quantity` and cannot be
        converted into a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If ``U`` is not in appropriate units.

    :exc:`~plasmapy.utils.exceptions.RelativityError`
        If ``U`` is greater than the speed of light.

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
    Re: `~astropy.units.Quantity`
        Dimensionless quantity.

    """
    return abs(rho * U * L / mu)


Re_ = Reynolds_number
"""Alias to `~plasmapy.formulary.dimensionless.Reynolds_number`."""


@validate_quantities(U={"can_be_negative": True})
def Mag_Reynolds(U: u.m / u.s, L: u.m, sigma: u.S / u.m) -> u.dimensionless_unscaled:
    r"""
    Compute the magnetic Reynolds number

    The magnetic Reynolds number is a dimensionless quantity that
    estimates the relative contributions of advection and induction
    to magnetic diffusion in a conducting medium.

    .. math::

        Rm = \frac{U L}{η}

    where :math:`η = \frac{1}{μ_0 σ}`
    and :math:`μ_0` is the permeability of free space.

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
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Raises
    ------
    `TypeError`
        If ``U`` is not a `~astropy.units.Quantity` and cannot be
        converted into a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If ``U`` is not in appropriate units.

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
    Rm : `~astropy.units.Quantity`
        The magnetic Reynolds number.

    """
    eta = 1 / (mu0 * sigma)
    return abs(U * L / eta)


Rm_ = Mag_Reynolds
"""Alias to `~plasmapy.formulary.dimensionless.Mag_Reynolds`."""
