"""Functions to calculate fundamental plasma parameters."""

__all__ = [
    "Bohm_diffusion",
    "Debye_number",
    "Hall_parameter",
    "kappa_thermal_speed",
    "magnetic_energy_density",
    "magnetic_pressure",
    "mass_density",
    "thermal_pressure",
]
__aliases__ = [
    "betaH_",
    "DB_",
    "nD_",
    "pmag_",
    "pth_",
    "rho_",
    "ub_",
    "vth_kappa_",
]
__lite_funcs__ = ["thermal_speed_lite"]

import astropy.units as u
import numbers
import numpy as np

from astropy.constants.si import e, eps0, k_B, mu0
from typing import Optional, Union

from plasmapy import particles
from plasmapy.formulary.parameters.frequencies import gyrofrequency
from plasmapy.formulary.parameters.lengths import Debye_length
from plasmapy.formulary.parameters.speeds import thermal_speed
from plasmapy.particles import Particle
from plasmapy.utils.decorators import (
    check_relativistic,
    validate_quantities,
)

__all__ += __aliases__ + __lite_funcs__

e_si_unitless = e.value
eps0_si_unitless = eps0.value
k_B_si_unitless = k_B.value


def _grab_charge(ion: Particle, z_mean=None):
    """
    Merge two possible inputs for particle charge.

    Parameters
    ----------
    ion : `~plasmapy.particles.Particle`
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

    particle : `~plasmapy.particles.Particle`
        The particle for which the mass density is being calculated for.  Must
        be a `~plasmapy.particles.Particle` or a value convertible to
        a `~plasmapy.particles.Particle` (e.g., ``'p'`` for protons,
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
        `~plasmapy.particles.Particle`.

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
"""Alias to `~plasmapy.formulary.parameters.parameters_.mass_density`."""


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
"""Alias to `~plasmapy.formulary.parameters.parameters_.thermal_pressure`."""


@check_relativistic
@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()}
)
def kappa_thermal_speed(
    T: u.K, kappa, particle: Particle, method="most_probable"
) -> u.m / u.s:
    r"""Return the most probable speed for a particle within a Kappa
    distribution.

    **Aliases:** `vth_kappa_`

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        The particle temperature in either kelvin or energy per particle

    kappa: `float`
        The ``kappa`` parameter is a dimensionless number which sets the slope
        of the energy spectrum of suprathermal particles forming the tail
        of the Kappa velocity distribution function. ``kappa`` must be greater
        than 3/2.

    particle : `~plasmapy.particles.Particle`
        Representation of the particle species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4). If no
        charge state information is provided, then the particles are
        assumed to be singly charged.

    method : `str`, optional
        Method to be used for calculating the thermal speed. Options are
        ``'most_probable'`` (default), ``'rms'``, and ``'mean_magnitude'``.

    Returns
    -------
    V : `~astropy.units.Quantity`
        Particle thermal speed.

    Raises
    ------
    `TypeError`
        The particle temperature is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the particle temperature is not in units of temperature or
        energy per particle.

    `ValueError`
        The particle temperature is invalid or particle cannot be used to
        identify an isotope or particle.

    Warns
    -----
    : `~plasmapy.utils.exceptions.RelativityWarning`
        If the particle thermal speed exceeds 5% of the speed of light.

    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    The particle thermal speed is given by:

    .. math::
        V_{th,i} = \sqrt{(2 κ - 3)\frac{2 k_B T_i}{κ m_i}}

    For more discussion on the ``'mean_magnitude'`` calculation method,
    see [1]_.


    Examples
    --------
    >>> from astropy import units as u
    >>> kappa_thermal_speed(5*u.eV, 4, 'p') # defaults to most probable
    <Quantity 24467.87... m / s>
    >>> kappa_thermal_speed(5*u.eV, 4, 'p', 'rms')
    <Quantity 37905.47... m / s>
    >>> kappa_thermal_speed(5*u.eV, 4, 'p', 'mean_magnitude')
    <Quantity 34922.98... m / s>

    References
    ----------
    .. [1] PlasmaPy Issue #186, https://github.com/PlasmaPy/PlasmaPy/issues/186

    See Also
    --------
    ~plasmapy.formulary.kappa_thermal_speed
    ~plasmapy.formulary.kappa_velocity_1D
    """
    # Checking thermal units
    if kappa <= 3 / 2:
        raise ValueError(
            f"Must have kappa > 3/2, instead of {kappa}, for "
            "kappa distribution function to be valid."
        )
    # different methods, as per https://en.wikipedia.org/wiki/Thermal_velocity
    vTh = thermal_speed(T=T, particle=particle, method=method)

    if method == "most_probable":
        # thermal velocity of Kappa distribution function is just Maxwellian
        # thermal speed modulated by the following factor.
        # This is only true for "most probable" case. RMS and mean
        # magnitude velocities are same as Maxwellian.
        coeff = np.sqrt((kappa - 3 / 2) / kappa)
    else:
        coeff = 1

    return vTh * coeff


vth_kappa_ = kappa_thermal_speed
"""Alias to `~plasmapy.formulary.parameters.parameters_.kappa_thermal_speed`."""


@validate_quantities(
    n={"can_be_negative": False},
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
def Hall_parameter(
    n: u.m ** -3,
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

    ion : `~plasmapy.particles.Particle`
        The type of ion ``particle`` is colliding with.

    particle : `~plasmapy.particles.Particle`
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
        (`~plasmapy.formulary.parameters.parameters_.thermal_speed`).

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
    ~plasmapy.formulary.parameters.parameters_.gyrofrequency
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
    >>> from astropy import units as u
    >>> Hall_parameter(1e10 * u.m**-3, 2.8e3 * u.eV, 2.3 * u.T, 'He-4 +1', 'e-')
    <Quantity 7.26446...e+16>
    >>> Hall_parameter(1e10 * u.m**-3, 5.8e3 * u.eV, 2.3 * u.T, 'He-4 +1', 'e-')
    <Quantity 2.11158...e+17>
    """
    from plasmapy.formulary.collisions import (
        fundamental_electron_collision_freq,
        fundamental_ion_collision_freq,
    )

    gyro_frequency = gyrofrequency(B, particle)
    gyro_frequency = gyro_frequency / u.radian
    if particles.Particle(particle).symbol == "e-":
        coll_rate = fundamental_electron_collision_freq(
            T, n, ion, coulomb_log, V, coulomb_log_method=coulomb_log_method
        )
    else:
        coll_rate = fundamental_ion_collision_freq(T, n, ion, coulomb_log, V)
    return gyro_frequency / coll_rate


betaH_ = Hall_parameter
"""Alias to `~plasmapy.formulary.parameters.parameters_.Hall_parameter`."""


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def Debye_number(T_e: u.K, n_e: u.m ** -3) -> u.dimensionless_unscaled:
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
    Debye_length

    Examples
    --------
    >>> from astropy import units as u
    >>> Debye_number(5e6*u.K, 5e9*u.cm**-3)
    <Quantity 2.17658...e+08>

    """

    lambda_D = Debye_length(T_e, n_e)
    N_D = (4 / 3) * np.pi * n_e * lambda_D ** 3

    return N_D


nD_ = Debye_number
"""Alias to `~plasmapy.formulary.parameters.parameters_.Debye_number`."""


@validate_quantities
def magnetic_pressure(B: u.T) -> u.Pa:
    r"""
    Calculate the magnetic pressure.

    **Aliases:** `pmag_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field in units convertible to tesla.

    Returns
    -------
    p_B : `~astropy.units.Quantity`
        The magnetic pressure in units in pascals (newtons per square meter).

    Raises
    ------
    `TypeError`
        If the input is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the input is not in units convertible to tesla.

    `ValueError`
        If the magnetic field strength is not a real number between
        :math:`±∞`\ .

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    The magnetic pressure is given by:

    .. math::
        p_B = \frac{B^2}{2 μ_0}

    The motivation behind having two separate functions for magnetic
    pressure and magnetic energy density is that it allows greater
    insight into the physics that are being considered by the user and
    thus more readable code.

    See Also
    --------
    magnetic_energy_density : returns an equivalent `~astropy.units.Quantity`,
        except in units of joules per cubic meter.

    Examples
    --------
    >>> from astropy import units as u
    >>> magnetic_pressure(0.1*u.T).to(u.Pa)
    <Quantity 3978.87... Pa>

    """
    return (B ** 2) / (2 * mu0)


pmag_ = magnetic_pressure
"""Alias to `~plasmapy.formulary.parameters.parameters_.magnetic_pressure`."""


@validate_quantities
def magnetic_energy_density(B: u.T) -> u.J / u.m ** 3:
    r"""
    Calculate the magnetic energy density.

    **Aliases:** `ub_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field in units convertible to tesla.

    Returns
    -------
    E_B : `~astropy.units.Quantity`
        The magnetic energy density in units of joules per cubic meter.

    Raises
    ------
    `TypeError`
        If the input is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the input is not in units convertible to tesla.

    `ValueError`
        If the magnetic field strength does not have an appropriate.
        value.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed

    Notes
    -----
    The magnetic energy density is given by:

    .. math::
        E_B = \frac{B^2}{2 μ_0}

    The motivation behind having two separate functions for magnetic
    pressure and magnetic energy density is that it allows greater
    insight into the physics that are being considered by the user and
    thus more readable code.

    See Also
    --------
    magnetic_pressure : Returns an equivalent `~astropy.units.Quantity`,
        except in units of pascals.

    Examples
    --------
    >>> from astropy import units as u
    >>> magnetic_energy_density(0.1*u.T)
    <Quantity 3978.87... J / m3>

    """
    return magnetic_pressure(B)


ub_ = magnetic_energy_density
"""Alias to `~plasmapy.formulary.parameters.parameters_.magnetic_energy_density`."""


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    B={"can_be_negative": False},
)
def Bohm_diffusion(T_e: u.K, B: u.T) -> u.m ** 2 / u.s:
    r"""
    Return the Bohm diffusion coefficient.

    The Bohm diffusion coefficient was conjectured to follow Bohm model
    of the diffusion of plasma across a magnetic field and describe the
    diffusion of early fusion energy machines :cite:p:`bohm:1949`. The
    rate predicted by Bohm diffusion is much higher than classical
    diffusion, and if there were no exceptions, magnetically confined
    fusion would be impractical.

    .. math::

        D_B = \frac{1}{16} \frac{k_B T}{e B}

    where :math:`k_B` is the Boltzmann constant
    and :math:`e` is the fundamental charge.

    **Aliases:** `DB_`

    Parameters
    ----------
    T_e : `~astropy.units.Quantity`
        The electron temperature.

    B : `~astropy.units.Quantity`
        The magnitude of the magnetic field in the plasma.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Raises
    ------
    `TypeError`
        ``T_e`` is not a `~astropy.units.Quantity` and cannot be
        converted into one.

    `~astropy.units.UnitConversionError`
        If ``T_e`` is not in appropriate units.

    Examples
    --------
    >>> import astropy.units as u
    >>> T_e = 5000 * u.K
    >>> B = 10 * u.T
    >>> Bohm_diffusion(T_e, B)
    <Quantity 0.00269292 m2 / s>
    >>> T_e = 50 * u.eV
    >>> B = 10 * u.T
    >>> Bohm_diffusion(T_e, B)
    <Quantity 0.3125 m2 / s>

    Returns
    -------
    D_B : `~astropy.units.Quantity`
        The Bohm diffusion coefficient in meters squared per second.

    """
    D_B = k_B * T_e / (16 * e * B)
    return D_B


DB_ = Bohm_diffusion
"""Alias to `~plasmapy.formulary.parameters.parameters_.Bohm_diffusion`."""
