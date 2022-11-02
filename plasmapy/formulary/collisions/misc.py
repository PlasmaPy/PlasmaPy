"""
Module of miscellaneous parameters related to collisions.
"""
__all__ = [
    "mobility",
    "Spitzer_resistivity",
]

import astropy.units as u
import numpy as np

from astropy.constants.si import e
from numbers import Real

from plasmapy import particles, utils
from plasmapy.formulary.collisions import frequencies
from plasmapy.formulary.speeds import thermal_speed
from plasmapy.utils.decorators import validate_quantities
from plasmapy.utils.decorators.checks import _check_relativistic


@validate_quantities(T={"equivalencies": u.temperature_energy()})
@particles.particle_input
def _process_inputs(T: u.K, species: (particles.Particle, particles.Particle), V):
    """
    Helper function for processing inputs to functionality contained
    in `plasmapy.formulary.collisions`.

    Also obtains the reduced mass in a 2 particle collision system
    along with thermal velocity.
    """
    masses = [p.mass for p in species]
    charges = [np.abs(p.charge) for p in species]

    # obtaining reduced mass of 2 particle collision system
    reduced_mass = particles.reduced_mass(*species)

    # getting thermal velocity of system if no velocity is given
    V = _replace_nan_velocity_with_thermal_velocity(V, T, reduced_mass)

    _check_relativistic(V, "V")

    return T, masses, charges, reduced_mass, V


# TODO: Remove redundant mass parameter
def _replace_nan_velocity_with_thermal_velocity(
    V, T, m, species=particles.Particle("e-")
):
    """
    Get thermal velocity of system if no velocity is given, for a given
    mass.  Handles vector checks for ``V``, you must already know that
    ``T`` and ``m`` are okay.
    """
    if np.any(V == 0):
        raise utils.PhysicsError("Collisions are not possible with a zero velocity.")

    if V is None:
        return thermal_speed(T, species, mass=m)

    if not np.any(np.isnan(V)):
        return V

    if not (np.isscalar(V.value) or np.isscalar(T.value)):
        V = V.copy()
        V[np.isnan(V)] = thermal_speed(T[np.isnan(V)], species, mass=m)
        return V

    if np.isscalar(V.value):
        return thermal_speed(T, species, mass=m)

    if np.isscalar(T.value):
        V[np.isnan(V)] = thermal_speed(T, species, mass=m)

    return V


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def mobility(
    T: u.K,
    n_e: u.m**-3,
    species,
    z_mean: Real = np.nan,
    V: u.m / u.s = np.nan * u.m / u.s,
    method="classical",
) -> u.m**2 / (u.V * u.s):
    r"""
    Return the electrical mobility.

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and the
        target particle.

    n_e : `~astropy.units.Quantity`
        The electron number density in units convertible to m\ :sup:`-3`.

    species : `tuple`
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second).

    z_mean : `~astropy.units.Quantity`, optional
        The average ionization (arithmetic mean) of a plasma for which a
        macroscopic description is valid. This parameter is used to
        compute the average ion density (given the average ionization
        and electron density) for calculating the ion sphere radius for
        non-classical impact parameters. It is also used to obtain the
        average mobility of a plasma with multiple charge state
        species. When ``z_mean`` is not given, the average charge
        between the two particles is used instead. ``z_mean`` is a
        required parameter if ``method`` is ``"ls_full_interp"``,
        ``"hls_max_interp"``, or ``"hls_full_interp"``.

    V : `~astropy.units.Quantity`, optional
        The relative velocity between particles. If not provided,
        thermal velocity is assumed: :math:`μ V^2 \sim 2 k_B T` where
        :math:`μ` is the reduced mass.

    method : `str`, optional
        The method by which to compute the Coulomb logarithm.  The
        default method is the classical straight-line Landau-Spitzer
        method (``"classical"`` or ``"ls"``). The other 6 supported
        methods are ``"ls_min_interp"``, ``"ls_full_interp"``,
        ``"ls_clamp_mininterp"``, ``"hls_min_interp"``,
        ``"hls_max_interp"``, and ``"hls_full_interp"``.  Please refer
        to the docstring of
        `~plasmapy.formulary.collisions.coulomb.Coulomb_logarithm` for
        more information about these methods.

    Returns
    -------
    mobility_value : `float` or `numpy.ndarray`
        The electrical mobility of particles in a collisional plasma.

    Raises
    ------
    `ValueError`
        If the mass or charge of either particle cannot be found, or any
        of the inputs contain incorrect values.

    `~astropy.units.UnitConversionError`
        If the units on any of the inputs are incorrect.

    `TypeError`
        If any of ``n_e``, ``T``, or ``V`` is not a
        `~astropy.units.Quantity`.

    `~plasmapy.utils.exceptions.RelativityError`
        If the input velocity is same or greater than the speed of
        light.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    : `~plasmapy.utils.exceptions.RelativityWarning`
        If the input velocity is greater than 5% of the speed of light.

    Notes
    -----
    The `mobility
    <https://en.wikipedia.org/wiki/Electrical_mobility#Mobility_in_gas_phase>`_
    is given by

    .. math::

        μ = \frac{q}{m ν}

    where :math:`q` is the particle charge, :math:`m` is the particle
    mass and :math:`ν` is the collisional frequency of the particle in
    the plasma.

    The mobility describes the forced diffusion of a particle in a
    collisional plasma which is under the influence of an electric
    field. The mobility is essentially the ratio of drift velocity due
    to collisions and the electric field driving the forced diffusion.

    Examples
    --------
    >>> import astropy.units as u
    >>> n = 1e19 * u.m ** -3
    >>> T = 1e6 * u.K
    >>> species = ('e', 'p')
    >>> mobility(T, n, species)
    <Quantity 250505... m2 / (s V)>
    >>> mobility(T, n, species, V=1e6 * u.m / u.s)
    <Quantity 1921.2784... m2 / (s V)>
    """
    freq = frequencies.collision_frequency(
        T=T, n=n_e, species=species, z_mean=z_mean, V=V, method=method
    )
    # we do this after collision_frequency since collision_frequency
    # already has a _process_inputs check and we are doing this just
    # to recover the charges, mass, etc.
    T, masses, charges, reduced_mass, V = _process_inputs(T=T, species=species, V=V)
    z_val = (charges[0] + charges[1]) / 2 if np.isnan(z_mean) else z_mean * e
    return z_val / (reduced_mass * freq)


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n={"can_be_negative": False},
)
def Spitzer_resistivity(
    T: u.K,
    n: u.m**-3,
    species,
    z_mean: Real = np.nan,
    V: u.m / u.s = np.nan * u.m / u.s,
    method="classical",
) -> u.Ohm * u.m:
    r"""
    Spitzer resistivity of a plasma.

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        Temperature in units of temperature.  This should be the
        electron temperature for electron-electron and electron-ion
        collisions, and the ion temperature for ion-ion collisions.

    n : `~astropy.units.Quantity`
        The density in units convertible to per cubic meter.  This
        should be the electron density for electron-electron collisions,
        and the ion density for electron-ion and ion-ion collisions.

    z_mean : `~astropy.units.Quantity`, optional
        The average ionization (arithmetic mean) of a plasma for which a
        macroscopic description is valid. This parameter is used to
        compute the average ion density (given the average ionization
        and electron density) for calculating the ion sphere radius for
        non-classical impact parameters. ``z_mean`` is a required
        parameter if ``method`` is ``"ls_full_interp"``,
        ``"hls_max_interp"``, or ``"hls_full_interp"``.

    species : `tuple`
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second).

    V : `~astropy.units.Quantity`, optional
        The relative velocity between particles. If not provided,
        thermal velocity is assumed: :math:`μ V^2 \sim 2 k_B T` where
        :math:`μ` is the reduced mass.

    method : `str`, optional
        The method by which to compute the Coulomb logarithm.  The
        default method is the classical straight-line Landau-Spitzer
        method (``"classical"`` or ``"ls"``). The other 6 supported
        methods are ``"ls_min_interp"``, ``"ls_full_interp"``,
        ``"ls_clamp_mininterp"``, ``"hls_min_interp"``,
        ``"hls_max_interp"``, and ``"hls_full_interp"``.  Please refer
        to the docstring of
        `~plasmapy.formulary.collisions.coulomb.Coulomb_logarithm` for
        more information about these methods.

    Returns
    -------
    spitzer : `float` or `numpy.ndarray`
        The resistivity of the plasma in ohm meters.

    Raises
    ------
    `ValueError`
        If the mass or charge of either particle cannot be found, or any
        of the inputs contain incorrect values.

    `~astropy.units.UnitConversionError`
        If the units on any of the inputs are incorrect.

    `TypeError`
        If any of ``n_e``, ``T``, or ``V`` are not of type
        `~astropy.units.Quantity`.

    `~plasmapy.utils.exceptions.RelativityError`
        If the input velocity is same or greater than the speed of
        light.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    : `~plasmapy.utils.exceptions.RelativityWarning`
        If the input velocity is greater than 5% of the speed of light.

    Notes
    -----
    The Spitzer resistivity (see Ch. 5 of :cite:t:`chen:2016`) is given
    by:

    .. math::

        η = \frac{m}{n Z_1 Z_2 q_e^2} ν_{1,2}

    where :math:`m` is the ion mass or the reduced mass, :math:`n` is
    the ion density, :math:`Z` is the particle charge state, :math:`q_e`
    is the charge of an electron, :math:`ν_{1,2}` is the collisional
    frequency between particle species 1 and 2.

    Typically, particle species 1 and 2 are selected to be an electron
    and an ion, since electron-ion collisions are inelastic and
    therefore produce resistivity in the plasma.

    Examples
    --------
    >>> import astropy.units as u
    >>> n = 1e19 * u.m ** -3
    >>> T = 1e6 * u.K
    >>> species = ('e', 'p')
    >>> Spitzer_resistivity(T, n, species)
    <Quantity 2.4915...e-06 m Ohm>
    >>> Spitzer_resistivity(T, n, species, V=1e6 * u.m / u.s)
    <Quantity 0.000324... m Ohm>
    """
    # collisional frequency
    freq = frequencies.collision_frequency(
        T=T, n=n, species=species, z_mean=z_mean, V=V, method=method
    )
    # fetching additional parameters
    T, masses, charges, reduced_mass, V = _process_inputs(T=T, species=species, V=V)
    return (
        freq * reduced_mass / (n * charges[0] * charges[1])
        if np.isnan(z_mean)
        else freq * reduced_mass / (n * (z_mean * e) ** 2)
    )
