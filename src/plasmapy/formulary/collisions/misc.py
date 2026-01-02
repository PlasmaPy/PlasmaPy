"""
Module of miscellaneous parameters related to collisions.
"""

__all__ = [
    "mobility",
    "Bethe_stopping",
    "Moliere_scattering",
    "Highland_scattering",
    "Spitzer_resistivity",
]
__lite_funcs__ = ["Bethe_stopping_lite"]

from collections.abc import Callable
from typing import Any

import astropy.constants as const
import astropy.units as u
import numpy as np
import numpy.typing as npt
from scipy.integrate import quad
from scipy.interpolate import make_interp_spline
from scipy.optimize import fsolve
from scipy.special import factorial, jv

from plasmapy.formulary.collisions import frequencies
from plasmapy.formulary.speeds import thermal_speed
from plasmapy.particles.atomic import reduced_mass, transmitted_energy_from_thickness
from plasmapy.particles.atomic import (
    stopping_power as NIST_stopping_power,  # noqa: N812
)
from plasmapy.particles.atomic import (
    stopping_range as NIST_stopping_range,  # noqa: N812
)
from plasmapy.particles.decorators import particle_input
from plasmapy.particles.particle_class import Particle
from plasmapy.particles.particle_collections import ParticleListLike
from plasmapy.utils.decorators import bind_lite_func, validate_quantities
from plasmapy.utils.decorators.checks import _check_relativistic
from plasmapy.utils.exceptions import PhysicsError

__all__ += __lite_funcs__

_a0 = const.a0
_alpha = const.alpha
_c = const.c
_e = const.e.si
_eps0 = const.eps0
_hbar = const.hbar
_m_e = const.m_e
_N_A = const.N_A

_ϑ_tabulated = [
    0.0,
    0.2,
    0.4,
    0.6,
    0.8,
    1.0,
    1.2,
    1.4,
    1.6,
    1.8,
    2.0,
    2.2,
    2.4,
    2.6,
    2.8,
    3.0,
    3.2,
    3.4,
    3.6,
    3.8,
    4.0,
    4.5,
    5.0,
    5.5,
    6.0,
    7.0,
    8.0,
    9.0,
    10.0,
]

# TODO: Validate the values here by comparing with actual integration
_f_mol_tabulated = [
    [
        +2.0000,
        +0.8456,
        +2.4929,
    ],
    [
        +1.9216,
        +0.7038,
        +2.0694,
    ],
    [
        +1.7043,  # Corrected, old value +1.7214
        +0.3437,
        +1.0488,
    ],
    [
        +1.3954,  # Corrected, old value +1.4094
        -0.0777,
        -0.0044,
    ],
    [
        +1.0546,
        -0.3981,
        -0.6068,
    ],
    [
        +0.7338,
        -0.5285,
        -0.6359,
    ],
    [
        +0.4738,
        -0.4770,
        -0.3086,
    ],
    [
        +0.2817,
        -0.3183,
        +0.0525,
    ],
    [
        +0.1546,
        -0.1396,
        +0.2423,
    ],
    [
        +0.0783,
        -6.1e-4,  # Corrected, old value -6e-4
        +0.2386,
    ],
    [
        +0.03660,
        +0.07820,
        +0.1316,
    ],
    [
        +0.01581,
        +0.10540,
        +0.0196,
    ],
    [
        +0.00630,
        +0.10080,
        -0.0467,
    ],
    [
        +0.00232,
        +0.08262,
        -0.0649,
    ],
    [
        +0.00079,
        +0.06247,
        -0.0546,
    ],
    [+0.000247, +0.04550, -0.03568],  # Corrected f0, old value +0.000250
    [+7.1e-5, +0.03288, -0.01923],  # Corrected f0, old value +7.3e-5
    [+1.9e-5, +0.02402, -0.00847],
    [+4.7e-6, +0.01791, -0.00264],
    [+1.07e-6, +0.01366, +4.6e-5],  # Corrected f0 and f2, old values +1.1e-6 and +5e-5
    # 4
    [1e-3 * 2.3e-4, 1e-3 * 10.638, 1e-3 * 1.0741],
    [1e-3 * 3e-6, 1e-3 * 6.140, 1e-3 * 1.2294],
    # 5
    [1e-3 * 2e-8, 1e-3 * 3.831, 1e-3 * 0.8326],
    [1e-3 * 2e-10, 1e-3 * 2.527, 1e-3 * 0.5368],
    [1e-3 * 5e-13, 1e-3 * 1.739, 1e-3 * 0.3495],
    [1e-3 * 1e-18, 1e-3 * 0.9080, 1e-3 * 0.1584],
    [1e-3 * 3e-25, 1e-3 * 0.5211, 1e-3 * 0.0783],
    [1e-3 * 1e-32, 1e-3 * 0.3208, 1e-3 * 0.0417],
    [1e-3 * 1e-40, 1e-3 * 0.2084, 1e-3 * 0.0237],
]

_f_mol_spline = make_interp_spline(_ϑ_tabulated, _f_mol_tabulated)


@validate_quantities(T={"equivalencies": u.temperature_energy()})
@particle_input
def _process_inputs(T: u.Quantity[u.K], species: (Particle, Particle), V):
    """
    Helper function for processing inputs to functionality contained
    in `plasmapy.formulary.collisions`.

    Also obtains the reduced mass in a 2 particle collision system
    along with thermal velocity.
    """
    masses = [p.mass for p in species]
    charges = [np.abs(p.charge) for p in species]

    # obtaining reduced mass of 2 particle collision system
    reduced_mass_ = reduced_mass(*species)

    # getting thermal velocity of system if no velocity is given
    V = _replace_nan_velocity_with_thermal_velocity(V, T, reduced_mass_)

    _check_relativistic(V, "V")

    return T, masses, charges, reduced_mass_, V


# TODO: Remove redundant mass parameter
def _replace_nan_velocity_with_thermal_velocity(
    V,
    T,
    m,
    species=Particle("e-"),  # noqa: B008
):
    """
    Get thermal velocity of system if no velocity is given, for a given
    mass.  Handles vector checks for ``V``, you must already know that
    ``T`` and ``m`` are okay.
    """
    if np.any(V == 0):
        raise PhysicsError("Collisions are not possible with a zero velocity.")

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
    T: u.Quantity[u.K],
    n_e: u.Quantity[u.m**-3],
    species,
    z_mean: float = np.nan,
    V: u.Quantity[u.m / u.s] = np.nan * u.m / u.s,
    method: str = "classical",
) -> u.Quantity[u.m**2 / (u.V * u.s)]:
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
    .. autolink-skip:: section

    >>> import astropy.units as u
    >>> n = 1e19 * u.m**-3
    >>> T = 1e6 * u.K
    >>> species = ("e", "p")
    >>> mobility(T, n, species)  # doctest: +SKIP
    <Quantity 250505... m2 / (V s)>
    >>> mobility(T, n, species, V=1e6 * u.m / u.s)  # doctest: +SKIP
    <Quantity 1921.2784... m2 / (V s)>
    """
    freq = frequencies.collision_frequency(
        T=T, n=n_e, species=species, z_mean=z_mean, V=V, method=method
    )
    # we do this after collision_frequency since collision_frequency
    # already has a _process_inputs check and we are doing this just
    # to recover the charges, mass, etc.
    T, _masses, charges, reduced_mass_, V = _process_inputs(T=T, species=species, V=V)
    z_val = (charges[0] + charges[1]) / 2 if np.isnan(z_mean) else z_mean * _e
    return z_val / (reduced_mass_ * freq)


def Bethe_stopping_lite(
    I: npt.NDArray[np.integer[Any] | np.floating[Any]],  # noqa: E741
    n: npt.NDArray[np.integer[Any] | np.floating[Any]],
    v: npt.NDArray[np.integer[Any] | np.floating[Any]],
    z: int,
) -> npt.NDArray[np.integer[Any] | np.floating[Any]]:
    r"""
    The :term:`lite-function` version of `~plasmapy.formulary.collisions.misc.Bethe_stopping`. Performs the same
    calculations as `~plasmapy.formulary.collisions.misc.Bethe_stopping`, but is intended for computational use
    and thus has data conditioning safeguards removed.

    The theoretical electronic stopping power for swift charged particles
    calculated from the Bethe formula.

    The Bethe formula should only be used for high energy particles, as
    higher order corrections become non-negligible for smaller energies.

    By convention, this function returns a positive value for the stopping
    power.

    Parameters
    ----------
    I: `~numpy.ndarray`
        The mean excitation energy for the material in which the particle is
        being stopped. Expressed in units of energy.

    n: `~numpy.ndarray`
        The electron number density of the material. Expressed in units of number density.

    v: `~numpy.ndarray`
        The velocity of the particle being stopped. Expressed in units of speed.

    z: `int`
        The charge of the charged particle in multiples of the electron charge.
        Expressed only as an integer.

    Returns
    -------
    dEdx : `~numpy.ndarray`
        The stopping power of the material given the particle's energy.

    """

    beta = v / _c.si.value

    return -np.asarray(
        4
        * np.pi
        * n
        * z**2
        / (_m_e.si.value * _c.si.value**2 * beta**2)
        * (_e.si.value**2 / (4 * np.pi * _eps0.si.value)) ** 2
        * (
            np.log(2 * _m_e.si.value * _c.si.value**2 * beta**2 / (I * (1 - beta**2)))
            - beta**2
        )
    )


@bind_lite_func(Bethe_stopping_lite)
@validate_quantities()
def Bethe_stopping(
    I: u.Quantity[u.J],  # noqa: E741
    n: u.Quantity[1 / u.m**3],
    v: u.Quantity[u.m / u.s],
    z: int,
) -> u.Quantity[u.J / u.m]:
    r"""
    The theoretical electronic stopping power for swift charged particles
    calculated from the Bethe formula.

    The Bethe formula should only be used for high energy particles, as
    higher order corrections become non-negligible for smaller energies.

    By convention, this function returns a positive value for the stopping
    power.

    Parameters
    ----------
    I: `~astropy.units.Quantity`
        The mean excitation energy for the material in which the particle is
        being stopped. Expressed in units of energy.

    n: `~astropy.units.Quantity`
        The electron number density of the material. Expressed in units of number density.

    v: `~astropy.units.Quantity`
        The velocity of the particle being stopped. Expressed in units of speed.

    z: `int`
        The charge of the charged particle in multiples of the electron charge.
        Expressed only as an integer.

    Returns
    -------
    dEdx : `~astropy.units.Quantity`
        The stopping power of the material given the particle's energy.

    """

    return Bethe_stopping_lite(I.si.value, n.si.value, v.si.value, z) * u.J / u.m


def Highland_scattering(
    m: u.Quantity[u.kg],
    v: u.Quantity[u.m / u.s],
    L: u.Quantity[u.m],
    L_rad: u.Quantity[u.kg / u.m**2],
):
    r"""Calculate the rms scattering angle using the Highland formula.

    Parameters
    ----------
    m : `~astropy.Units.Quantity`
        The mass of the projectile particles streaming through the target.
    v : `~astropy.units.Quantity`
        The speed of the projectile particles streaming through the target.
    L : `~astropy.units.Quantity`
        The thickness of the target in units of length.
    L_rad : `~astropy.units.Quantity`
        The radiation length of the target material in units of length.
        See notes for more details.

    Returns
    -------
    theta_rms : `~astropy.units.Quantity`
        The rms scattering angle in radians.

    Notes
    -----
    The root-mean-square (rms) scattering angle is given in :cite:t:`highland:1975` as:

    .. math::
        \theta_{1/e} = \frac{17.5 \; \text{MeV}}{p\beta c}\sqrt{\frac{L}{L_R}}
        \left(1 + 0.125\log_{10}\left(\frac{L}{0.1L_R}\right)\right)

    where :math:`p` is the momentum of the projectile particles,
    :math:`\beta` is the relativistic beta, :math:`L` is the thickness
    of the target, and :math:`L_R` is the radiation length--a characteristic
    distance scale over which energy loss to bremsstrahlung radiation is
    relevant.

    The Highland formula is an approximation that works best for high Z targets.
    For low Z targets, the number of scattering events may be underestimated
    by the Highland formula, and a different model should be used.
    """
    # Fitting constant, value provided by Highland
    E_s = 17.5 * u.MeV

    # Eq (2) in Highland
    epsilon = 0.125 * np.log10(10 * L / L_rad)

    # Eq (4) in Highland
    return E_s / (m * v**2) * np.sqrt(L / L_rad) * (1 + epsilon)


# Moliere scattering utility functions
def _preprocess_stopping_arguments(
    beam,
    target,
    Rho,
    NIST_material,
    stopping_power,
):
    """
    Validate the arguments passed into `Moliere_scattering`. If the user has
    provided sufficient arguments, return the stopping interpolator.

    The stopping interpolator can be either user-provided or loaded from PSTAR.
    """
    if stopping_power:
        # Convert output of provided interpolator to areal length units
        _stopping_power = stopping_power

        def stopping_power(x):
            return _stopping_power(x) / Rho

        return stopping_power

    # Interpolate PSTAR (or equivalent) NIST data if no stopping power curve
    # is provided
    if not NIST_material:
        # If we're working with an elemental target, take the name of the element
        # in all uppercase
        if len(target) == 1:
            NIST_material = target[0].element_name.upper()
        else:
            raise ValueError("Please provide one of ")

    # TODO: Add error handling to handle the case we can't find the material
    #  the user is looking for.
    return NIST_stopping_power(beam, NIST_material, return_interpolator=True)


def _beta_and_p(beam: Particle, T: u.Quantity[u.MeV]):
    beta = np.sqrt(1 - 1 / (T / (beam.mass * const.c**2) + 1) ** 2)
    p = beam.mass * beta * const.c / np.sqrt(1 - beta**2)

    return beta, p


def _x_a_squared(beam, T, c_1, c_2):
    """Eq. 6."""
    beta, p = _beta_and_p(beam, T)

    # Eq. 7
    x_0_squared = c_2 / (p * const.c) ** 2

    # Eq. 8
    alpha_squared = c_1 / beta**2

    return x_0_squared * (1.13 + 3.76 * alpha_squared)


def _x_a_bar_integrand(beam, T: u.MeV, target, stopping_power, c_1, c_2):
    """
    Eq 20, for a single value of `i`.
    That is, the provided `target` must be of a pure element.
    """
    beta, p = _beta_and_p(beam, T)
    x_a_squared = _x_a_squared(beam, T, c_1, c_2)
    # Eq. 22, Fano's correction
    # TODO: Do we need this if it's for "scattering by atomic electrons"?
    D = (
        np.log(1130 / (target.atomic_number ** (4 / 3) * (1 / beta**2 - 1)))
        + 5
        - 1 / 2 * beta**2
    )

    return (np.log(x_a_squared) - D / target.atomic_number) / (
        -stopping_power(T) * (p * beta * const.c) ** 2
    )


def _calculate_characteristic_angles(
    beam,
    incident_energy,
    target,
    Rho,
    t,
    NIST_material,
    stopping_power,
    use_constant_energy_approximation,
    c_1,
    c_2,
    c_3,
):
    # TODO: Update this function to work with user-provided range data
    if use_constant_energy_approximation:
        beta, p = _beta_and_p(beam, incident_energy)
        # Eq. 3
        x_c_squared = c_3 * Rho * t / (p * beta * const.c) ** 2
        # Eq. 6
        x_a_squared = _x_a_squared(beam, incident_energy, c_1, c_2)

        return x_c_squared, x_a_squared

    stopping_power = _preprocess_stopping_arguments(
        beam, target, Rho, NIST_material, stopping_power
    )

    # The bounds of the integrals are functions of the incident and
    # transmitted kinetic energies
    _energy, beam_areal_range = NIST_stopping_range(
        beam, NIST_material, incident_energy
    )
    _thickness, transmitted_energy = transmitted_energy_from_thickness(
        beam, NIST_material, thickness=(beam_areal_range - t * Rho)
    )

    def x_c_integrand(T: u.Quantity[u.MeV]):
        """Eq. 18."""
        beta, p = _beta_and_p(beam, T)

        return 1 / (-stopping_power(T) * (p * beta * const.c) ** 2)

    c_3_summed = np.sum(c_3, axis=0)
    x_c_squared = (
        c_3_summed
        * quad(
            lambda x: x_c_integrand(x * u.MeV).cgs.value,
            incident_energy.to(u.MeV).value,
            transmitted_energy.to(u.MeV).value,
        )[0]
    )
    # Apply units associated with the integration
    x_c_squared *= x_c_integrand(1 * u.MeV).cgs.unit * u.MeV
    ln_x_a_bar_squared = np.sum(
        u.Quantity(
            [
                c_3_i
                / x_c_squared
                * quad(
                    lambda T,
                    target_i=target_i,
                    c_1_i=c_1_i,
                    c_2_i=c_2_i: _x_a_bar_integrand(
                        beam, T * u.MeV, target_i, stopping_power, c_1_i, c_2_i
                    ).cgs.value,
                    incident_energy.to(u.MeV).value,
                    transmitted_energy.to(u.MeV).value,
                )[0]
                for (target_i, c_1_i, c_2_i, c_3_i) in zip(
                    target, c_1, c_2, c_3, strict=True
                )
            ]
        )
    )

    cgs_integrand_units = _x_a_bar_integrand(
        beam, 1 * u.MeV, target[0], stopping_power, c_1[0], c_2[0]
    ).cgs.unit
    ln_x_a_bar_squared *= cgs_integrand_units * u.MeV
    x_a_squared = np.exp(ln_x_a_bar_squared.cgs)

    return x_c_squared, x_a_squared


def _f_n_mol_integrand(
    u: float,
    ϑ: float,
    n: int,
):
    """Eq. 26 of Bethe."""
    # TODO: Move this into the body of `Bethe_Ferrari_Moliere_scattering`
    return u * jv(0, ϑ * u) * np.exp(-(u**2) / 4) * (u**2 / 4 * np.log(u**2 / 4)) ** n


def _f_mol_n(
    ϑ: u.Quantity[u.dimensionless_unscaled],
    n: int,
):
    integral = quad(_f_n_mol_integrand, 0, np.inf, args=(ϑ, n))[0]

    return integral / factorial(n)


def Moliere_scattering_B_residual(B, b):
    """Eq. 23 of Bethe."""
    return B - np.log(B) - b


def _wrapped_Moliere_angular_distribution(x_c_squared, B, use_f_mol_interpolator):
    def Moliere_angular_distribution(theta):
        theta_prime = (theta / np.sqrt(x_c_squared * B)).cgs
        B_coefficients = np.asarray([1 / B**i for i in range(3)]).T

        if use_f_mol_interpolator:
            f_n = _f_mol_spline(theta_prime)
        else:
            f_n = np.asarray([_f_mol_n(theta_prime, i) for i in range(3)])

        result = np.sum(B_coefficients * f_n, axis=-1)
        return theta_prime * result

    return Moliere_angular_distribution


@particle_input
def Moliere_scattering(
    # Projectile parameters
    beam: Particle,
    incident_energy: u.Quantity[u.MeV],
    # Target parameters
    target: ParticleListLike,
    Rho: u.Quantity[u.g / u.cm**3],
    t: u.Quantity[u.m],
    f_i: npt.NDArray[np.float64] | None = None,
    # Miscellaneous parameters
    use_constant_energy_approximation: bool = False,
    NIST_material: str | None = None,
    stopping_power: Callable[[u.Quantity[u.MeV]], u.Quantity[u.MeV / u.g]]
    | None = None,
    use_f_mol_interpolator: bool = True,
    return_rms: bool = False,
):
    """Calculate the angular scattering distribution for the provided particle species.

    Parameters
    ----------
    beam : `~plasmapy.particles.particle_class.Particle`
        The particle specie that is incident upon the target.

    incident_energy : `~astropy.units.Quantity`
        The kinetic energy of the beam particles in units convertible to MeV.

    target : `~plasmapy.particles.particle_collections.ParticleListLike`
        A |Particle| representation of the element(s) in the target.

    Rho : `~astropy.units.Quantity`
        The mass density of the target element(s) in units convertible to grams
        per cubic centimeter.

    t : `~astropy.units.Quantity`
        The thickness of the target in units convertible to centimeters.

    f_i : array of `float`
        The fractional weight of each constituent element's atoms in the target.

    use_constant_energy_approximation : bool, optional
        Whether to calculate the angular probability distribution using the
        approximation that the beam particles do not lose energy to the
        material. Defaults to `False`.

    NIST_material : str, optional
        The name of the material in the NIST PSTAR database. Providing this
        argument is only necessary if `stopping_power` is `None`.

    stopping_power : `Callable[[u.Quantity[u.MeV]], u.MeV / u.g]`
        Function used to calculate the stopping power. Takes in kinetic energy
        and returns the stopping power in units of energy per
        unit length.

    use_f_mol_interpolator : `bool`
        Whether the integrals involved in the calculation of the angular
        distribution should be estimated using interpolation or numerical
        integration using `~scipy.integration.quad`. Defaults to `True`.

    return_rms : `bool`
        Whether to return the rms scattering angle in radians, or return the
        distribution function. Defaults to `False`.

    """
    # In the case that we are dealing with a compound target, iterate over the constituent atoms
    c_1 = u.Quantity(
        [
            (
                const.e.esu**2
                / (const.hbar * const.c)
                * beam.atomic_number
                * element.atomic_number
            )
            ** 2
            for element in target
        ]
    )
    c_2 = u.Quantity(
        [
            (
                1
                / 0.885
                * (const.e.esu**2 / (const.hbar * const.c))
                * (const.m_e * const.c**2)
                * element.atomic_number ** (1 / 3)
            )
            ** 2
            for element in target
        ]
    )
    c_3 = u.Quantity(
        [
            (
                4
                * np.pi
                * const.N_A
                * const.e.esu**4
                * beam.atomic_number**2
                * element.atomic_number**2
                * f
                / (element.mass * const.N_A)
            ).cgs
            for (element, f) in zip(target, f_i, strict=True)
        ]
    )
    x_c_squared, x_a_squared = _calculate_characteristic_angles(
        beam,
        incident_energy,
        target,
        Rho,
        t,
        NIST_material,
        stopping_power,
        use_constant_energy_approximation,
        c_1,
        c_2,
        c_3,
    )

    b = np.log(x_c_squared / (1.167 * x_a_squared))

    # The transcendental equation associated with `B` yields two solutions for
    # every `b`. We want to solve for values where B > 1, this corresponds to
    # our initial guess satisfying b > 1.
    B = fsolve(Moliere_scattering_B_residual, x0=b, args=(b,))

    angular_distribution = _wrapped_Moliere_angular_distribution(
        x_c_squared, B, use_f_mol_interpolator
    )

    if not return_rms:
        return angular_distribution
    else:
        # Algebraic manipulation of Eq. 30
        theta_prime_unit = np.sqrt(x_c_squared * B)

        # Calculate the value of theta associated with theta_prime = 5
        # At this value, the value of the distributions have dropped to nearly 0
        upper_bound = 5 * theta_prime_unit
        total_scattering_probability = quad(angular_distribution, 0, upper_bound)[0]

        mean_squared = quad(
            lambda theta: theta**2
            * angular_distribution(theta)
            / total_scattering_probability,
            0,
            upper_bound,
        )[0]

        return np.sqrt(mean_squared)


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n={"can_be_negative": False},
)
def Spitzer_resistivity(
    T: u.Quantity[u.K],
    n: u.Quantity[u.m**-3],
    species,
    z_mean: float = np.nan,
    V: u.Quantity[u.m / u.s] = np.nan * u.m / u.s,
    method: str = "classical",
) -> u.Quantity[u.Ohm * u.m]:
    r"""
    Spitzer resistivity of a plasma.

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        Temperature in units of Kelvin or eV.  This should be the
        electron temperature for electron-electron and electron-ion
        collisions, and the ion temperature for ion-ion collisions. An
        example of temperature given in eV can be found below.

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
    .. autolink-skip:: section

    >>> import astropy.units as u
    >>> n = 1e19 * u.m**-3
    >>> T = 1e6 * u.K
    >>> species = ("e", "p")
    >>> Spitzer_resistivity(T, n, species)  # doctest: +SKIP
    <Quantity 2.4915...e-06 Ohm m>
    >>> Spitzer_resistivity(T, n, species, V=1e6 * u.m / u.s)  # doctest: +SKIP
    <Quantity 0.000324... Ohm m>
    >>> T_eV = 86.173 * u.eV
    >>> T_K = (T_eV).to("K", equivalencies=u.temperature_energy())
    >>> Spitzer_resistivity(T_K, n, species)
    <Quantity 2.49158...e-06 Ohm m>
    """
    # collisional frequency
    freq = frequencies.collision_frequency(
        T=T, n=n, species=species, z_mean=z_mean, V=V, method=method
    )
    # fetching additional parameters
    T, _masses, charges, reduced_mass_, V = _process_inputs(T=T, species=species, V=V)
    return (
        freq * reduced_mass_ / (n * charges[0] * charges[1])
        if np.isnan(z_mean)
        else freq * reduced_mass_ / (n * (z_mean * _e) ** 2)
    )
