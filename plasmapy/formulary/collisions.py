"""Functions to calculate transport coefficients.

This module includes a number of functions for handling Coulomb
collisions spanning weakly coupled (low density) to strongly coupled
(high density) regimes.

Coulomb collisions
==================

Coulomb collisions are collisions where the interaction force is
conveyed via the electric field, instead of any kind of contact
force. They usually result in relatively small deflections of particle
trajectories. However, given that there are many charged particles in a
plasma, one has to take into account the cumulative effects of many such
collisions.

Coulomb logarithms
==================

Please read the documentation of
`~plasmapy.formulary.collisions.Coulomb_logarithm` below for an
explanation of the seven PlasmaPy-supported methods of computing the
Coulomb logarithm.

Collision rates
===============

The module gathers a few functions helpful for calculating collision
rates between particles. The most general of these is
`~plasmapy.formulary.collisions.collision_frequency`.

Macroscopic properties
======================

These include:

* `~plasmapy.formulary.collisions.Spitzer_resistivity`
* `~plasmapy.formulary.collisions.mobility`
* `~plasmapy.formulary.collisions.Knudsen_number`
* `~plasmapy.formulary.collisions.coupling_parameter`
"""
__all__ = [
    "Coulomb_logarithm",
    "impact_parameter_perp",
    "impact_parameter",
    "collision_frequency",
    "Coulomb_cross_section",
    "fundamental_electron_collision_freq",
    "fundamental_ion_collision_freq",
    "mean_free_path",
    "Spitzer_resistivity",
    "mobility",
    "Knudsen_number",
    "coupling_parameter",
]

import astropy.units as u
import numpy as np
import warnings

from astropy.constants.si import e, eps0, hbar, k_B, m_e
from numbers import Real
from numpy import pi

from plasmapy import particles, utils
from plasmapy.formulary.lengths import Debye_length
from plasmapy.formulary.mathematics import Fermi_integral
from plasmapy.formulary.quantum import (
    chemical_potential,
    thermal_deBroglie_wavelength,
    Wigner_Seitz_radius,
)
from plasmapy.formulary.speeds import thermal_speed
from plasmapy.utils.decorators import validate_quantities
from plasmapy.utils.decorators.checks import _check_relativistic


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    V={"none_shall_pass": True},
)
@particles.particle_input
def Coulomb_logarithm(
    T: u.K,
    n_e: u.m**-3,
    species: (particles.Particle, particles.Particle),
    z_mean: Real = np.nan,
    V: u.m / u.s = np.nan * u.m / u.s,
    method="classical",
):
    r"""
    Compute the Coulomb logarithm.

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
        non-classical impact parameters. ``z_mean`` is a required
        parameter if ``method`` is ``"ls_full_interp"``,
        ``"hls_max_interp"``, or ``"hls_full_interp"``.

    V : `~astropy.units.Quantity`, optional
        The relative velocity between particles. If not provided,
        thermal velocity is assumed: :math:`μ V^2 \sim 2 k_B T` where
        :math:`μ` is the reduced mass.

    method : `str`, optional
        The method by which to compute the Coulomb logarithm.  The
        default method is the classical straight-line Landau-Spitzer
        method (``"classical"`` or ``"ls"``). The other six supported
        methods are ``"ls_min_interp"``, ``"ls_full_interp"``,
        ``"ls_clamp_mininterp"``, ``"hls_min_interp"``,
        ``"hls_max_interp"``, and ``"hls_full_interp"``.  Please refer
        to the Notes section of this docstring for more information,
        including about abbreviated aliases of these names.

    Returns
    -------
    ln_Lambda : `float` or `numpy.ndarray`
        The dimensionless Coulomb logarithm.

    Raises
    ------
    `ValueError`
        If the mass or charge of either particle cannot be found, or any
        of the inputs contain incorrect values.

    `~astropy.units.UnitConversionError`
        If the units on any of the inputs are incorrect, or if any of
        ``n_e``, ``T``, or ``V`` is not a `~astropy.units.Quantity`.

    `~plasmapy.utils.exceptions.PhysicsError`
        If the result is smaller than 1.

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
    **Summary of Supported Methods of Computing the Coulomb Logarithm**

    PlasmaPy supports seven methods of computing the Coulomb logarithm:

    1. ``"classical"`` or ``"ls"``
    2. ``"ls_min_interp"`` or ``"GMS-1"``
    3. ``"ls_full_interp"`` or ``"GMS-2"``
    4. ``"ls_clamp_mininterp"`` or ``"GMS-3"``
    5. ``"hls_min_interp"`` or ``"GMS-4"``
    6. ``"hls_max_interp"`` or ``"GMS-5"``
    7. ``"hls_full_interp"`` or ``"GMS-6"``

    Options 1–4 are straight-line Landau-Spitzer (``"ls..."``) methods
    in which the trajectory of a Coulomb collision is modeled as a
    straight line. For the straight-line Landau-Spitzer methods, the
    Coulomb logarithm (:math:`\ln{Λ}`) is defined to be:

    .. math::

        \ln{Λ} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    Options 5–7 are hyperbolic Landau-Spitzer (``"hls..."``) methods in
    which the trajectory of a Coulomb collision is modeled as a
    hyperbola. For the hyperbolic Landau-Spitzer methods, the Coulomb
    logarithm (:math:`\ln{Λ}`) is defined to be:

    .. math::

        \ln{Λ} \equiv \frac{1}{2} \ln\left(1 + \frac{b_{max}^2}{b_{min}^2} \right)

    For all 7 methods, :math:`b_{min}` and :math:`b_{max}` are the inner
    impact parameter and the outer impact parameter, respectively, for
    Coulomb collisions :cite:p:`spitzer:1962`\ ; :math:`b_{min}` and
    :math:`b_{max}` are each computed by
    `~plasmapy.formulary.collisions.impact_parameter`.

    The abbreviations of Options 2–7 (``"GMS-..."``) refer to the first
    initials of the three authors of :cite:t:`gericke:2002`.

    .. note::

        For strongly-coupled plasma, PlasmaPy recommends Option 7,
        ``"hls_full_interp"`` or ``"GMS-6"``, because of its high
        accuracy regardless of a plasma's strength of coupling.

    **Explanation of Supported Methods of Computing the Coulomb Logarithm**

    In this section, further information about each method, such as
    about interpolation and other special features, is documented.
    Please refer to :cite:t:`spitzer:1962` and :cite:t:`gericke:2002`
    for additional information about these methods.

    Option 1: ``"classical"`` or ``"ls"`` (Landau-Spitzer)

        The classical straight-line Landau-Spitzer method in which
        :math:`b_{min}` is defined to be the higher of the de Broglie
        wavelength (:math:`λ_{de Broglie}`) and the distance of closest
        approach (:math:`ρ_⟂`) if they are not equal (and either of the
        two if they are equal) and :math:`b_{max}` is defined to be the
        Debye length (:math:`λ_{Debye}`).

        .. math::

            \ln{Λ} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

        .. math::

            b_{min} \equiv
            \left\{
                \begin{array}{ll}
                           λ_{de Broglie} & \mbox{if } λ_{de Broglie} \geq ρ_⟂ \\
                           ρ_⟂         & \mbox{if } ρ_⟂ \geq λ_{de Broglie}
                \end{array}
            \right.

        .. math::

            b_{max} \equiv λ_{Debye}

        The inner impact parameter (:math:`b_{min}`) is the higher of
        :math:`λ_{de Broglie}` and :math:`ρ_⟂` because for impact
        parameters lower than :math:`λ_{de Broglie}`, quantum effects
        cause the collision to be non-Coulombic
        :cite:p:`chen:2016,fundamenski:2007`.

        The outer impact parameter (:math:`b_{max}`) is defined to be
        the Debye length (:math:`λ_{Debye}`) because at distances higher
        than the Debye length, the electric fields created by other
        particles are screened out by the electrons rearranging
        themselves.

        The uncertainty of the classical straight-line Landau-Spitzer
        method is on the order of its reciprocal.

        This method is invalid if :math:`\ln{Λ} < 2` because of the
        uncertainty of this method and is invalid if :math:`\ln{Λ} < 0`,
        which may be true if the coupling parameter is high (such as for
        nonideal, dense, cold plasmas).

        Please refer to :cite:t:`spitzer:1962` for additional
        information about this method.

    Option 2: ``"ls_min_interp"`` or ``"GMS-1"`` (Landau-Spitzer,
    interpolation of :math:`b_{min}`)

        A straight-line Landau-Spitzer method in which :math:`b_{min}`
        is interpolated between the de Broglie wavelength (:math:`λ_{de
        Broglie}`) and the distance of closest approach (:math:`ρ_⟂`)
        and :math:`b_{max}` is defined to be the Debye length
        (:math:`λ_{Debye}`).

        .. math::

            \ln{Λ} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

        .. math::

            b_{min} \equiv \sqrt{λ_{de Broglie}^2 + ρ_⟂^2}

        .. math::

            b_{max} \equiv λ_{Debye}

        This method is invalid if :math:`\ln{Λ} < 0`, which may be true
        if the coupling parameter is high (such as for nonideal, dense,
        cold plasmas).

        Note: This is the first method in Table 1 of
        :cite:t:`gericke:2002`.

    Option 3: ``"ls_full_interp"`` or ``"GMS-2"`` (Landau-Spitzer,
    interpolation of :math:`b_{min}` and :math:`b_{max}`)

        A straight-line Landau-Spitzer method in which :math:`b_{min}`
        and :math:`b_{max}` are each interpolated. :math:`b_{min}` is
        interpolated between the de Broglie wavelength (:math:`λ_{de
        Broglie}`) and the distance of closest approach (:math:`ρ_⟂`).
        :math:`b_{max}` is interpolated between the Debye length
        (:math:`λ_{Debye}`) and the ion sphere radius (:math:`a_i`).

        .. math::

            \ln{Λ} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

        .. math::

            b_{min} \equiv \sqrt{λ_{de Broglie}^2 + ρ_⟂^2}

        .. math::

            b_{max} \equiv \sqrt{λ_{Debye}^2 + a_i^2}

        This method is invalid if :math:`\ln{Λ} < 0`, which may be true
        if the coupling parameter is high (such as for nonideal, dense,
        cold plasmas).

        Note: This is the second method in Table 1 of
        :cite:t:`gericke:2002`.

    Option 4: ``"ls_clamp_mininterp"`` or ``"GMS-3"`` (Landau-Spitzer
    with a clamp, interpolation of :math:`b_{min}`)

        A straight-line Landau-Spitzer method in which the value of
        :math:`\ln{Λ}` is clamped at a minimum of :math:`2` so that it
        is impossible for :math:`\ln{Λ} < 0` (which is possible by the
        classical Landau-Spitzer method). :math:`b_{min}` is
        interpolated between the de Broglie wavelength (:math:`λ_{de
        Broglie}`) and the distance of closest approach
        (:math:`ρ_⟂`). :math:`b_{max}` is defined to be the Debye length
        (:math:`λ_{Debye}`).

        .. math::

            \ln{Λ} \equiv
            \left\{
               \begin{array}{ll}
                  \ln\left(\frac{b_{max}}{b_{min}} \right)
                     & \mbox{if } \ln\left( \frac{b_{max}}{b_{min}}\right)
                  \geq 2 \\
                  2 & \mbox{if } \ln\left( \frac{b_{max}}{b_{min}} \right)
                  \leq 2
                \end{array}
            \right.

        .. math::

            b_{min} \equiv \sqrt{λ_{de Broglie}^2 + ρ_⟂^2}

        .. math::

            b_{max} \equiv λ_{Debye}

        This method is valid for any plasma because it is impossible for
        :math:`\ln{Λ} < 0` by this method, even if the coupling
        parameter is high (such as for nonideal, dense, cold plasmas).

        Note: This is the third method in Table 1 of
        :cite:t:`gericke:2002`.

    Option 5: ``"hls_min_interp"`` or ``"GMS-4"`` (Hyperbolic
    Landau-Spitzer, interpolation of :math:`b_{min}`)

        A hyperbolic Landau-Spitzer method in which :math:`b_{min}` is
        interpolated between the de Broglie wavelength (:math:`λ_{de
        Broglie}`) and the distance of closest approach (:math:`ρ_⟂`)
        and :math:`b_{max}` is defined to be the Debye length
        (:math:`λ_{Debye}`).

        .. math::

            \ln{Λ} \equiv \frac{1}{2} \ln\left(1 +
            \frac{b_{max}^2}{b_{min}^2} \right)

        .. math::

            b_{min} \equiv \sqrt{λ_{de Broglie}^2 + ρ_⟂^2}

        .. math::

            b_{max} \equiv λ_{Debye}

        This method is valid for any plasma because it is impossible for
        :math:`\ln{Λ} < 0` by this method, even if the coupling
        parameter is high (such as for nonideal, dense, cold plasmas).

        Note: This is the fourth method in Table 1 of
        :cite:t:`gericke:2002`.

    Option 6: ``"hls_max_interp"`` or ``"GMS-5"`` (Hyperbolic
    Landau-Spitzer, interpolation of :math:`b_{max}`)

        A hyperbolic Landau-Spitzer method in which :math:`b_{max}` is
        interpolated between the Debye length (:math:`λ_{Debye}`) and
        the ion sphere radius (:math:`a_i`) and :math:`b_{min}` is
        defined to be the distance of closest approach (:math:`ρ_⟂`).

        .. math::

            \ln{Λ} \equiv \frac{1}{2} \ln\left(1 +
            \frac{b_{max}^2}{b_{min}^2} \right)

        .. math::

            b_{min} \equiv ρ_⟂

        .. math::

            b_{max} \equiv \sqrt{λ_{Debye}^2 + a_i^2}

        This method is valid for any plasma because it is impossible for
        :math:`\ln{Λ} < 0` by this method, even if the coupling
        parameter is high (such as for nonideal, dense, cold plasmas).

        Caution: This method overestimates :math:`\ln{Λ}` at high
        temperatures.

        Note: This is the fifth method in Table 1 of
        :cite:t:`gericke:2002`.

    Option 7: ``"hls_full_interp"`` or ``"GMS-6"`` (Hyperbolic
    Landau-Spitzer, interpolation of :math:`b_{min}` and
    :math:`b_{max}`)

        A hyperbolic Landau-Spitzer method in which :math:`b_{min}` and
        :math:`b_{max}` are each interpolated. :math:`b_{min}` is
        interpolated between the de Broglie wavelength (:math:`λ_{de
        Broglie}`) and the distance of closest approach (:math:`ρ_⟂`).
        :math:`b_{max}` is interpolated between the Debye length
        (:math:`λ_{Debye}`) and the ion sphere radius (:math:`a_i`).

        .. math::

            \ln{Λ} \equiv \frac{1}{2} \ln\left(1 +
            \frac{b_{max}^2}{b_{min}^2} \right)

        .. math::

            b_{min} \equiv \sqrt{λ_{de Broglie}^2 + ρ_⟂^2}

        .. math::

            b_{max} \equiv \sqrt{λ_{Debye}^2 + a_i^2}

        This method is valid for any plasma because it is impossible for
        :math:`\ln{Λ} < 0` by this method, even if the coupling
        parameter is high (such as for nonideal, dense, cold plasmas).

        Note: This is the sixth method in Table 1 of
        :cite:t:`gericke:2002`.

    Examples
    --------
    >>> import astropy.units as u
    >>> n = 1e19 * u.m**-3
    >>> T = 1e6 * u.K
    >>> Coulomb_logarithm(T, n, ('e-', 'p+'))
    14.545527...
    >>> Coulomb_logarithm(T, n, ('e-', 'p+'), V = 1e6 * u.m / u.s)
    11.363478...

    See Also
    --------
    ~plasmapy.formulary.collisions.impact_parameter : Computes
        :math:`b_{min}` and :math:`b_{max}`.
    """
    # fetching impact min and max impact parameters
    bmin, bmax = impact_parameter(
        T=T, n_e=n_e, species=species, z_mean=z_mean, V=V, method=method
    )

    if method in (
        "classical",
        "ls",
        "ls_min_interp",
        "GMS-1",
        "ls_full_interp",
        "GMS-2",
    ):
        ln_Lambda = np.log(bmax / bmin)
    elif method in ("ls_clamp_mininterp", "GMS-3"):
        ln_Lambda = np.log(bmax / bmin)
        if np.any(ln_Lambda < 2):
            if np.isscalar(ln_Lambda.value):
                ln_Lambda = 2 * u.dimensionless_unscaled
            else:
                ln_Lambda[ln_Lambda < 2] = 2 * u.dimensionless_unscaled
    elif method in (
        "hls_min_interp",
        "GMS-4",
        "hls_max_interp",
        "GMS-5",
        "hls_full_interp",
        "GMS-6",
    ):
        ln_Lambda = 0.5 * np.log(1 + bmax**2 / bmin**2)
    else:
        raise ValueError(
            'Unknown method. Choose from "classical", "ls_min_interp", '
            '"ls_full_interp", "ls_clamp_mininterp", "hls_min_interp", '
            '"hls_max_interp", "hls_full_interp", and their aliases. '
            "Please refer to the documentation of this function for "
            "more information."
        )

    # applying dimensionless units
    ln_Lambda = ln_Lambda.to(u.dimensionless_unscaled).value

    # Allow NaNs through the < checks without warning
    with np.errstate(invalid="ignore"):
        if np.any(ln_Lambda < 2) and method in [
            "classical",
            "ls",
            "ls_min_interp",
            "GMS-1",
            "ls_full_interp",
            "GMS-2",
        ]:
            warnings.warn(
                f"The Coulomb logarithm is {ln_Lambda}, and the specified "
                f'method, "{method}", depends on weak coupling.',
                utils.CouplingWarning,
            )
        elif np.any(ln_Lambda < 4):
            warnings.warn(
                f"The Coulomb logarithm is {ln_Lambda}, so strong "
                "coupling effects may exist for the plasma.",
                utils.CouplingWarning,
            )

    return ln_Lambda


@validate_quantities(T={"equivalencies": u.temperature_energy()})
@particles.particle_input
def _process_inputs(T: u.K, species: (particles.Particle, particles.Particle), V):
    """
    Check the inputs to functions in ``collisions.py``.  Also obtains
    the reduced mass in a 2 particle collision system along with thermal
    velocity.
    """
    masses = [p.mass for p in species]
    charges = [np.abs(p.charge) for p in species]

    # obtaining reduced mass of 2 particle collision system
    reduced_mass = particles.reduced_mass(*species)

    # getting thermal velocity of system if no velocity is given
    V = _replace_nan_velocity_with_thermal_velocity(V, T, reduced_mass)

    _check_relativistic(V, "V")

    return T, masses, charges, reduced_mass, V


def _replace_nan_velocity_with_thermal_velocity(V, T, m):
    """
    Get thermal velocity of system if no velocity is given, for a given
    mass.  Handles vector checks for ``V``, you must already know that
    ``T`` and ``m`` are okay.
    """
    if np.any(V == 0):
        raise utils.PhysicsError("Collisions are not possible with a zero velocity.")

    if V is None:
        return thermal_speed(T, "e-", mass=m)

    if not np.any(np.isnan(V)):
        return V

    if not (np.isscalar(V.value) or np.isscalar(T.value)):
        V = V.copy()
        V[np.isnan(V)] = thermal_speed(T[np.isnan(V)], "e-", mass=m)
        return V

    if np.isscalar(V.value):
        return thermal_speed(T, "e-", mass=m)

    if np.isscalar(T.value):
        V[np.isnan(V)] = thermal_speed(T, "e-", mass=m)

    return V


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()}
)
@particles.particle_input
def impact_parameter_perp(
    T: u.K,
    species: (particles.Particle, particles.Particle),
    V: u.m / u.s = np.nan * u.m / u.s,
) -> u.m:
    r"""
    Distance of closest approach for a 90° Coulomb collision.

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and the
        target particle

    species : `tuple`
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    V : `~astropy.units.Quantity`, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`μ V^2 \sim 2 k_B T` where
        :math:`μ` is the reduced mass.

    Returns
    -------
    impact_parameter_perp : `float` or `numpy.ndarray`
        The distance of closest approach for a 90° Coulomb collision.

    Raises
    ------
    `ValueError`
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    `~astropy.units.UnitConversionError`
        If the units on any of the inputs are incorrect.

    `TypeError`
        If either of ``T`` or ``V`` is not a `~astropy.units.Quantity`.

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
    The distance of closest approach, impact_parameter_perp, is given by
    (see Ch. 5 of :cite:t:`chen:2016`):

    .. math::

        b_⟂ = \frac{Z_1 Z_2}{4 π \epsilon_0 m v^2}

    Examples
    --------
    >>> import astropy.units as u
    >>> T = 1e6 * u.K
    >>> species = ('e', 'p')
    >>> impact_parameter_perp(T, species)
    <Quantity 8.3550...e-12 m>
    """
    # Note: This formulation corresponds to collisions that result in a deflection of 90°s,
    #       which is valid when classical effects dominate.
    # TODO: need to incorporate an average ionization parameter

    T, masses, charges, reduced_mass, V = _process_inputs(T=T, species=species, V=V)

    return charges[0] * charges[1] / (4 * pi * eps0 * reduced_mass * V**2)


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
    V={"none_shall_pass": True},
)
def impact_parameter(
    T: u.K,
    n_e: u.m**-3,
    species,
    z_mean: Real = np.nan,
    V: u.m / u.s = np.nan * u.m / u.s,
    method="classical",
):
    r"""
    Impact parameters for classical and quantum Coulomb collision.

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and the
        target particle.

    n_e : `~astropy.units.Quantity`
        The electron number density in units convertible to per cubic
        meter.

    species : `tuple`
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second).

    z_mean : `~astropy.units.Quantity`, optional
        The average ionization (arithmetic mean) of a plasma for which a
        macroscopic description is valid. This parameter is used to
        compute the average ion density (given the average ionization
        and electron density) for calculating the ion sphere radius for
        non-classical impact parameters. ``z_mean`` is a required
        parameter if ``method`` is ``"ls_full_interp"``,
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
        to the docstring of `Coulomb_logarithm` for more information
        about these methods.

    Returns
    -------
    bmin, bmax : `tuple` of floats
        The minimum and maximum impact parameters (distances) for a
        Coulomb collision.

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
    The minimum and maximum impact parameters may be calculated in a
    variety of ways. The maximum impact parameter is typically the Debye
    length.

    For quantum plasmas the maximum impact parameter can be the
    quadratic sum of the debye length and ion radius (Wigner_Seitz)
    :cite:p:`gericke:2002`

    .. math::

        b_{max} = \left(λ_{De}^2 + a_i^2\right)^{1/2}

    The minimum impact parameter is typically some combination of the
    thermal de Broglie wavelength and the distance of closest approach
    for a 90° Coulomb collision. A quadratic sum is used for all GMS
    methods, except for GMS-5, where ``b_min`` is simply set to the
    distance of closest approach :cite:p:`gericke:2002`.

    .. math::

        b_{min} = \left(Λ_{de Broglie}^2 + ρ_⟂^2\right)^{1/2}

    Examples
    --------
    >>> import astropy.units as u
    >>> n = 1e19 * u.m**-3
    >>> T = 1e6 * u.K
    >>> species = ('e', 'p')
    >>> impact_parameter(T, n, species)
    (<Quantity 1.051...e-11 m>, <Quantity 2.182...e-05 m>)
    >>> impact_parameter(T, n, species, V=1e6 * u.m / u.s)
    (<Quantity 2.534...e-10 m>, <Quantity 2.182...e-05 m>)
    """
    T, masses, charges, reduced_mass, V = _process_inputs(T=T, species=species, V=V)
    # catching error where mean charge state is not given for non-classical
    # methods that require the ion density
    if method in (
        "ls_full_interp",
        "GMS-2",
        "hls_max_interp",
        "GMS-5",
        "hls_full_interp",
        "GMS-6",
    ) and np.isnan(z_mean):
        raise ValueError(
            'Must provide a z_mean for "ls_full_interp", '
            '"hls_max_interp", and "hls_full_interp" methods.'
        )
    # Debye length
    lambdaDe = Debye_length(T, n_e)
    # de Broglie wavelength
    lambdaBroglie = hbar / (2 * reduced_mass * V)
    # distance of closest approach in 90° Coulomb collision
    bPerp = impact_parameter_perp(T=T, species=species, V=V)

    # obtaining minimum and maximum impact parameters depending on which
    # method is requested
    if method in ["classical", "ls"]:
        bmax = lambdaDe
        # Coulomb-style collisions will not happen for impact parameters
        # shorter than either of these two impact parameters, so we choose
        # the larger of these two possibilities. That is, between the
        # de Broglie wavelength and the distance of closest approach.
        # ARRAY NOTES
        # T and V should be guaranteed to be same size inputs from _process_inputs
        # therefore, lambdaBroglie and bPerp are either both scalar or both array
        # if np.isscalar(bPerp.value) and np.isscalar(lambdaBroglie.value):  # both scalar
        try:  # assume both scalar
            bmin = bPerp if bPerp > lambdaBroglie else lambdaBroglie
        except ValueError:  # both lambdaBroglie and bPerp are arrays
            bmin = lambdaBroglie
            bmin[bPerp > lambdaBroglie] = bPerp[bPerp > lambdaBroglie]
    elif method in ["ls_min_interp", "GMS-1"]:
        # 1st method listed in Table 1 of reference [1]
        # This is just another form of the classical Landau-Spitzer
        # approach, but bmin is interpolated between the de Broglie
        # wavelength and distance of closest approach.
        bmax = lambdaDe
        bmin = (lambdaBroglie**2 + bPerp**2) ** (1 / 2)
    elif method in ["ls_full_interp", "GMS-2"]:
        # 2nd method listed in Table 1 of reference [1]
        # Another Landau-Spitzer like approach, but now bmax is also
        # being interpolated. The interpolation is between the Debye
        # length and the ion sphere radius, allowing for descriptions
        # of dilute plasmas.
        # Mean ion density.
        n_i = n_e / z_mean
        # mean ion sphere radius.
        ionRadius = Wigner_Seitz_radius(n_i)
        bmax = (lambdaDe**2 + ionRadius**2) ** (1 / 2)
        bmin = (lambdaBroglie**2 + bPerp**2) ** (1 / 2)
    elif method in ["ls_clamp_mininterp", "GMS-3"]:
        # 3rd method listed in Table 1 of reference [1]
        # same as GMS-1, but not Lambda has a clamp at Lambda_min = 2
        # where Lambda is the argument to the Coulomb logarithm.
        bmax = lambdaDe
        bmin = (lambdaBroglie**2 + bPerp**2) ** (1 / 2)
    elif method in ["hls_min_interp", "GMS-4"]:
        # 4th method listed in Table 1 of reference [1]
        bmax = lambdaDe
        bmin = (lambdaBroglie**2 + bPerp**2) ** (1 / 2)
    elif method in ["hls_max_interp", "GMS-5"]:
        # 5th method listed in Table 1 of reference [1]
        # Mean ion density.
        n_i = n_e / z_mean
        # mean ion sphere radius.
        ionRadius = Wigner_Seitz_radius(n_i)
        bmax = (lambdaDe**2 + ionRadius**2) ** (1 / 2)
        bmin = bPerp
    elif method in ["hls_full_interp", "GMS-6"]:
        # 6th method listed in Table 1 of reference [1]
        # Mean ion density.
        n_i = n_e / z_mean
        # mean ion sphere radius.
        ionRadius = Wigner_Seitz_radius(n_i)
        bmax = (lambdaDe**2 + ionRadius**2) ** (1 / 2)
        bmin = (lambdaBroglie**2 + bPerp**2) ** (1 / 2)
    else:
        raise ValueError(f"Method {method} not found!")

    # ARRAY NOTES
    # it could be that bmin and bmax have different sizes. If Te is a scalar,
    # T and V will be scalar from _process_inputs, so bmin will scalar. However
    # if n_e is an array, then bmax will be an array. if this is the case,
    # we want to extend the scalar bmin to match the dimensions of bmax.
    if bmin.size == 1 and bmax.size != 1:
        bmin = bmin * np.ones(bmax.shape)

    return bmin.to(u.m), bmax.to(u.m)


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n={"can_be_negative": False},
)
def collision_frequency(
    T: u.K,
    n: u.m**-3,
    species,
    z_mean: Real = np.nan,
    V: u.m / u.s = np.nan * u.m / u.s,
    method="classical",
) -> u.Hz:
    r"""
    Collision frequency of particles in a plasma.

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

    species : `tuple`
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second).

    z_mean : `~astropy.units.Quantity`, optional
        The average ionization (arithmetic mean) of a plasma for which a
        macroscopic description is valid. This parameter is used to
        compute the average ion density (given the average ionization
        and electron density) for calculating the ion sphere radius for
        non-classical impact parameters. ``z_mean`` is a required
        parameter if ``method`` is ``"ls_full_interp"``,
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
        `~plasmapy.formulary.collisions.Coulomb_logarithm` for more
        information about these methods.

    Returns
    -------
    freq : float or numpy.ndarray
        The collision frequency of particles in a plasma.

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
        If the input velocity is same or greater than the speed
        of light.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed

    : `~plasmapy.utils.exceptions.RelativityWarning`
        If the input velocity is greater than 5% of the speed of light.

    Notes
    -----
    The collision frequency (see Ch. 5 of :cite:t:`chen:2016`) is given
    by

    .. math::

        ν = n σ v \ln{Λ}

    where :math:`n` is the particle density, :math:`σ` is the
    collisional cross-section, :math:`v` is the inter-particle velocity
    (typically taken as the thermal velocity), and :math:`\ln{Λ}` is the
    Coulomb logarithm accounting for small angle collisions.

    See Equation (2.14) in :cite:t:`callen:unpublished`.

    Examples
    --------
    >>> import astropy.units as u
    >>> n = 1e19 * u.m**-3
    >>> T = 1e6 * u.K
    >>> species = ('e', 'p')
    >>> collision_frequency(T, n, species)
    <Quantity 70249... Hz>
    """
    T, masses, charges, reduced_mass, V_r = _process_inputs(T=T, species=species, V=V)
    # using a more descriptive name for the thermal velocity using
    # reduced mass
    V_reduced = V_r

    if species[0] in ("e", "e-") and species[1] in ("e", "e-"):
        # electron-electron collision
        # if a velocity was passed, we use that instead of the reduced
        # thermal velocity
        V = _replace_nan_velocity_with_thermal_velocity(V, T, reduced_mass)
        # impact parameter for 90° collision
        bPerp = impact_parameter_perp(T=T, species=species, V=V_reduced)
    elif species[0] in ("e", "e-") or species[1] in ("e", "e-"):
        # electron-ion collision
        # Need to manually pass electron thermal velocity to obtain
        # correct perpendicular collision radius
        # we ignore the reduced velocity and use the electron thermal
        # velocity instead
        V = _replace_nan_velocity_with_thermal_velocity(V, T, m_e)
        # need to also correct mass in collision radius from reduced
        # mass to electron mass
        bPerp = impact_parameter_perp(T=T, species=species, V=V) * reduced_mass / m_e
        # !!! may also need to correct Coulomb logarithm to be
        # electron-electron version !!!
    else:
        # ion-ion collision
        # if a velocity was passed, we use that instead of the reduced
        # thermal velocity
        V = _replace_nan_velocity_with_thermal_velocity(V, T, reduced_mass)
        bPerp = impact_parameter_perp(T=T, species=species, V=V)
    cou_log = Coulomb_logarithm(T, n, species, z_mean, V=V, method=method)
    # collisional cross section
    sigma = Coulomb_cross_section(bPerp)
    # collision frequency where Coulomb logarithm accounts for
    # small angle collisions, which are more frequent than large
    # angle collisions.
    return n * sigma * V * cou_log


@validate_quantities(impact_param={"can_be_negative": False})
def Coulomb_cross_section(impact_param: u.m) -> u.m**2:
    r"""
    Cross section for a large angle Coulomb collision.

    Parameters
    ----------
    impact_param : `~astropy.units.Quantity`
        Impact parameter for the collision.

    Returns
    -------
    sigma : `~astropy.units.Quantity`
        The Coulomb collision cross section area.

    Notes
    -----
    The `collisional cross-section
    <https://en.wikipedia.org/w/index.php?title=Cross_section_(physics)&oldid=1037954726#Collision_among_gas_particles>`__
    for a 90° Coulomb collision is obtained by

    .. math::

        σ = π (2 * ρ_⟂)^2

    where :math:`ρ_⟂` is the distance of closest approach for a 90°
    Coulomb collision. This function is a generalization of that
    calculation. Please note that it is not guaranteed to return the
    correct results for small angle collisions.

    Examples
    --------
    >>> Coulomb_cross_section(7e-10*u.m)
    <Quantity 6.157...e-18 m2>
    >>> Coulomb_cross_section(0.5*u.m)
    <Quantity 3.141... m2>
    """
    return np.pi * (2 * impact_param) ** 2


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def fundamental_electron_collision_freq(
    T_e: u.K,
    n_e: u.m**-3,
    ion,
    coulomb_log=None,
    V=None,
    coulomb_log_method="classical",
) -> u.s**-1:
    r"""
    Average momentum relaxation rate for a slowly flowing Maxwellian
    distribution of electrons.

    :cite:t:`braginskii:1965` provides a derivation of this as an
    average collision frequency between electrons and ions for a
    Maxwellian distribution. It is thus a special case of the collision
    frequency with an averaging factor, and is on many occasions
    in transport theory the most relevant collision frequency that has
    to be considered. It is heavily related to diffusion and resistivity
    in plasmas.

    Parameters
    ----------
    T_e : `~astropy.units.Quantity`
        The electron temperature of the Maxwellian test electrons.

    n_e : `~astropy.units.Quantity`
        The number density of the Maxwellian test electrons.

    ion : `str`
        String signifying a particle type of the field ions, including
        charge state information.

    V : `~astropy.units.Quantity`, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`μ V^2 \sim 2 k_B T` where
        :math:`μ` is the reduced mass.

    coulomb_log : `float` or dimensionless `~astropy.units.Quantity`, optional
        Option to specify a Coulomb logarithm of the electrons on the
        ions.  If not specified, the Coulomb log will is calculated
        using the `~plasmapy.formulary.collisions.Coulomb_logarithm` function.

    coulomb_log_method : `str`, optional
        The method by which to compute the Coulomb logarithm.  The
        default method is the classical straight-line Landau-Spitzer
        method (``"classical"`` or ``"ls"``). The other 6 supported
        methods are ``"ls_min_interp"``, ``"ls_full_interp"``,
        ``"ls_clamp_mininterp"``, ``"hls_min_interp"``,
        ``"hls_max_interp"``, and ``"hls_full_interp"``.  Please refer
        to the docstring of `Coulomb_logarithm` for more information
        about these methods.

    Returns
    -------
    nu_e : `~astropy.units.Quantity`

    Notes
    -----
    Equations (2.17) and (2.120) in :cite:t:`callen:unpublished` provide
    the original source used to implement this formula, however, the
    simplest form that connects our average collision frequency to the
    general collision frequency is this (from 2.17):

    .. math::

        ν_e = \frac{4}{3 \sqrt{π}} ν(v_{Te})

    Where :math:`ν` is the general collision frequency and
    :math:`v_{Te}` is the electron thermal velocity (the average, for a
    Maxwellian distribution).

    This implementation of the average collision frequency is equivalent to:

    * :math:`1/τ_e` from equation (2.5e) on page 215 of
      :cite:t:`braginskii:1965`
    * :math:`ν_e` from page 33 of :cite:t:`nrlformulary:2019`

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.constants import c
    >>> fundamental_electron_collision_freq(0.1 * u.eV, 1e6 / u.m ** 3, 'p')
    <Quantity 0.001801... 1 / s>
    >>> fundamental_electron_collision_freq(1e6 * u.K, 1e6 / u.m ** 3, 'p')
    <Quantity 1.07221...e-07 1 / s>
    >>> fundamental_electron_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p')
    <Quantity 3935958.7... 1 / s>
    >>> fundamental_electron_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p', coulomb_log_method = 'GMS-1')
    <Quantity 3872815.5... 1 / s>
    >>> fundamental_electron_collision_freq(0.1 * u.eV, 1e6 / u.m ** 3, 'p', V = c / 100)
    <Quantity 5.6589...e-07 1 / s>
    >>> fundamental_electron_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p', coulomb_log = 20)
    <Quantity 5812633... 1 / s>

    See Also
    --------
    ~plasmapy.formulary.collisions.collision_frequency
    ~plasmapy.formulary.collisions.fundamental_ion_collision_freq
    """
    # specify to use electron thermal velocity (most probable), not based on reduced mass
    V = _replace_nan_velocity_with_thermal_velocity(V, T_e, m_e)

    species = [ion, "e-"]
    Z_i = particles.charge_number(ion) * u.dimensionless_unscaled
    n_i = n_e / Z_i
    nu = collision_frequency(
        T_e, n_i, species, z_mean=Z_i, V=V, method=coulomb_log_method
    )
    coeff = 4 / np.sqrt(np.pi) / 3

    # accounting for when a Coulomb logarithm value is passed
    if np.any(coulomb_log):
        cLog = Coulomb_logarithm(
            T_e, n_e, species, z_mean=Z_i, V=V, method=coulomb_log_method
        )
        # dividing out by typical Coulomb logarithm value implicit in
        # the collision frequency calculation and replacing with
        # the user defined Coulomb logarithm value
        nu_mod = nu * coulomb_log / cLog
        return coeff * nu_mod
    else:
        return coeff * nu


@validate_quantities(
    T_i={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_i={"can_be_negative": False},
)
def fundamental_ion_collision_freq(
    T_i: u.K,
    n_i: u.m**-3,
    ion,
    coulomb_log=None,
    V=None,
    coulomb_log_method="classical",
) -> u.s**-1:
    r"""
    Average momentum relaxation rate for a slowly flowing Maxwellian
    distribution of ions.

    :cite:t:`braginskii:1965` provides a derivation of this as an
    average collision frequency between ions and ions for a Maxwellian
    distribution. It is thus a special case of the collision frequency
    with an averaging factor.

    Parameters
    ----------
    T_i : `~astropy.units.Quantity`
        The electron temperature of the Maxwellian test ions.

    n_i : `~astropy.units.Quantity`
        The number density of the Maxwellian test ions.

    ion : `str`
        String signifying a particle type of the test and field ions,
        including charge state information. This function assumes the
        test and field ions are the same species.

    V : `~astropy.units.Quantity`, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`μ V^2 \sim 2 k_B T` where
        :math:`μ` is the reduced mass.

    coulomb_log : `float` or dimensionless `~astropy.units.Quantity`, optional
        Option to specify a Coulomb logarithm of the electrons on the
        ions.  If not specified, the Coulomb log will is calculated
        using the ~plasmapy.formulary.Coulomb_logarithm function.

    coulomb_log_method : `str`, optional
        The method by which to compute the Coulomb logarithm.  The
        default method is the classical straight-line Landau-Spitzer
        method (``"classical"`` or ``"ls"``). The other 6 supported
        methods are ``"ls_min_interp"``, ``"ls_full_interp"``,
        ``"ls_clamp_mininterp"``, ``"hls_min_interp"``,
        ``"hls_max_interp"``, and ``"hls_full_interp"``.  Please refer
        to the docstring of `Coulomb_logarithm` for more information
        about these methods.

    Returns
    -------
    nu_i : `~astropy.units.Quantity`

    Notes
    -----
    Equations (2.36) and (2.122) in :cite:t:`callen:unpublished` provide
    the original source used to implement this formula, however, in our
    implementation we use the very same process that leads to the
    fundamental electron collision rate (2.17), gaining simply a
    different coefficient:

    .. math::

        ν_i = \frac{8}{3 * 4 * \sqrt{π}} ν(v_{Ti})

    Where :math:`ν` is the general collision frequency and
    :math:`v_{Ti}` is the ion thermal velocity (the average, for a
    Maxwellian distribution).

    Note that in the derivation, it is assumed that electrons are
    present in such numbers as to establish quasineutrality, but the
    effects of the test ions colliding with them are not considered
    here. This is a very typical approximation in transport theory.

    This result is an ion momentum relaxation rate, and is used in many
    classical transport expressions. It is equivalent to:

    * :math:`1/τ_i` from equation (2.5i) on page 215 of
      :cite:t:`braginskii:1965`
    * :math:`ν_i` from page 33 of :cite:t:`nrlformulary:2019`

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.constants import c
    >>> fundamental_ion_collision_freq(0.1 * u.eV, 1e6 / u.m ** 3, 'p')
    <Quantity 2.868...e-05 1 / s>
    >>> fundamental_ion_collision_freq(1e6 * u.K, 1e6 / u.m ** 3, 'p')
    <Quantity 1.741...e-09 1 / s>
    >>> fundamental_ion_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p')
    <Quantity 63087.5... 1 / s>
    >>> fundamental_ion_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p', coulomb_log_method='GMS-1')
    <Quantity 63085.1... 1 / s>
    >>> fundamental_ion_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p', V = c / 100)
    <Quantity 9.111... 1 / s>
    >>> fundamental_ion_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p', coulomb_log=20)
    <Quantity 95918.7... 1 / s>

    See Also
    --------
    ~plasmapy.formulary.collisions.collision_frequency
    ~plasmapy.formulary.collisions.fundamental_electron_collision_freq
    """
    m_i = particles.particle_mass(ion)
    species = [ion, ion]

    # specify to use ion thermal velocity (most probable), not based on reduced mass
    V = _replace_nan_velocity_with_thermal_velocity(V, T_i, m_i)

    Z_i = particles.charge_number(ion) * u.dimensionless_unscaled

    nu = collision_frequency(
        T_i, n_i, species, z_mean=Z_i, V=V, method=coulomb_log_method
    )
    # factor of 4 due to reduced mass in bperp and the rest is
    # due to differences in definitions of collisional frequency
    coeff = np.sqrt(8 / np.pi) / 3 / 4

    # accounting for when a Coulomb logarithm value is passed
    if np.any(coulomb_log):
        cLog = Coulomb_logarithm(
            T_i, n_i, species, z_mean=Z_i, V=V, method=coulomb_log_method
        )
        # dividing out by typical Coulomb logarithm value implicit in
        # the collision frequency calculation and replacing with
        # the user defined Coulomb logarithm value
        nu_mod = nu * coulomb_log / cLog
        return coeff * nu_mod
    else:
        return coeff * nu


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def mean_free_path(
    T: u.K,
    n_e: u.m**-3,
    species,
    z_mean: Real = np.nan,
    V: u.m / u.s = np.nan * u.m / u.s,
    method="classical",
) -> u.m:
    r"""
    Collisional mean free path (m).

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and the
        target particle.

    n_e : `~astropy.units.Quantity`
        The electron number density in units convertible to per cubic
        meter.

    species : `tuple`
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second).

    z_mean : `~astropy.units.Quantity`, optional
        The average ionization (arithmetic mean) of a plasma for which a
        macroscopic description is valid. This parameter is used to
        compute the average ion density (given the average ionization
        and electron density) for calculating the ion sphere radius for
        non-classical impact parameters. ``z_mean`` is a required
        parameter if ``method`` is ``"ls_full_interp"``,
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
        `~plasmapy.formulary.collisions.Coulomb_logarithm` for more
        information about these methods.

    Returns
    -------
    mfp : `float` or `numpy.ndarray`
        The collisional mean free path for particles in a plasma.

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
    The collisional mean free path (see :cite:t:`chen:2016`) is given by:

    .. math::

        λ_{mfp} = \frac{v}{ν}

    where :math:`v` is the inter-particle velocity (typically taken to
    be the thermal velocity) and :math:`ν` is the collision frequency.

    Examples
    --------
    >>> import astropy.units as u
    >>> n = 1e19 * u.m ** -3
    >>> T = 1e6 * u.K
    >>> mean_free_path(T, n, ('e-', 'p+'))
    <Quantity 7.839... m>
    >>> mean_free_path(T, n, ('e-', 'p+'), V=1e6 * u.m / u.s)
    <Quantity 0.0109... m>
    """
    # collisional frequency
    freq = collision_frequency(
        T=T, n=n_e, species=species, z_mean=z_mean, V=V, method=method
    )
    # boiler plate to fetch velocity
    # this has been moved to after collision_frequency to avoid use of
    # reduced mass thermal velocity in electron-ion collision case.
    # Should be fine since collision_frequency has its own _process_inputs
    # check, and we are only using this here to get the velocity.
    T, masses, charges, reduced_mass, V = _process_inputs(T=T, species=species, V=V)
    return V / freq


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
        `~plasmapy.formulary.collisions.Coulomb_logarithm` for more
        information about these methods.

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
    freq = collision_frequency(
        T=T, n=n, species=species, z_mean=z_mean, V=V, method=method
    )
    # fetching additional parameters
    T, masses, charges, reduced_mass, V = _process_inputs(T=T, species=species, V=V)
    return (
        freq * reduced_mass / (n * charges[0] * charges[1])
        if np.isnan(z_mean)
        else freq * reduced_mass / (n * (z_mean * e) ** 2)
    )


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
        `~plasmapy.formulary.collisions.Coulomb_logarithm` for more
        information about these methods.

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
    freq = collision_frequency(
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
    n_e={"can_be_negative": False},
)
def Knudsen_number(
    characteristic_length,
    T: u.K,
    n_e: u.m**-3,
    species,
    z_mean: Real = np.nan,
    V: u.m / u.s = np.nan * u.m / u.s,
    method="classical",
) -> u.dimensionless_unscaled:
    r"""
    Knudsen number (dimensionless).

    Parameters
    ----------
    characteristic_length : `~astropy.units.Quantity`
        Rough order-of-magnitude estimate of the relevant size of the
        system.

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
        non-classical impact parameters. ``z_mean`` is a required
        parameter if ``method`` is ``"ls_full_interp"``,
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
        `~plasmapy.formulary.collisions.Coulomb_logarithm` for more
        information about these methods.

    Returns
    -------
    knudsen_param : `float` or `numpy.ndarray`
        The dimensionless Knudsen number.

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
    The `Knudsen number <https://en.wikipedia.org/wiki/Knudsen_number>`_
    is given by

    .. math::

        Kn = \frac{λ_{mfp}}{L}

    where :math:`λ_{mfp}` is the collisional mean free path for
    particles in a plasma and :math:`L` is the characteristic scale
    length of interest.

    The characteristic scale length is typically the plasma size or the
    size of a diagnostic (such as the length or radius of a Langmuir
    probe tip). The Knudsen number tells us whether collisional effects
    are important on this scale length.

    Examples
    --------
    >>> import astropy.units as u
    >>> L = 1e-3 * u.m
    >>> n = 1e19 * u.m ** -3
    >>> T = 1e6 * u.K
    >>> species = ('e', 'p')
    >>> Knudsen_number(L, T, n, species)
    <Quantity 7839.5...>
    >>> Knudsen_number(L, T, n, species, V=1e6 * u.m / u.s)
    <Quantity 10.91773...>
    """
    path_length = mean_free_path(
        T=T, n_e=n_e, species=species, z_mean=z_mean, V=V, method=method
    )
    return path_length / characteristic_length


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def coupling_parameter(
    T: u.K,
    n_e: u.m**-3,
    species,
    z_mean: Real = np.nan,
    V: u.m / u.s = np.nan * u.m / u.s,
    method="classical",
) -> u.dimensionless_unscaled:
    r"""
    Ratio of the Coulomb energy to the kinetic (usually thermal) energy.

    Classical plasmas are weakly coupled (:math:`Γ ≪ 1`, where :math:`Γ`
    is the coupling parameter).  Dense plasmas tend to have significant
    to strong coupling (:math:`Γ ≥ 1`\ ).  For more details, see the
    notes section below.

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and the
        target particle.

    n_e : `~astropy.units.Quantity`
        The electron number density in units convertible to per cubic
        meter.

    species : `tuple`
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second).

    z_mean : `~astropy.units.Quantity`, optional
        The average ionization (arithmetic mean) of a plasma for which a
        macroscopic description is valid. This parameter is used to
        compute the average ion density (given the average ionization
        and electron density) for calculating the ion sphere radius for
        non-classical impact parameters. ``z_mean`` is a required
        parameter if ``method`` is ``"ls_full_interp"``,
        ``"hls_max_interp"``, or ``"hls_full_interp"``.

    V : `~astropy.units.Quantity`, optional
        The relative velocity between particles. If not provided,
        thermal velocity is assumed: :math:`μ V^2 \sim 2 k_B T` where
        :math:`μ` is the reduced mass.

    method : `str`, optional
        The method by which to compute the coupling parameter: either
        ``"classical"`` or ``"quantum"``. The default method is ``"classical"``.
        The Notes section of this docstring has more information about
        these two methods.

    Returns
    -------
    coupling : `float` or `~numpy.ndarray`
        The coupling parameter for a plasma.

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
    The coupling parameter is given by

    .. math::

        Γ = \frac{E_{Coulomb}}{E_{Kinetic}}

    The Coulomb energy is given by

    .. math::

        E_{Coulomb} = \frac{Z_1 Z_2 q_e^2}{4 π ε_0 r}

    where :math:`r` is the Wigner-Seitz radius, and 1 and 2 refer to
    particle species 1 and 2 between which we want to determine the
    coupling.

    In the classical case the kinetic energy is the thermal energy:

    .. math::

        E_{kinetic} = k_B T_e

    The quantum case is more complex. The kinetic energy is dominated by
    the Fermi energy, modulated by a correction factor based on the
    ideal chemical potential. This is obtained more precisely by taking
    the thermal kinetic energy and dividing by the degeneracy
    parameter, modulated by the Fermi integral :cite:p:`gericke:2002`\ :

    .. math::

        E_{kinetic} = 2 k_B T_e / χ f_{3/2} (μ_{ideal} / k_B T_e)

    where :math:`χ` is the degeneracy parameter, :math:`f_{3/2}` is the
    Fermi integral, and :math:`μ_{ideal}` is the ideal chemical
    potential.

    The degeneracy parameter is given by

    .. math::

        χ = n_e Λ_{de Broglie} ^ 3

    where :math:`n_e` is the electron density and :math:`Λ_{de Broglie}`
    is the thermal de Broglie wavelength.

    See equations 1.2, 1.3 and footnote 5 in :cite:t:`bonitz:1998` for
    details on the ideal chemical potential.

    Examples
    --------
    >>> import astropy.units as u
    >>> n = 1e19 * u.m**-3
    >>> T = 1e6 * u.K
    >>> species = ('e', 'p')
    >>> coupling_parameter(T, n, species)
    <Quantity 5.8033...e-05>
    >>> coupling_parameter(T, n, species, V=1e6 * u.m / u.s)
    <Quantity 5.8033...e-05>
    """
    T, masses, charges, reduced_mass, V = _process_inputs(T=T, species=species, V=V)

    if np.isnan(z_mean):
        # using mean charge to get average ion density.
        # If you are running this, you should strongly consider giving
        # a value of z_mean as an argument instead.
        Z1 = np.abs(particles.charge_number(species[0]))
        Z2 = np.abs(particles.charge_number(species[1]))
        Z = (Z1 + Z2) / 2
        # getting ion density from electron density
        n_i = n_e / Z
    else:
        # getting ion density from electron density
        n_i = n_e / z_mean
    # getting Wigner-Seitz radius based on ion density
    radius = Wigner_Seitz_radius(n_i)
    # Coulomb potential energy between particles
    if np.isnan(z_mean):
        coulomb_energy = charges[0] * charges[1] / (4 * np.pi * eps0 * radius)
    else:
        coulomb_energy = (z_mean * e) ** 2 / (4 * np.pi * eps0 * radius)

    if method == "classical":
        # classical thermal kinetic energy
        kinetic_energy = k_B * T
    elif method == "quantum":
        # quantum kinetic energy for dense plasmas
        lambda_deBroglie = thermal_deBroglie_wavelength(T)
        chem_potential = chemical_potential(n_e, T)
        fermi_integral = Fermi_integral(chem_potential.si.value, 1.5)
        denominator = (n_e * lambda_deBroglie**3) * fermi_integral
        kinetic_energy = 2 * k_B * T / denominator
        if np.all(np.imag(kinetic_energy) == 0):
            kinetic_energy = np.real(kinetic_energy)
        else:  # coverage: ignore
            raise ValueError(
                "Kinetic energy should not be imaginary."
                "Something went horribly wrong."
            )
    else:
        raise ValueError(
            f"Keyword 'method' must be either 'classical' or "
            f"'quantum', instead of '{method}'."
        )

    return coulomb_energy / kinetic_energy
