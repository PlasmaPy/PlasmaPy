"""Functions to calculate transport coefficients.

This module includes a number of functions for handling Coulomb collisions
spanning weakly coupled (low density) to strongly coupled (high density)
regimes.

Coulomb collisions
==================

Coulomb collisions are collisions where the interaction force is conveyed
via the electric field, instead of any kind of contact force. They usually
result in relatively small deflections of particle trajectories. However,
given that there are many charged particles in a plasma, one has to take
into account the cumulative effects of many such collisions.

Coulomb logarithms
==================

Please see the documentation for the `Coulomb logarithm <Coulomb_logarithm>`_
for a review of the many ways in which one can define and calculate
that quantity.

Collision rates
===============

The module gathers a few functions helpful for calculating collision
rates between particles. The most general of these is `collision_frequency`,
while if you need average values for a Maxwellian distribution, try
out `collision_rate_electron_ion` and `collision_rate_ion_ion`. These
use `collision_frequency` under the hood.

Macroscopic properties
======================

These include:
* `Spitzer_resistivity`
* `mobility`
* `Knudsen_number`
* `coupling_parameter`

"""


# python modules
from astropy import units as u
import numpy as np
import warnings

# plasmapy modules
from plasmapy import utils
from plasmapy.constants import (c, m_e, k_B, e, eps0, pi, hbar)
from plasmapy import atomic
from plasmapy.physics import parameters
from plasmapy.physics.quantum import (Wigner_Seitz_radius,
                                      thermal_deBroglie_wavelength,
                                      chemical_potential)
from plasmapy.mathematics import Fermi_integral
from plasmapy.utils import check_quantity, _check_relativistic

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


@utils.check_quantity({"T": {"units": u.K, "can_be_negative": False},
                       "n_e": {"units": u.m ** -3}
                       })
def Coulomb_logarithm(T,
                      n_e,
                      particles,
                      z_mean=np.nan * u.dimensionless_unscaled,
                      V=np.nan * u.m / u.s,
                      method="classical"):
    r"""
    Estimates the Coulomb logarithm.

    Parameters
    ----------

    T : ~astropy.units.Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle.

    n_e : ~astropy.units.Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second).

    z_mean : ~astropy.units.Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    V : ~astropy.units.Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
        where `mu` is the reduced mass.

    method: str, optional
        Selects which theory to use when calculating the Coulomb
        logarithm. Defaults to classical method.

    Returns
    -------
    lnLambda : float or numpy.ndarray
        An estimate of the Coulomb logarithm that is accurate to
        roughly its reciprocal.

    Raises
    ------
    ValueError
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    UnitConversionError
        If the units on any of the inputs are incorrect.

        If the n_e, T, or V are not Quantities.

    PhysicsError
        If the result is smaller than 1.

    RelativityError
        If the input velocity is same or greater than the speed
        of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed

    ~plasmapy.utils.RelativityWarning
        If the input velocity is greater than 5% of the speed of
        light.

    Notes
    -----
    The classical Coulomb logarithm is given by

    .. math::
        \ln{\Lambda} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    where :math:`b_{min}` and :math:`b_{max}` are the inner and outer
    impact parameters for Coulomb collisions [1]_.

    The outer impact parameter is given by the Debye length:
    :math:`b_{min} = \lambda_D` which is a function of electron
    temperature and electron density.  At distances greater than the
    Debye length, electric fields from other particles will be
    screened out due to electrons rearranging themselves.

    The choice of inner impact parameter is either the distance of closest
    approach for a 90 degree Coulomb collision or the thermal deBroglie
    wavelength, whichever is larger. This is because Coulomb-style collisions
    cannot occur for impact parameters shorter than the deBroglie
    wavelength because quantum effects will change the fundamental
    nature of the collision [2]_, [3]_.

    Errors associated with the classical Coulomb logarithm are of order its
    inverse. If the Coulomb logarithm is of order unity, then the
    assumptions made in the standard analysis of Coulomb collisions
    are invalid.

    For dense plasmas where the classical Coulomb logarithm breaks
    down there are various extended methods. These can be found
    in D.O. Gericke et al's paper, which has a table summarizing
    the methods [4]_. The GMS-1 through GMS-6 methods correspond
    to the methods found it that table.

    It should be noted that GMS-4 thru GMS-6 modify the Coulomb
    logarithm to the form:

    .. math::
        \ln{\Lambda} \equiv 0.5 \ln\left(1 + \frac{b_{max}^2}{b_{min}^2} \right)

    This means the Coulomb logarithm will not break down for Lambda < 0,
    which occurs for dense, cold plasmas.

    Methods
    ---
    Classical
        classical Landau-Spitzer approach. Fails for large coupling
        parameter where Lambda can become less than zero.
    GMS-1
        1st method listed in Table 1 of reference [3]
        Landau-Spitzer, but with interpolated bmin instead of bmin
        selected between deBroglie wavelength and distance of closest
        approach. Fails for large coupling
        parameter where Lambda can become less than zero.
    GMS-2
        2nd method listed in Table 1 of reference [3]
        Another Landau-Spitzer like approach, but now bmax is also
        being interpolated. The interpolation is between the Debye
        length and the ion sphere radius, allowing for descriptions
        of dilute plasmas. Fails for large coupling
        parameter where Lambda can become less than zero.
        3rd method listed in Table 1 of reference [3]
        classical Landau-Spitzer fails for argument of Coulomb logarithm
        Lambda < 0, therefore a clamp is placed at Lambda_min = 2
    GMS-4
        4th method listed in Table 1 of reference [3]
        Spitzer-like extension to Coulomb logarithm by noting that
        Coulomb collisions take hyperbolic trajectories. Removes
        divergence for small bmin issue in classical Landau-Spitzer
        approach, so bmin can be zero. Also doesn't break down as
        Lambda < 0 is now impossible, even when coupling parameter is large.
    GMS-5
        5th method listed in Table 1 of reference [3]
        Similar to GMS-4, but setting bmin as distance of closest approach
        and bmax interpolated between Debye length and ion sphere radius.
        Lambda < 0 impossible.
    GMS-6
        6th method listed in Table 1 of reference [3]
        Similar to GMS-4 and GMS-5, but using interpolation methods
        for both bmin and bmax.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> Coulomb_logarithm(T, n, particles)
    14.545527226436974
    >>> Coulomb_logarithm(T, n, particles, V=1e6 * u.m / u.s)
    11.363478214139432

    References
    ----------
    .. [1] Physics of Fully Ionized Gases, L. Spitzer (1962)

    .. [2] Francis, F. Chen. Introduction to plasma physics and controlled
       fusion 3rd edition. Ch 5 (Springer 2015).

    .. [3] Comparison of Coulomb Collision Rates in the Plasma Physics
       and Magnetically Confined Fusion Literature, W. Fundamenski and
       O.E. Garcia, EFDA–JET–R(07)01
       (http://www.euro-fusionscipub.org/wp-content/uploads/2014/11/EFDR07001.pdf)

    .. [4] Dense plasma temperature equilibration in the binary collision
       approximation. D. O. Gericke et. al. PRE,  65, 036418 (2002).
       DOI: 10.1103/PhysRevE.65.036418
    """
    # fetching impact min and max impact parameters
    bmin, bmax = impact_parameter(T=T,
                                  n_e=n_e,
                                  particles=particles,
                                  z_mean=z_mean,
                                  V=V,
                                  method=method)
    if method in ["classical", "GMS-1", "GMS-2"]:
        ln_Lambda = np.log(bmax / bmin)
    elif method == "GMS-3":
        ln_Lambda = np.log(bmax / bmin)
        if ln_Lambda < 2:
            ln_Lambda = 2 * u.dimensionless_unscaled
    elif method in ["GMS-4", "GMS-5", "GMS-6"]:
        ln_Lambda = 0.5 * np.log(1 + bmax ** 2 / bmin ** 2)
    else:
        raise ValueError("Unknown method! Choose from 'classical' and 'GMS-N', N from 1 to 6.")
    # applying dimensionless units
    ln_Lambda = ln_Lambda.to(u.dimensionless_unscaled).value
    if ln_Lambda < 2 and method in ["classical", "GMS-1", "GMS-2"]:
        warnings.warn(f"Coulomb logarithm is {ln_Lambda} and {method} relies on weak coupling.",
                      utils.CouplingWarning)
    elif ln_Lambda < 4:
        warnings.warn(f"Coulomb logarithm is {ln_Lambda}, you might have strong coupling effects",
                      utils.CouplingWarning)
    return ln_Lambda


def _boilerPlate(T, particles, V):
    """
    Some boiler plate code for checking if inputs to functions in
    collisions.py are good. Also obtains reduced in mass in a
    2 particle collision system along with thermal velocity.
    """
    # checking temperature is in correct units
    T = T.to(u.K, equivalencies=u.temperature_energy())
    # extracting particle information
    if not isinstance(particles, (list, tuple)) or len(particles) != 2:
        raise ValueError("Particles input must be a "
                         "list or tuple containing representations of two  "
                         f"charged particles. Got {particles} instead.")

    particles = [atomic.Particle(p) for p in particles]
    masses = [p.mass for p in particles]
    charges = [np.abs(p.charge) for p in particles]

    # obtaining reduced mass of 2 particle collision system
    reduced_mass = atomic.reduced_mass(*particles)

    if V == 0:
        raise utils.exceptions.PhysicsError("You cannot have a collision for zero velocity!")
    # getting thermal velocity of system if no velocity is given
    if np.isnan(V):
        V = parameters.thermal_speed(T, mass=reduced_mass)
    _check_relativistic(V, 'V')
    return T, masses, charges, reduced_mass, V


@check_quantity({"T": {"units": u.K, "can_be_negative": False}
                 })
def impact_parameter_perp(T,
                          particles,
                          V=np.nan * u.m / u.s):
    r"""Distance of closest approach for a 90 degree Coulomb collision.

    Parameters
    ----------

    T : ~astropy.units.Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    V : ~astropy.units.Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
        where `mu` is the reduced mass.

    Returns
    -------
    impact_parameter_perp : float or numpy.ndarray
        The distance of closest approach for a 90 degree Coulomb collision.

    Raises
    ------
    ValueError
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    UnitConversionError
        If the units on any of the inputs are incorrect

    TypeError
        If T, or V are not Quantities.

    RelativityError
        If the input velocity is same or greater than the speed
        of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed

    ~plasmapy.utils.RelativityWarning
        If the input velocity is greater than 5% of the speed of
        light.

    Notes
    -----
    The distance of closest approach, impact_parameter_perp, is given by [1]_

    .. math::
        b_{\perp} = \frac{Z_1 Z_2}{4 \pi \epsilon_0 m v^2}


    Examples
    --------
    >>> from astropy import units as u
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> impact_parameter_perp(T, particles)
    <Quantity 8.35505011e-12 m>


    References
    ----------
    .. [1] Francis, F. Chen. Introduction to plasma physics and controlled
       fusion 3rd edition. Ch 5 (Springer 2015).

    """
    # boiler plate checks
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    # Corresponds to a deflection of 90 degrees, which is valid when
    # classical effects dominate.
    # !!!Note: an average ionization parameter will have to be
    # included here in the future
    bPerp = (charges[0] * charges[1] / (4 * pi * eps0 * reduced_mass * V ** 2))
    return bPerp.to(u.m)


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n_e": {"units": u.m ** -3}
                 })
def impact_parameter(T,
                     n_e,
                     particles,
                     z_mean=np.nan * u.dimensionless_unscaled,
                     V=np.nan * u.m / u.s,
                     method="classical"):
    r"""Impact parameters for classical and quantum Coulomb collision

    Parameters
    ----------

    T : ~astropy.units.Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    n_e : ~astropy.units.Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    z_mean : ~astropy.units.Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    V : ~astropy.units.Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
        where `mu` is the reduced mass.

    method: str, optional
        Selects which theory to use when calculating the Coulomb
        logarithm. Defaults to classical method.

    Returns
    -------
    bmin, bmax : tuple of floats
        The minimum and maximum impact parameters (distances) for a
        Coulomb collision.

    Raises
    ------
    ValueError
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    UnitConversionError
        If the units on any of the inputs are incorrect

    TypeError
        If the n_e, T, or V are not Quantities.

    RelativityError
        If the input velocity is same or greater than the speed
        of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed

    ~plasmapy.utils.RelativityWarning
        If the input velocity is greater than 5% of the speed of
        light.

    Notes
    -----
    The minimum and maximum impact parameters may be calculated in a
    variety of ways. The maximum impact parameter is typically
    the Debye length.

    For quantum plasmas the maximum impact parameter can be the
    quadratic sum of the debye length and ion radius (Wigner_Seitz) [1]_

    .. math::
        b_{max} = \left(\lambda_{De}^2 + a_i^2\right)^{1/2}

    The minimum impact parameter is typically some combination of the
    thermal deBroglie wavelength and the distance of closest approach
    for a 90 degree Coulomb collision. A quadratic sum is used for
    all GMS methods, except for GMS-5, where b_min is simply set to
    the distance of closest approach [1]_.

    .. math::
        b_{min} = \left(\Lambda_{deBroglie}^2 + \rho_{\perp}^2\right)^{1/2}

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> impact_parameter(T, n, particles)
    (<Quantity 1.05163088e-11 m>, <Quantity 2.18225522e-05 m>)
    >>> impact_parameter(T, n, particles, V=1e6 * u.m / u.s)
    (<Quantity 2.53401778e-10 m>, <Quantity 2.18225522e-05 m>)

    References
    ----------
    .. [1] Dense plasma temperature equilibration in the binary collision
       approximation. D. O. Gericke et. al. PRE,  65, 036418 (2002).
       DOI: 10.1103/PhysRevE.65.036418
    """
    # boiler plate checks
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    # catching error where mean charge state is not given for non-classical
    # methods that require the ion density
    if method == "GMS-2" or method == "GMS-5" or method == "GMS-6":
        if np.isnan(z_mean):
            raise ValueError("Must provide a z_mean for GMS-2, GMS-5, and "
                             "GMS-6 methods.")
    # Debye length
    lambdaDe = parameters.Debye_length(T, n_e)
    # deBroglie wavelength
    lambdaBroglie = hbar / (2 * reduced_mass * V)
    # distance of closest approach in 90 degree Coulomb collision
    bPerp = impact_parameter_perp(T=T,
                                  particles=particles,
                                  V=V)
    # obtaining minimum and maximum impact parameters depending on which
    # method is requested
    if method == "classical":
        bmax = lambdaDe
        # Coulomb-style collisions will not happen for impact parameters
        # shorter than either of these two impact parameters, so we choose
        # the larger of these two possibilities. That is, between the
        # deBroglie wavelength and the distance of closest approach.
        if bPerp > lambdaBroglie:
            bmin = bPerp
        else:
            bmin = lambdaBroglie
    elif method == "GMS-1":
        # 1st method listed in Table 1 of reference [1]
        # This is just another form of the classical Landau-Spitzer
        # approach, but bmin is interpolated between the deBroglie
        # wavelength and distance of closest approach.
        bmax = lambdaDe
        bmin = (lambdaBroglie ** 2 + bPerp ** 2) ** (1 / 2)
    elif method == "GMS-2":
        # 2nd method listed in Table 1 of reference [1]
        # Another Landau-Spitzer like approach, but now bmax is also
        # being interpolated. The interpolation is between the Debye
        # length and the ion sphere radius, allowing for descriptions
        # of dilute plasmas.
        # Mean ion density.
        n_i = n_e / z_mean
        # mean ion sphere radius.
        ionRadius = Wigner_Seitz_radius(n_i)
        bmax = (lambdaDe ** 2 + ionRadius ** 2) ** (1 / 2)
        bmin = (lambdaBroglie ** 2 + bPerp ** 2) ** (1 / 2)
    elif method == "GMS-3":
        # 3rd method listed in Table 1 of reference [1]
        # same as GMS-1, but not Lambda has a clamp at Lambda_min = 2
        # where Lambda is the argument to the Coulomb logarithm.
        bmax = lambdaDe
        bmin = (lambdaBroglie ** 2 + bPerp ** 2) ** (1 / 2)
    elif method == "GMS-4":
        # 4th method listed in Table 1 of reference [1]
        bmax = lambdaDe
        bmin = (lambdaBroglie ** 2 + bPerp ** 2) ** (1 / 2)
    elif method == "GMS-5":
        # 5th method listed in Table 1 of reference [1]
        # Mean ion density.
        n_i = n_e / z_mean
        # mean ion sphere radius.
        ionRadius = Wigner_Seitz_radius(n_i)
        bmax = (lambdaDe ** 2 + ionRadius ** 2) ** (1 / 2)
        bmin = bPerp
    elif method == "GMS-6":
        # 6th method listed in Table 1 of reference [1]
        # Mean ion density.
        n_i = n_e / z_mean
        # mean ion sphere radius.
        ionRadius = Wigner_Seitz_radius(n_i)
        bmax = (lambdaDe ** 2 + ionRadius ** 2) ** (1 / 2)
        bmin = (lambdaBroglie ** 2 + bPerp ** 2) ** (1 / 2)
    else:
        raise ValueError(f"Method {method} not found!")
    return bmin.to(u.m), bmax.to(u.m)


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n": {"units": u.m ** -3}
                 })
def collision_frequency(T,
                        n,
                        particles,
                        z_mean=np.nan * u.dimensionless_unscaled,
                        V=np.nan * u.m / u.s,
                        method="classical"):
    r"""Collision frequency of particles in a plasma.

    Parameters
    ----------

    T : ~astropy.units.Quantity
        Temperature in units of temperature.
        This should be the electron temperature for electron-electron
        and electron-ion collisions, and the ion temperature for
        ion-ion collisions.


    n : ~astropy.units.Quantity
        The density in units convertible to per cubic meter.
        This should be the electron density for electron-electron collisions,
        and the ion density for electron-ion and ion-ion collisions.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    z_mean : ~astropy.units.Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    V : ~astropy.units.Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
        where `mu` is the reduced mass.

    method: str, optional
        Selects which theory to use when calculating the Coulomb
        logarithm. Defaults to classical method.

    Returns
    -------
    freq : float or numpy.ndarray
        The collision frequency of particles in a plasma.

    Raises
    ------
    ValueError
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    UnitConversionError
        If the units on any of the inputs are incorrect

    TypeError
        If the n_e, T, or V are not Quantities.

    RelativityError
        If the input velocity is same or greater than the speed
        of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed

    ~plasmapy.utils.RelativityWarning
        If the input velocity is greater than 5% of the speed of
        light.

    Notes
    -----
    The collision frequency is given by [1]_

    .. math::
        \nu = n \sigma v \ln{\Lambda}

    where n is the particle density, :math:`\sigma` is the collisional
    cross-section, :math:`v` is the inter-particle velocity (typically
    taken as the thermal velocity), and :math:`\ln{\Lambda}` is the Coulomb
    logarithm accounting for small angle collisions.

    See eq (2.14) in [2]_.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> collision_frequency(T, n, particles)
    <Quantity 702505.15998601 Hz>

    References
    ----------
    .. [1] Francis, F. Chen. Introduction to plasma physics and controlled
       fusion 3rd edition. Ch 5 (Springer 2015).
    .. [2] http://homepages.cae.wisc.edu/~callen/chap2.pdf

    """
    # boiler plate checks
    T, masses, charges, reduced_mass, V_r = _boilerPlate(T=T,
                                                         particles=particles,
                                                         V=V)
    # using a more descriptive name for the thermal velocity using
    # reduced mass
    V_reduced = V_r
    if particles[0] in ('e','e-') and particles[1] in ('e','e-'):
        # if a velocity was passed, we use that instead of the reduced
        # thermal velocity
        if np.isnan(V):
            V = V_reduced
        # electron-electron collision
        # impact parameter for 90 degree collision
        bPerp = impact_parameter_perp(T=T,
                                      particles=particles,
                                      V=V_reduced)
        # Coulomb logarithm
        cou_log = Coulomb_logarithm(T,
                                    n,
                                    particles,
                                    z_mean,
                                    V=np.nan * u.m / u.s,
                                    method=method)
    elif particles[0] in ('e','e-') or particles[1] in ('e','e-'):
        # electron-ion collision
        # Need to manually pass electron thermal velocity to obtain
        # correct perpendicular collision radius
        if np.isnan(V):
            # we ignore the reduced velocity and use the electron thermal
            # velocity instead
            V = np.sqrt(2 * k_B * T / m_e)
        # need to also correct mass in collision radius from reduced
        # mass to electron mass
        bPerp = impact_parameter_perp(T=T,
                                      particles=particles,
                                      V=V) * reduced_mass / m_e
        # Coulomb logarithm
        # !!! may also need to correct Coulomb logarithm to be
        # electron-electron version !!!
        cou_log = Coulomb_logarithm(T,
                                    n,
                                    particles,
                                    z_mean,
                                    V=np.nan * u.m / u.s,
                                    method=method)
    else:
        # if a velocity was passed, we use that instead of the reduced
        # thermal velocity
        if np.isnan(V):
            V = V_reduced
        # ion-ion collision
        bPerp = impact_parameter_perp(T=T,
                                      particles=particles,
                                      V=V)
        # Coulomb logarithm
        cou_log = Coulomb_logarithm(T,
                                    n,
                                    particles,
                                    z_mean,
                                    V=np.nan * u.m / u.s,
                                    method=method)
    # collisional cross section
    sigma = Coulomb_cross_section(bPerp)
    # collision frequency where Coulomb logarithm accounts for
    # small angle collisions, which are more frequent than large
    # angle collisions.
    freq = n * sigma * V * cou_log
    return freq.to(u.Hz)


@check_quantity({
    'impact_param': {'units': u.m,
                    'can_be_negative': False}
    })
def Coulomb_cross_section(impact_param: u.m):
    r"""Cross section for a large angle Coulomb collision.

    Parameters
    ----------
    impact_param : ~astropy.units.Quantity
        Impact parameter for the collision.

    Examples
    --------
    >>> Coulomb_cross_section(7e-10*u.m)
    <Quantity 6.1575216e-18 m2>
    >>> Coulomb_cross_section(0.5*u.m)
    <Quantity 3.14159265 m2>

    Notes
    -----
    The collisional cross-section (see [1]_ for a graphical demonstration)
    for a 90 degree Coulomb collision is obtained by

    .. math::
        \sigma = \pi (2 * \rho_{\perp})^2

    where :math:`\rho_{\perp}` is the distance of closest approach for
    a 90 degree Coulomb collision. This function is a generalization of that
    calculation. Please note that it is not guaranteed to return the correct
    results for small angle collisions.

    Returns
    -------
    ~astropy.units.Quantity
        The Coulomb collision cross section area.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Cross_section_(physics)#Collision_among_gas_particles
    """
    sigma = np.pi * (2 * impact_param) ** 2
    return sigma

@utils.check_quantity({
    'T_e': {'units': u.K, 'can_be_negative': False},
    'n_e': {'units': u.m ** -3, 'can_be_negative': False}
    })
def fundamental_electron_collision_freq(T_e,
                                        n_e,
                                        ion_particle,
                                        coulomb_log=None,
                                        V=None,
                                        coulomb_log_method="classical"):
    r"""
    Average momentum relaxation rate for a slowly flowing Maxwellian distribution of electrons.

    [3]_ provides a derivation of this as an average collision frequency between electrons
    and ions for a Maxwellian distribution. It is thus a special case of the collision
    frequency with an averaging factor, and is on many occasions in transport theory
    the most relevant collision frequency that has to be considered. It is heavily
    related to diffusion and resistivity in plasmas.

    Parameters
    ----------
    T_e : ~astropy.units.Quantity
        The electron temperature of the Maxwellian test electrons

    n_e : ~astropy.units.Quantity
        The number density of the Maxwellian test electrons

    ion_particle: str
        String signifying a particle type of the field ions, including charge
        state information.

    V : ~astropy.units.Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
        where `mu` is the reduced mass.

    coulomb_log : float or dimensionless ~astropy.units.Quantity, optional
        Option to specify a Coulomb logarithm of the electrons on the ions.
        If not specified, the Coulomb log will is calculated using the
        `~plasmapy.physics.transport.Coulomb_logarithm` function.

    coulomb_log_method : string, optional
        Method used for Coulomb logarithm calculation (see that function
        for more documentation). Choose from "classical" or "GMS-1" to "GMS-6".

    Notes
    -----
    Equations (2.17) and (2.120) in [3]_ provide the original source used
    to implement this formula, however, the simplest form that connects our average
    collision frequency to the general collision frequency is is this (from 2.17):

    .. math::
        \nu_e = \frac{4}{3 \sqrt{\pi}} \nu(v_{Te})

    Where :math:`\nu` is the general collision frequency and :math:`v_{Te}`
    is the electron thermal velocity (the average, for a Maxwellian distribution).

    This implementation of the average collision frequency is is equivalent to:
    * 1/tau_e from ref [1]_ eqn (2.5e) pp. 215,
    * nu_e from ref [2]_ pp. 33,

    References
    ----------
    .. [1] Braginskii, S. I. "Transport processes in a plasma." Reviews of
       plasma physics 1 (1965): 205.

    .. [2] Huba, J. D. "NRL (Naval Research Laboratory) Plasma Formulary,
       revised." Naval Research Lab. Report NRL/PU/6790-16-614 (2016).
       https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary

    .. [3] J.D. Callen, Fundamentals of Plasma Physics draft material,
       Chapter 2, http://homepages.cae.wisc.edu/~callen/chap2.pdf

    Examples
    --------
    >>> from astropy import units as u
    >>> fundamental_electron_collision_freq(0.1 * u.eV, 1e6 / u.m ** 3, 'p')
    <Quantity 0.00180172 1 / s>
    >>> fundamental_electron_collision_freq(1e6 * u.K, 1e6 / u.m ** 3, 'p')
    <Quantity 1.07222852e-07 1 / s>
    >>> fundamental_electron_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p')
    <Quantity 3936037.8595928 1 / s>
    >>> fundamental_electron_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p', coulomb_log_method = 'GMS-1')
    <Quantity 3872922.52743562 1 / s>
    >>> fundamental_electron_collision_freq(0.1 * u.eV, 1e6 / u.m ** 3, 'p', V = c / 100)
    <Quantity 4.41166015e-07 1 / s>
    >>> fundamental_electron_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p', coulomb_log = 20)
    <Quantity 5812633.74935003 1 / s>

    See Also
    --------
    collision_frequency
    fundamental_ion_collision_freq
    """
    T_e = T_e.to(u.K, equivalencies=u.temperature_energy())
    if not V:
        # electron thermal velocity (most probable)
        V = np.sqrt(2 * k_B * T_e / m_e)

    particles = [ion_particle, 'e-']
    Z_i = atomic.integer_charge(ion_particle)
    nu = collision_frequency(T_e,
                             n_e,
                             particles,
                             z_mean=Z_i,
                             V=V,
                             method=coulomb_log_method
                             )
    coeff = 4 / np.sqrt(np.pi) / 3


    # accounting for when a Coulomb logarithm value is passed
    if coulomb_log:
        cLog = Coulomb_logarithm(T_e,
                                 n_e,
                                 particles,
                                 z_mean=Z_i,
                                 V=np.nan, # probably needs to be enabled!
                                 method=coulomb_log_method)
        # dividing out by typical Coulomb logarithm value implicit in
        # the collision frequency calculation and replacing with
        # the user defined Coulomb logarithm value
        nu_mod = nu * coulomb_log / cLog
        nu_e = coeff * nu_mod
    else:
        nu_e = coeff * nu
    return nu_e.to(1 / u.s)


@utils.check_quantity({
    'T_i': {'units': u.K, 'can_be_negative': False},
    'n_i': {'units': u.m ** -3, 'can_be_negative': False}
    })
def fundamental_ion_collision_freq(T_i,
                                   n_i,
                                   ion_particle,
                                   coulomb_log=None,
                                   V=None,
                                   coulomb_log_method="classical"):
    r"""
    Average momentum relaxation rate for a slowly flowing Maxwellian distribution of ions.

    [3]_ provides a derivation of this as an average collision frequency between ions
    and ions for a Maxwellian distribution. It is thus a special case of the collision
    frequency with an averaging factor.

    Parameters
    ----------
    T_i : ~astropy.units.Quantity
        The electron temperature of the Maxwellian test ions

    n_i : ~astropy.units.Quantity
        The number density of the Maxwellian test ions

    ion_particle: str
        String signifying a particle type of the test and field ions,
        including charge state information. This function assumes the test
        and field ions are the same species.

    V : ~astropy.units.Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
        where `mu` is the reduced mass.

    coulomb_log : float or dimensionless ~astropy.units.Quantity, optional
        Option to specify a Coulomb logarithm of the electrons on the ions.
        If not specified, the Coulomb log will is calculated using the
        ~plasmapy.physics.transport.Coulomb_logarithm function.

    coulomb_log_method : str, optional
        Method used for Coulomb logarithm calculation (see that function
        for more documentation). Choose from "classical" or "GMS-1" to "GMS-6".

    Notes
    -----
    Equations (2.36) and (2.122) in [3]_ provide the original source used
    to implement this formula, however, in our implementation we use the very
    same process that leads to the fundamental electron collison rate (2.17),
    gaining simply a different coefficient:

    .. math::
        \nu_i = \frac{8}{3 * 4 * \sqrt{\pi}} \nu(v_{Ti})

    Where :math:`\nu` is the general collision frequency and :math:`v_{Ti}`
    is the ion thermal velocity (the average, for a Maxwellian distribution).

    Note that in the derivation, it is assumed that electrons are present
    in such numbers as to establish quasineutrality, but the effects of the
    test ions colliding with them are not considered here. This is a very
    typical approximation in transport theory.

    This result is an ion momentum relaxation rate, and is used in many
    classical transport expressions. It is equivalent to:
    * 1/tau_i from ref [1]_, equation (2.5i) pp. 215,
    * nu_i from ref [2]_ pp. 33,


    References
    ----------
    .. [1] Braginskii, S. I. "Transport processes in a plasma." Reviews of
       plasma physics 1 (1965): 205.

    .. [2] Huba, J. D. "NRL (Naval Research Laboratory) Plasma Formulary,
       revised." Naval Research Lab. Report NRL/PU/6790-16-614 (2016).
       https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary

    .. [3] J.D. Callen, Fundamentals of Plasma Physics draft material,
       Chapter 2, http://homepages.cae.wisc.edu/~callen/chap2.pdf

    Examples
    --------
    >>> from astropy import units as u
    >>> fundamental_ion_collision_freq(0.1 * u.eV, 1e6 / u.m ** 3, 'p')
    <Quantity 2.97315582e-05 1 / s>
    >>> fundamental_ion_collision_freq(1e6 * u.K, 1e6 / u.m ** 3, 'p')
    <Quantity 1.78316012e-09 1 / s>
    >>> fundamental_ion_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p')
    <Quantity 66411.80316364 1 / s>
    >>> fundamental_ion_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p', coulomb_log_method='GMS-1')
    <Quantity 66407.00859126 1 / s>
    >>> fundamental_ion_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p', V = c / 100)
    <Quantity 6.53577473 1 / s>
    >>> fundamental_ion_collision_freq(100 * u.eV, 1e20 / u.m ** 3, 'p', coulomb_log=20)
    <Quantity 95918.76240877 1 / s>

    See Also
    --------
    collision_frequency
    fundamental_electron_collision_freq
    """
    T_i = T_i.to(u.K, equivalencies=u.temperature_energy())
    m_i = atomic.particle_mass(ion_particle)
    particles = [ion_particle, ion_particle]
    if not V:
        # ion thermal velocity (most probable)
        V = np.sqrt(2 * k_B * T_i / m_i)
    Z_i = atomic.integer_charge(ion_particle)
    nu = collision_frequency(T_i,
                             n_i,
                             particles,
                             z_mean=Z_i,
                             V=V,
                             method=coulomb_log_method)
    # factor of 4 due to reduced mass in bperp and the rest is
    # due to differences in definitions of collisional frequency
    coeff = np.sqrt(8 / np.pi) / 3 / 4

    # accounting for when a Coulomb logarithm value is passed
    if coulomb_log:
        cLog = Coulomb_logarithm(T_i,
                                 n_i,
                                 particles,
                                 z_mean=Z_i,
                                 V=np.nan * u.m / u.s,
                                 method=coulomb_log_method)
        # dividing out by typical Coulomb logarithm value implicit in
        # the collision frequency calculation and replacing with
        # the user defined Coulomb logarithm value
        nu_mod = nu * coulomb_log / cLog
        nu_i = coeff * nu_mod
    else:
        nu_i = coeff * nu
    return nu_i.to(1 / u.s)


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n_e": {"units": u.m ** -3}
                 })
def mean_free_path(T,
                   n_e,
                   particles,
                   z_mean=np.nan * u.dimensionless_unscaled,
                   V=np.nan * u.m / u.s,
                   method="classical"):
    r"""Collisional mean free path (m)

    Parameters
    ----------

    T : ~astropy.units.Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    n_e : ~astropy.units.Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    z_mean : ~astropy.units.Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    V : ~astropy.units.Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
        where `mu` is the reduced mass.

    method: str, optional
        Selects which theory to use when calculating the Coulomb
        logarithm. Defaults to classical method.

    Returns
    -------
    mfp : float or numpy.ndarray
        The collisional mean free path for particles in a plasma.

    Raises
    ------
    ValueError
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    UnitConversionError
        If the units on any of the inputs are incorrect

    TypeError
        If the n_e, T, or V are not Quantities.

    RelativityError
        If the input velocity is same or greater than the speed
        of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed

    ~plasmapy.utils.RelativityWarning
        If the input velocity is greater than 5% of the speed of
        light.

    Notes
    -----
    The collisional mean free path is given by [1]_

    .. math::
        \lambda_{mfp} = \frac{v}{\nu}

    where :math:`v` is the inter-particle velocity (typically taken to be
    the thermal velocity) and :math:`\nu` is the collision frequency.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> mean_free_path(T, n, particles)
    <Quantity 7.8393631 m>
    >>> mean_free_path(T, n, particles, V=1e6 * u.m / u.s)
    <Quantity 0.00852932 m>

    References
    ----------
    .. [1] Francis, F. Chen. Introduction to plasma physics and controlled
       fusion 3rd edition. Ch 5 (Springer 2015).
    """
    # collisional frequency
    freq = collision_frequency(T=T,
                               n=n_e,
                               particles=particles,
                               z_mean=z_mean,
                               V=V,
                               method=method)
    # boiler plate to fetch velocity
    # this has been moved to after collision_frequency to avoid use of
    # reduced mass thermal velocity in electron-ion collision case.
    # Should be fine since collision_frequency has its own boiler_plate
    # check, and we are only using this here to get the velocity.
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    # mean free path length
    mfp = V / freq
    return mfp.to(u.m)


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n": {"units": u.m ** -3}
                 })
def Spitzer_resistivity(T,
                        n,
                        particles,
                        z_mean=np.nan * u.dimensionless_unscaled,
                        V=np.nan * u.m / u.s,
                        method="classical"):
    r"""Spitzer resistivity of a plasma

    Parameters
    ----------

    T : ~astropy.units.Quantity
        Temperature in units of temperature.
        This should be the electron temperature for electron-electron
        and electron-ion collisions, and the ion temperature for
        ion-ion collisions.

    n : ~astropy.units.Quantity
        The density in units convertible to per cubic meter.
        This should be the electron density for electron-electron collisions,
        and the ion density for electron-ion and ion-ion collisions.

    z_mean : ~astropy.units.Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    V : ~astropy.units.Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
        where `mu` is the reduced mass.

    method: str, optional
        Selects which theory to use when calculating the Coulomb
        logarithm. Defaults to classical method.

    Returns
    -------
    spitzer : float or numpy.ndarray
        The resistivity of the plasma in Ohm meters.

    Raises
    ------
    ValueError
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    UnitConversionError
        If the units on any of the inputs are incorrect

    TypeError
        If the n_e, T, or V are not Quantities.

    RelativityError
        If the input velocity is same or greater than the speed
        of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed

    ~plasmapy.utils.RelativityWarning
        If the input velocity is greater than 5% of the speed of
        light.

    Notes
    -----
    The Spitzer resistivity is given by [1]_ [2]_

    .. math::
        \eta = \frac{m}{n Z_1 Z_2 q_e^2} \nu_{1,2}

    where :math:`m` is the ion mass or the reduced mass, :math:`n` is the
    ion density, :math:`Z` is the particle charge state, :math:`q_e` is the
    charge of an electron, :math:`\nu_{1,2}` is the collisional frequency
    between particle species 1 and 2.

    Typically, particle species 1 and 2 are selected to be an electron
    and an ion, since electron-ion collisions are inelastic and therefore
    produce resistivity in the plasma.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> Spitzer_resistivity(T, n, particles)
    <Quantity 2.4916169e-06 m Ohm>
    >>> Spitzer_resistivity(T, n, particles, V=1e6 * u.m / u.s)
    <Quantity 0.00041583 m Ohm>

    References
    ----------
    .. [1] Francis, F. Chen. Introduction to plasma physics and controlled
       fusion 3rd edition. Ch 5 (Springer 2015).
    .. [2] http://homepages.cae.wisc.edu/~callen/chap2.pdf

    """
    # collisional frequency
    freq = collision_frequency(T=T,
                               n=n,
                               particles=particles,
                               z_mean=z_mean,
                               V=V,
                               method=method)
    # boiler plate checks
    # fetching additional parameters
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    if np.isnan(z_mean):
        spitzer = freq * reduced_mass / (n * charges[0] * charges[1])
    else:
        spitzer = freq * reduced_mass / (n * (z_mean * e) ** 2)
    return spitzer.to(u.Ohm * u.m)


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n_e": {"units": u.m ** -3}
                 })
def mobility(T,
             n_e,
             particles,
             z_mean=np.nan * u.dimensionless_unscaled,
             V=np.nan * u.m / u.s,
             method="classical"):
    r"""Electrical mobility (m^2/(V s))

    Parameters
    ----------

    T : ~astropy.units.Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    n_e : ~astropy.units.Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    z_mean : ~astropy.units.Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters. It is also used the obtain the average mobility
        of a plasma with multiple charge state species. When z_mean
        is not given, the average charge between the two particles is
        used instead.

    V : ~astropy.units.Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
        where `mu` is the reduced mass.

    method: str, optional
        Selects which theory to use when calculating the Coulomb
        logarithm. Defaults to classical method.

    Returns
    -------
    mobility_value : float or numpy.ndarray
        The electrical mobility of particles in a collisional plasma.

    Raises
    ------
    ValueError
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    UnitConversionError
        If the units on any of the inputs are incorrect

    TypeError
        If the n_e, T, or V are not Quantities.

    RelativityError
        If the input velocity is same or greater than the speed
        of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed

    ~plasmapy.utils.RelativityWarning
        If the input velocity is greater than 5% of the speed of
        light.

    Notes
    -----
    The mobility is given by [1]_

    .. math::
        \mu = \frac{q}{m \nu}

    where :math:`q` is the particle charge, :math:`m` is the particle mass
    and :math:`\nu` is the collisional frequency of the particle in the
    plasma.

    The mobility describes the forced diffusion of a particle in a collisional
    plasma which is under the influence of an electric field. The mobility
    is essentially the ratio of drift velocity due to collisions and the
    electric field driving the forced diffusion.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> mobility(T, n, particles)
    <Quantity 250500.35318738 m2 / (s V)>
    >>> mobility(T, n, particles, V=1e6 * u.m / u.s)
    <Quantity 1500.97042427 m2 / (s V)>

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Electrical_mobility#Mobility_in_gas_phase
    """
    freq = collision_frequency(T=T,
                               n=n_e,
                               particles=particles,
                               z_mean=z_mean,
                               V=V,
                               method=method)
    # boiler plate checks
    # we do this after collision_frequency since collision_frequency
    # already has a boiler_plate check and we are doing this just
    # to recover the charges, mass, etc.
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    if np.isnan(z_mean):
        z_val = (charges[0] + charges[1]) / 2
    else:
        z_val = z_mean * e
    mobility_value = z_val / (reduced_mass * freq)
    return mobility_value.to(u.m ** 2 / (u.V * u.s))


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n_e": {"units": u.m ** -3}
                 })
def Knudsen_number(characteristic_length,
                   T,
                   n_e,
                   particles,
                   z_mean=np.nan * u.dimensionless_unscaled,
                   V=np.nan * u.m / u.s,
                   method="classical"):
    r"""Knudsen number (dimless)

    Parameters
    ----------

    characteristic_length : ~astropy.units.Quantity
        Rough order-of-magnitude estimate of the relevant size of the system.

    T : ~astropy.units.Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    n_e : ~astropy.units.Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    z_mean : ~astropy.units.Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    V : ~astropy.units.Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
        where `mu` is the reduced mass.

    method: str, optional
        Selects which theory to use when calculating the Coulomb
        logarithm. Defaults to classical method.

    Returns
    -------
    knudsen_param : float or numpy.ndarray
        The dimensionless Knudsen number.

    Raises
    ------
    ValueError
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    UnitConversionError
        If the units on any of the inputs are incorrect

    TypeError
        If the n_e, T, or V are not Quantities.

    RelativityError
        If the input velocity is same or greater than the speed
        of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed

    ~plasmapy.utils.RelativityWarning
        If the input velocity is greater than 5% of the speed of
        light.

    Notes
    -----
    The Knudsen number is given by [1]_

    .. math::
        Kn = \frac{\lambda_{mfp}}{L}

    where :math:`\lambda_{mfp}` is the collisional mean free path for
    particles in a plasma and :math`L` is the characteristic scale
    length of interest.

    Typically the characteristic scale length is the plasma size or the
    size of a diagnostic (such a the length or radius of a Langmuir
    probe tip). The Knudsen number tells us whether collisional effects
    are important on this scale length.

    Examples
    --------
    >>> from astropy import units as u
    >>> L = 1e-3 * u.m
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> Knudsen_number(L, T, n, particles)
    <Quantity 7839.36310417>
    >>> Knudsen_number(L, T, n, particles, V=1e6 * u.m / u.s)
    <Quantity 8.52931736>

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Knudsen_number

    """
    path_length = mean_free_path(T=T,
                                 n_e=n_e,
                                 particles=particles,
                                 z_mean=z_mean,
                                 V=V,
                                 method=method)
    knudsen_param = path_length / characteristic_length
    return knudsen_param.to(u.dimensionless_unscaled)


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n_e": {"units": u.m ** -3}
                 })
def coupling_parameter(T,
                       n_e,
                       particles,
                       z_mean=np.nan * u.dimensionless_unscaled,
                       V=np.nan * u.m / u.s,
                       method="classical"):
    r"""Coupling parameter.
    Coupling parameter compares Coulomb energy to kinetic energy (typically)
    thermal. Classical plasmas are weakly coupled Gamma << 1, whereas dense
    plasmas tend to have significant to strong coupling Gamma >= 1.

    Parameters
    ----------

    T : ~astropy.units.Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    n_e : ~astropy.units.Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    z_mean : ~astropy.units.Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    V : ~astropy.units.Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
        where `mu` is the reduced mass.

    method: str, optional
        Selects which theory to use when calculating the Coulomb
        logarithm. Defaults to classical method.

    Returns
    -------
    coupling : float or numpy.ndarray
        The coupling parameter for a plasma.

    Raises
    ------
    ValueError
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    UnitConversionError
        If the units on any of the inputs are incorrect

    TypeError
        If the n_e, T, or V are not Quantities.

    RelativityError
        If the input velocity is same or greater than the speed
        of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed

    ~plasmapy.utils.RelativityWarning
        If the input velocity is greater than 5% of the speed of
        light.

    Notes
    -----
    The coupling parameter is given by

    .. math::
        \Gamma = \frac{E_{Coulomb}}{E_{Kinetic}}

    The Coulomb energy is given by

    .. math::
        E_{Coulomb} = \frac{Z_1 Z_2 q_e^2}{4 \pi \epsilon_0 r}

    where :math:`r` is the Wigner-Seitz radius, and 1 and 2 refer to
    particle species 1 and 2 between which we want to determine the
    coupling.

    In the classical case the kinetic energy is simply the thermal energy

    .. math::
        E_{kinetic} = k_B T_e

    The quantum case is more complex. The kinetic energy is dominated by
    the Fermi energy, modulated by a correction factor based on the
    ideal chemical potential. This is obtained more precisely
    by taking the the thermal kinetic energy and dividing by
    the degeneracy parameter, modulated by the Fermi integral [1]_

    .. math::
        E_{kinetic} = 2 k_B T_e / \chi f_{3/2} (\mu_{ideal} / k_B T_e)

    where :math:`\chi` is the degeneracy parameter, :math:`f_{3/2}` is the
    Fermi integral, and :math:`\mu_{ideal}` is the ideal chemical
    potential.

    The degeneracy parameter is given by

    .. math::
        \chi = n_e \Lambda_{deBroglie} ^ 3

    where :math:`n_e` is the electron density and :math:`\Lambda_{deBroglie}`
    is the thermal deBroglie wavelength.

    See equations 1.2, 1.3 and footnote 5 in [2]_ for details on the ideal
    chemical potential.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> coupling_parameter(T, n, particles)
    <Quantity 5.80330315e-05>
    >>> coupling_parameter(T, n, particles, V=1e6 * u.m / u.s)
    <Quantity 5.80330315e-05>

    References
    ----------
    .. [1] Dense plasma temperature equilibration in the binary collision
       approximation. D. O. Gericke et. al. PRE,  65, 036418 (2002).
       DOI: 10.1103/PhysRevE.65.036418
    .. [2] Bonitz, Michael. Quantum kinetic theory. Stuttgart: Teubner, 1998.


    """
    # boiler plate checks
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    if np.isnan(z_mean):
        # using mean charge to get average ion density.
        # If you are running this, you should strongly consider giving
        # a value of z_mean as an argument instead.
        Z1 = np.abs(atomic.integer_charge(particles[0]))
        Z2 = np.abs(atomic.integer_charge(particles[1]))
        Z = (Z1 + Z2) / 2
        # getting ion density from electron density
        n_i = n_e / Z
        # getting Wigner-Seitz radius based on ion density
        radius = Wigner_Seitz_radius(n_i)
    else:
        # getting ion density from electron density
        n_i = n_e / z_mean
        # getting Wigner-Seitz radius based on ion density
        radius = Wigner_Seitz_radius(n_i)
    # Coulomb potential energy between particles
    if np.isnan(z_mean):
        coulombEnergy = charges[0] * charges[1] / (4 * np.pi * eps0 * radius)
    else:
        coulombEnergy = (z_mean * e) ** 2 / (4 * np.pi * eps0 * radius)
    if method == "classical":
        # classical thermal kinetic energy
        kineticEnergy = k_B * T
    elif method == "quantum":
        # quantum kinetic energy for dense plasmas
        lambda_deBroglie = thermal_deBroglie_wavelength(T)
        chemicalPotential = chemical_potential(n_e, T)
        fermiIntegral = Fermi_integral(chemicalPotential.si.value, 1.5)
        denom = (n_e * lambda_deBroglie ** 3) * fermiIntegral
        kineticEnergy = 2 * k_B * T / denom
        if np.imag(kineticEnergy) == 0:
            kineticEnergy = np.real(kineticEnergy)
        else:
            raise ValueError("Kinetic energy should not be imaginary."
                             "Something went horribly wrong.")
    coupling = coulombEnergy / kineticEnergy
    return coupling.to(u.dimensionless_unscaled)
