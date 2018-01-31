# coding=utf-8
"""Functions to calculate transport coefficients."""

from astropy import units as u
import numpy as np
import plasmapy.atomic as atomic
from plasmapy import utils
from plasmapy.utils.checks import check_quantity, check_relativistic
from plasmapy.utils.exceptions import PhysicsError, PhysicsWarning
from plasmapy.constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, h, hbar)
from ..atomic import (ion_mass, integer_charge)
from plasmapy.atomic.atomic import _is_electron
from .parameters import (Debye_length, Hall_parameter,
                         collision_rate_electron_ion, collision_rate_ion_ion)
from .quantum import (Wigner_Seitz_radius,
                      thermal_deBroglie_wavelength,
                      chemical_potential)
from ..mathematics import Fermi_integral
from inspect import stack
from copy import copy
import warnings


@utils.check_quantity({"T": {"units": u.K, "can_be_negative": False},
                       "n_e": {"units": u.m**-3}
                       })
def Coulomb_logarithm(T,
                      n_e,
                      particles,
                      z_mean=np.nan*u.dimensionless_unscaled,
                      V=np.nan*u.m/u.s,
                      method="classical"):
    r"""
    Estimates the Coulomb logarithm.

    Parameters
    ----------

    T : Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle.

    n_e : Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second).

    z_mean : Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    V : Quantity, optional
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

    UserWarning
        If the input velocity is greater than 80% of the speed of
        light.

    TypeError
        If the n_e, T, or V are not Quantities.

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

    The choice of inner impact parameter is the distance of closest
    approach for a 90 degree Coulomb collision [2]_, [3]_.
    
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
        
    This means the Coulomb logarithm will not breakdown for Lambda < 0,
    which occurs for dense, cold plasmas.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> Coulomb_logarithm(T, n, particles)
    14.545527226436974
    >>> Coulomb_logarithm(T, n, particles, V=1e6*u.m/u.s)
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
    if method == "classical":
        # classical Landau-Spitzer approach. Fails for large coupling
        # parameter where Lambda can become less than zero.
        ln_Lambda = np.log(bmax / bmin)
    elif method == "GMS-1":
        # 1st method listed in Table 1 of reference [3]
        # Landau-Spitzer, but with interpolated bmin instead of bmin
        # selected between deBroglie wavelength and distance of closest
        # approach. Fails for large coupling
        # parameter where Lambda can become less than zero.
        ln_Lambda = np.log(bmax / bmin)
    elif method == "GMS-2":
        # 2nd method listed in Table 1 of reference [3]
        # Another Landau-Spitzer like approach, but now bmax is also
        # being interpolated. The interpolation is between the Debye
        # length and the ion sphere radius, allowing for descriptions
        # of dilute plasmas. Fails for large coupling
        # parameter where Lambda can become less than zero.
        ln_Lambda = np.log(bmax / bmin)
    elif method == "GMS-3":
        # 3rd method listed in Table 1 of reference [3]
        # classical Landau-Spitzer fails for argument of Coulomb logarithm
        # Lambda < 0, therefore a clamp is placed at Lambda_min = 2
        ln_Lambda = np.log(bmax / bmin)
        if ln_Lambda < 2:
            ln_Lambda = 2
    elif method == "GMS-4":
        # 4th method listed in Table 1 of reference [3]
        # Spitzer-like extension to Coulomb logarithm by noting that
        # Coulomb collisions take hyperbolic trajectories. Removes
        # divergence for small bmin issue in classical Landau-Spitzer
        # approach, so bmin can be zero. Also doesn't breakdown as
        # Lambda < 0 is now impossible, even when coupling parameter is large.
        ln_Lambda = 0.5 * np.log(1 + bmax ** 2 / bmin ** 2)
    elif method == "GMS-5":
        # 5th method listed in Table 1 of reference [3]
        # Similar to GMS-4, but setting bmin as distance of closest approach
        # and bmax interpolated between Debye length and ion sphere radius.
        # Lambda < 0 impossible.
        ln_Lambda = 0.5 * np.log(1 + bmax ** 2 / bmin ** 2)
    elif method == "GMS-6":
        # 6th method listed in Table 1 of reference [3]
        # Similar to GMS-4 and GMS-5, but using interpolation methods
        # for both bmin and bmax.
        # Lambda < 0 impossible.
        ln_Lambda = 0.5 * np.log(1 + bmax ** 2 / bmin ** 2)
    # applying dimensionless units
    ln_Lambda = ln_Lambda.to(u.dimensionless_unscaled).value
    return ln_Lambda

#@check_relativistic
def _boilerPlate(T, particles, V):
    """
    Some boiler plate code for checking if inputs to functions in
    transport.py are good. Also obtains reduced in mass in a
    2 particle collision system along with thermal velocity.
    """
    # checking temperature is in correct units
    T = T.to(u.K, equivalencies=u.temperature_energy())
    # extracting particle information
    if not isinstance(particles, (list, tuple)) or len(particles) != 2:
        raise ValueError("Particles input must be a "
                         "list or tuple containing representations of two  "
                         f"charged particles. Got {particles} instead.")

    masses = np.zeros(2) * u.kg
    charges = np.zeros(2) * u.C

    for particle, i in zip(particles, range(2)):

        try:
            masses[i] = ion_mass(particles[i])
        except Exception:
            raise ValueError("Unable to find mass of particle: "
                             f"{particles[i]}.")
        try:
            charges[i] = np.abs(e * integer_charge(particles[i]))
            if charges[i] is None:
                raise ValueError("Unable to find charge of particle: "
                                 f"{particles[i]}.")
        except Exception:
            raise ValueError("Unable to find charge of particle: "
                             f"{particles[i]}.")
    # obtaining reduced mass of 2 particle collision system
    reduced_mass = masses[0] * masses[1] / (masses[0] + masses[1])
    # getting thermal velocity of system if no velocity is given
    if np.isnan(V):
        V = np.sqrt(2 * k_B * T / reduced_mass)
#    else:
#        _check_relativistic(V, 'impact_parameter', betafrac=0.8)
    return (T, masses, charges, reduced_mass, V)


@check_quantity({"T": {"units": u.K, "can_be_negative": False}
                 })
def b_perp(T,
           particles,
           V=np.nan*u.m/u.s):
    """
    Distance of closest approach for a 90 degree Coulomb collision
    """
    # boiler plate checks
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    # Corresponds to a deflection of 90 degrees, which is valid when
    # classical effects dominate.
    # !!!Note: an average ionization parameter will have to be
    # included here in the future
    bPerp = (charges[0] * charges[1] /
             (4 * pi * eps0 * reduced_mass * V ** 2))
    return bPerp.to(u.m)


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n_e": {"units": u.m**-3}
                 })
def impact_parameter(T,
                     n_e,
                     particles,
                     z_mean=np.nan*u.dimensionless_unscaled,
                     V=np.nan*u.m/u.s,
                     method="classical"):
    r"""Impact parameter for classical Coulomb collision

    Parameters
    ----------

    T : Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    n_e : Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    z_mean : Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    V : Quantity, optional
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
        If the units on any of the inputs are incorrect

    UserWarning
        If the inputted velocity is greater than 80% of the speed of
        light.

    TypeError
        If the n_e, T, or V are not Quantities.

    Notes
    -----
    The Coulomb logarithm is given by

    .. math::
        \ln{\Lambda} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    where :math:`b_{min}` and :math:`b_{max}` are the inner and outer
    impact parameters for Coulomb collisions [1]_.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> impact_parameter(T, n, particles)
    (<Quantity 1.05163088e-11 m>, <Quantity 2.18225522e-05 m>)
    >>> impact_parameter(T, n, particles, V=1e6*u.m/u.s)
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
    lambdaDe = Debye_length(T, n_e)
    # deBroglie wavelength
    lambdaDB = hbar / (2 * reduced_mass * V)
    # distance of closest approach in 90 degree Coulomb collision
    bPerp = b_perp(T=T,
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
        if bPerp > lambdaDB:
            bmin = bPerp
        else:
            bmin = lambdaDB
    elif method == "GMS-1":
        # 1st method listed in Table 1 of reference [1]
        # This is just another form of the classical Landau-Spitzer
        # approach, but bmin is interpolated between the deBroglie
        # wavelength and distance of closest approach.
        bmax = lambdaDe
        bmin = (lambdaDB ** 2 + bPerp ** 2) ** (1 / 2)
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
        bmin = (lambdaDB ** 2 + bPerp ** 2) ** (1 / 2)
    elif method == "GMS-3":
        # 3rd method listed in Table 1 of reference [1]
        # same as GMS-1, but not Lambda has a clamp at Lambda_min = 2
        # where Lambda is the argument to the Coulomb logarithm.
        bmax = lambdaDe
        bmin = (lambdaDB ** 2 + bPerp ** 2) ** (1 / 2)
    elif method == "GMS-4":
        # 4th method listed in Table 1 of reference [1]
        bmax = lambdaDe
        bmin = (lambdaDB ** 2 + bPerp ** 2) ** (1 / 2)
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
        bmin = (lambdaDB ** 2 + bPerp ** 2) ** (1 / 2)
    else:
        raise ValueError(f"Method {method} not found!")
    return (bmin.to(u.m), bmax.to(u.m))


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n": {"units": u.m**-3}
                 })
def collision_frequency(T,
                        n,
                        particles,
                        z_mean=np.nan*u.dimensionless_unscaled,
                        V=np.nan*u.m/u.s,
                        method="classical"):
    r"""Collision frequency of particles in a plasma.

    Parameters
    ----------

    T : Quantity
        Temperature in units of temperature.
        This should be the electron temperature for electron-electron
        and electron-ion collisions, and the ion temperature for
        ion-ion collisions.


    n : Quantity
        The density in units convertible to per cubic meter.
        This should be the electron density for electron-electron collisions,
        and the ion density for electron-ion and ion-ion collisions.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    z_mean : Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    V : Quantity, optional
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
        If the units on any of the inputs are incorrect

    UserWarning
        If the inputted velocity is greater than 80% of the speed of
        light.

    TypeError
        If the n_e, T, or V are not Quantities.

    Notes
    -----
    The Coulomb logarithm is given by

    .. math::
        \ln{\Lambda} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    where :math:`b_{min}` and :math:`b_{max}` are the inner and outer
    impact parameters for Coulomb collisions.

    See eq (2.14) in [1].

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


    """
    # boiler plate checks
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    if particles[0] == 'e' and particles[1] == 'e':
        # electron-electron collision
        # impact parameter for 90 degree collision
        bPerp = b_perp(T=T,
                       particles=particles,
                       V=V)
        # Coulomb logarithm
        cou_log = Coulomb_logarithm(T,
                                    n,
                                    particles,
                                    z_mean,
                                    V=np.nan*u.m/u.s,
                                    method=method)
    elif particles[0] == 'e' or particles[1] == 'e':
        # electron-ion collision
        # Need to manually pass electron thermal velocity to obtain
        # correct perpendicular collision radius
        V = np.sqrt(2 * k_B * T / m_e)
        # need to also correct mass in collision radius from reduced
        # mass to electron mass
        bPerp = b_perp(T=T,
                       particles=particles,
                       V=V) * reduced_mass / m_e
        # Coulomb logarithm
        #!!! may also need to correct Coulomb logarithm to be
        # electron-electron version !!!
        cou_log = Coulomb_logarithm(T,
                                    n,
                                    particles,
                                    z_mean,
                                    V=np.nan*u.m/u.s,
                                    method=method)
    else:
        # ion-ion collision
        bPerp = b_perp(T=T,
                       particles=particles,
                       V=V)
        # Coulomb logarithm
        cou_log = Coulomb_logarithm(T,
                                    n,
                                    particles,
                                    z_mean,
                                    V=np.nan*u.m/u.s,
                                    method=method)
    # collisional cross section
    sigma = np.pi * (2 * bPerp) ** 2
    # collision frequency where Coulomb logarithm accounts for
    # small angle collisions, which are more frequent than large
    # angle collisions.
    freq =  n * sigma * V * cou_log
    return freq.to(u.Hz)


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n_e": {"units": u.m**-3}
                 })
def mean_free_path(T,
                   n_e,
                   particles,
                   z_mean=np.nan*u.dimensionless_unscaled,
                   V=np.nan*u.m/u.s,
                   method="classical"):
    r"""Collisional mean free path (m)

    Parameters
    ----------

    T : Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    n_e : Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    z_mean : Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    V : Quantity, optional
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
        If the units on any of the inputs are incorrect

    UserWarning
        If the inputted velocity is greater than 80% of the speed of
        light.

    TypeError
        If the n_e, T, or V are not Quantities.

    Notes
    -----
    The Coulomb logarithm is given by

    .. math::
        \ln{\Lambda} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    where :math:`b_{min}` and :math:`b_{max}` are the inner and outer
    impact parameters for Coulomb collisions.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> mean_free_path(T, n, particles)
    <Quantity 7.8393631 m>
    >>> mean_free_path(T, n, particles, V=1e6*u.m/u.s)
    <Quantity 1.42347709 m>

    References
    ----------

    """
    # boiler plate checks
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    # collisional frequency
    freq = collision_frequency(T=T,
                               n=n_e,
                               particles=particles,
                               z_mean=z_mean,
                               V=V,
                               method=method)
    # mean free path length
    mfp = V / freq
    return mfp.to(u.m)


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n": {"units": u.m**-3}
                 })
def Spitzer_resistivity(T,
                        n,
                        particles,
                        z_mean=np.nan*u.dimensionless_unscaled,
                        V=np.nan*u.m/u.s,
                        method="classical"):
    r"""Spitzer resistivity of a plasma

    Parameters
    ----------

    T : Quantity
        Temperature in units of temperature.
        This should be the electron temperature for electron-electron
        and electron-ion collisions, and the ion temperature for
        ion-ion collisions.


    n : Quantity
        The density in units convertible to per cubic meter.
        This should be the electron density for electron-electron collisions,
        and the ion density for electron-ion and ion-ion collisions.

    z_mean : Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    V : Quantity, optional
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
        If the units on any of the inputs are incorrect

    UserWarning
        If the inputted velocity is greater than 80% of the speed of
        light.

    TypeError
        If the n_e, T, or V are not Quantities.

    Notes
    -----
    The Coulomb logarithm is given by

    .. math::
        \ln{\Lambda} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    where :math:`b_{min}` and :math:`b_{max}` are the inner and outer
    impact parameters for Coulomb collisions.

    See eq (2.14) in [1]_.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> Spitzer_resistivity(T, n, particles)
    <Quantity 2.4916169e-06 m Ohm>
    >>> Spitzer_resistivity(T, n, particles, V=1e6*u.m/u.s)
    <Quantity 2.4916169e-06 m Ohm>

    References
    ----------
    .. [1] http://homepages.cae.wisc.edu/~callen/chap2.pdf

    """
    # boiler plate checks
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    # collisional frequency
    freq = collision_frequency(T=T,
                               n=n,
                               particles=particles,
                               z_mean=z_mean,
                               V=V,
                               method=method)
    spitzer = freq * reduced_mass / (n * charges[0] * charges[1])
    return spitzer.to(u.Ohm * u.m)


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n_e": {"units": u.m**-3}
                 })
def mobility(T,
             n_e,
             particles,
             z_mean=np.nan*u.dimensionless_unscaled,
             V=np.nan*u.m/u.s,
             method="classical"):
    r"""Electrical mobility (m^2/(V s))

    Parameters
    ----------

    T : Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    n_e : Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    z_mean : Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters. It is also used the obtain the average mobility
        of a plasma with multiple charge state species. When z_mean
        is not given, the average charge between the two particles is
        used instead.

    V : Quantity, optional
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
        If the units on any of the inputs are incorrect

    UserWarning
        If the inputted velocity is greater than 80% of the speed of
        light.

    TypeError
        If the n_e, T, or V are not Quantities.

    Notes
    -----
    The Coulomb logarithm is given by

    .. math::
        \ln{\Lambda} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    where :math:`b_{min}` and :math:`b_{max}` are the inner and outer
    impact parameters for Coulomb collisions.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> mobility(T, n, particles)
    <Quantity 250500.35318738 m2 / (s V)>
    >>> mobility(T, n, particles, V=1e6*u.m/u.s)
    <Quantity 250500.35318738 m2 / (s V)>

    References
    ----------

    """
    # boiler plate checks
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    freq = collision_frequency(T=T,
                               n=n_e,
                               particles=particles,
                               z_mean=z_mean,
                               V=V,
                               method=method)
    if np.isnan(z_mean):
        z_val = (charges[0] + charges[1]) / 2
    else:
        z_val = z_mean
    mobility_value = z_val / (reduced_mass * freq)
    return mobility_value.to(u.m ** 2 / (u.V * u.s))


@check_quantity({"T": {"units": u.K, "can_be_negative": False},
                 "n_e": {"units": u.m**-3}
                 })
def knudsen(characteristic_length,
            T,
            n_e,
            particles,
            z_mean=np.nan*u.dimensionless_unscaled,
            V=np.nan*u.m/u.s,
            method="classical"):
    r"""Knudsen number (dimless)

    Parameters
    ----------

    T : Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    n_e : Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)

    z_mean : Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. This is used to recover the
        average ion density (given the average ionization and electron
        density) for calculating the ion sphere radius for non-classical
        impact parameters.

    V : Quantity, optional
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
        If the units on any of the inputs are incorrect

    UserWarning
        If the inputted velocity is greater than 80% of the speed of
        light.

    TypeError
        If the n_e, T, or V are not Quantities.

    Notes
    -----
    The Coulomb logarithm is given by

    .. math::
        \ln{\Lambda} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    where :math:`b_{min}` and :math:`b_{max}` are the inner and outer
    impact parameters for Coulomb collisions.

    Examples
    --------
    >>> from astropy import units as u
    >>> L = 1e-3 * u.m
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> knudsen(L, T, n, particles)
    <Quantity 7839.36310417>
    >>> knudsen(L, T, n, particles, V=1e6*u.m/u.s)
    <Quantity 1423.47708879>

    References
    ----------

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
                 "n_e": {"units": u.m**-3}
                 })
def coupling_parameter(T,
                       n_e,
                       particles,
                       V=np.nan*u.m/u.s,
                       method="classical"):
    r"""Coupling parameter.
    Coupling parameter compares Coulomb energy to kinetic energy (typically)
    thermal. Classical plasmas are weakly coupled Gamma << 1, whereas dense
    plasmas tend to have significant to strong coupling Gamma >= 1.

    Parameters
    ----------

    T : Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle

    n_e : Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second)
        
    V : Quantity, optional
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
        If the units on any of the inputs are incorrect

    UserWarning
        If the inputted velocity is greater than 80% of the speed of
        light.

    TypeError
        If the n_e, T, or V are not Quantities.

    Notes
    -----
    The Coulomb logarithm is given by

    .. math::
        \ln{\Lambda} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    where :math:`b_{min}` and :math:`b_{max}` are the inner and outer
    impact parameters for Coulomb collisions.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*u.m**-3
    >>> T = 1e6*u.K
    >>> particles = ('e', 'p')
    >>> coupling_parameter(T, n, particles)
    <Quantity 5.80330315e-05>
    >>> coupling_parameter(T, n, particles, V=1e6*u.m/u.s)
    <Quantity 5.80330315e-05>

    References
    ----------
    Dense plasma temperature equilibration in the binary collision
    approximation. D. O. Gericke et. al. PRE,  65, 036418 (2002).
    DOI: 10.1103/PhysRevE.65.036418


    """
    # boiler plate checks
    T, masses, charges, reduced_mass, V = _boilerPlate(T=T,
                                                       particles=particles,
                                                       V=V)
    # Wigner-Seitz radius
    radius = Wigner_Seitz_radius(n_e)
    # Coulomb potential energy between particles
    coulombEnergy = charges[0] * charges[1] / (4 * np.pi * eps0 * radius)
    if method == "classical":
        # classical thermal kinetic energy
        kineticEnergy = k_B * T
    elif method == "quantum":
        # quantum kinetic energy for dense plasmas
        lambda_deBroglie = thermal_deBroglie_wavelength(T)
        chemicalPotential = chemical_potential(n_e, T)
        fermiIntegral = Fermi_integral(chemicalPotential, 3 / 2)
        denom = (n_e * lambda_deBroglie ** 3) * fermiIntegral
        kineticEnergy = 2 * k_B * T / denom
    coupling = coulombEnergy / kineticEnergy
    return coupling.to(u.dimensionless_unscaled)


class classical_transport:
    r"""
    Braginskii-esque classical transport coefficients

    Notes
    -----
    Classical transport theory is derived by using kinetic theory to close the
    plasma two-fluid (electron and ion fluid) equations in the collisional
    limit. The first complete model in this form was done by S. I. Braginskii.

    This function uses fitting functions from literature to calculate
    the transport coefficients, which are the resistivity, thermoelectric
    conductivity, thermal conductivity, and viscosity.

    Note well the assumptions in the derivation of classical transport.
    It is assumed that the plasma is fully ionized, so only consisting of ions
    and electrons. Neutral atoms aren't considered. It's assumed that turbulent
    transport does not dominate, and the velocity distribution function is
    close to Maxwellian (no extremely strong gradients), which is equivalent
    to the following conditions:

    * collisional frequency >> gyrofrequency
    * collisional mean free path << gradient scale length along field
    * gyroradius << gradient scale length perpendicular to field

    When classical transport is not valid, e.g. due to the presence of strong
    gradients or turbulent transport, the transport is significantly increased
    by these other effects. Thus, classical transport often serves as a lower
    bound on the losses / transport encountered in a plasma.

    Models implemented:

    Braginskii [1]_
    The original Braginskii treatment as presented in the highly cited review
    paper from 1965. Coefficients are found from expansion of the kinetic
    equation in Laguerre / Sonine polynomials, truncated at K = 2. This theory
    allow for arbitrary Hall parameter and include results for Z = 1, 2, 3, 4,
    and infinity (Lorentz gas / stationary ion approximation).

    Spitzer-Harm [2]_ [3]_
    These coefficients were obtained from a numerical solution of the Fokker-
    Planck equation. They give one of the earliest and most accurate (in the
    Fokker-Planck sense) results for electron transport in simple plasma. They
    principally apply in the unmagnetized / parallel field case, although for
    resistivity Spitzer also calculated a famous result for a strong
    perpendicular magnetic field. Results are for Z = 1, 2, 4, 16,
    and infinity (Lorentz gas / stationary ion approximation).

    Epperlein-Haines [4]_
    Not yet implemented.

    Ji-Held [5]_
    This is a modern treatment of the classical transport problem that has been
    carried out with laudable care. It allows for arbitrary hall parameter and
    arbitrary Z for all coefficients. Similar to the Epperlein-Haines model,
    it corrects some known inaccuracies in the original Braginskii results,
    notably the asymptotic behavior of alpha-cross and beta_perp as Hall ->
    +infinity. It also studies effects of electron collisions in the ion
    terms, which all other treatments have not. To neglect electron-electron
    collisions, leave mu = 0. To consider them, specify mu and theta.

    Parameters
    ----------
    T_e : Quantity
        Temperature in units of temperature or energy per particle

    n_e : Quantity
        The number density in units convertible to per cubic meter.

    T_i : Quantity
        Temperature in units of temperature or energy per particle

    n_i : Quantity
        The number density in units convertible to per cubic meter.

    ion_particle : string
        Representation of the ion species (e.g., 'p' for protons,
        'e' for electrons, 'D+' for deuterium, or 'He-4 +1' for singly
        ionized helium-4). If no charge state information is provided,
        then the particles are assumed to be singly charged.

    Z : integer or np.inf, optional
        The ion charge state. Overrides particle charge state if included.
        Different theories support different values of Z. For the original
        Braginskii model, Z can be any of [1,2,3,4,infinity]. The Ji-Held
        model supports arbitrary Z. Average ionization states Z_mean can be
        input using this input and the Ji-Held model, although doing so may
        neglect effects caused by multiple ion populations.

    B : Quantity, optional
        The magnetic field strength in units convertible to Tesla. Defaults
        to zero.

    model: string
        Indication of whose formulation from literature to use. Allowed values
        are:
        'Braginskii',
        'Spitzer-Harm',
        'Epperlein-Haines' (not yet implemented),
        'Ji-Held'.
        See refs [1]_, [2]_, [3]_, [4]_ and [5]_.

    field_orientation : string
        Either of 'parallel', 'par', 'perpendicular', 'perp', 'cross', or
        'all', indicating the cardinal orientation of the magnetic field with
        respect to the transport direction of interest. Note that 'perp' refers
        to transport perpendicular to the field direction in the direction of
        the temperature gradient, while 'cross' refers to the B X grad(T)
        direction. The option 'all' will return a numpy array of all three,
        np.array((par, perp, cross)).

    coulomb_log_ei: float or dimensionless Quantity, optional
        Force a particular value to be used for the electron-ion Coulomb
        logarithm (test electrons on field ions). If None, the PlasmaPy
        function Coulomb_Logarithm() will be used. Useful for comparing
        calculations.

    V_ei: Quantity, optional
       Supplied to coulomb_logarithm() function, not otherwise used.
       The relative velocity between particles.  If not provided,
       thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
       where `mu` is the reduced mass.

    coulomb_log_ii: float or dimensionless Quantity, optional
        Force a particular value to be used for the ion-ion Coulomb logarithm
        (test ions on field ions). If None, the PlasmaPy function
        Coulomb_Logarithm() will be used. Useful for comparing calculations.

    V_ii: Quantity, optional
       Supplied to coulomb_logarithm() function, not otherwise used.
       The relative velocity between particles.  If not provided,
       thermal velocity is assumed: :math:`\mu V^2 \sim 2 k_B T`
       where `mu` is the reduced mass.

    hall_e: float or dimensionless Quantity, optional
        Force a particular value to be used for the electron Hall parameter. If
        None, the PlasmaPy function Hall_parameter() will be used. Useful
        for comparing calculations.

    hall_i: float or dimensionless Quantity, optional
        Force a particular value to be used for the ion Hall parameter. If
        None, the PlasmaPy function Hall_parameter() will be used. Useful
        for comparing calculations.

    mu: optional, float or dimensionless Quantity
        Ji-Held model only, may be used to include ion-electron effects
        on the ion transport coefficients. Defaults to zero, thus
        disabling these effects.

    theta: optional, float or dimensionless Quantity
        theta = T_e / T_i
        Ji-Held model only, may be used to include ion-electron effects
        on the ion transport coefficients. Defaults to T_e / T_i. Only
        has effect if mu is non-zero.


    Raises
    ------
    ValueError
        On incorrect or unknown values of arguments.
    plasmapy.utils.PhysicsError
        If input or calculated values for Coulomb logarithms are nonphysical.

    Examples
    --------
    >>> from astropy import units as u
    >>> t = classical_transport(1*u.eV, 1e20/u.m**3,
    ...                         1*u.eV, 1e20/u.m**3, 'p')
    >>> t.resistivity()
    <Quantity 0.00036701 m Ohm>
    >>> t.thermoelectric_conductivity()
    <Quantity 0.711084>
    >>> t.ion_thermal_conductivity()
    <Quantity 0.01552066 W / (K m)>
    >>> t.electron_thermal_conductivity()
    <Quantity 0.38064293 W / (K m)>
    >>> t.ion_viscosity()
    <Quantity [4.62129725e-07, 4.60724824e-07, 4.60724824e-07, 0.00000000e+00,
               0.00000000e+00] Pa s>
    >>> t.electron_viscosity()
    <Quantity [5.82273805e-09, 5.82082061e-09, 5.82082061e-09, 0.00000000e+00,
               0.00000000e+00] Pa s>

    References
    ----------
    .. [1] Braginskii, S. I. "Transport processes in a plasma." Reviews of
           plasma physics 1 (1965): 205. (1965)
    .. [2] Spitzer Jr, Lyman, and Richard Härm. "Transport phenomena in a
           completely ionized gas." Physical Review 89.5 (1953): 977. (1953)
    .. [3] Physics of Fully Ionized Gases, L. Spitzer (1962)
    .. [4] Epperlein, E. M., and M. G. Haines. "Plasma transport coefficients
           in a magnetic field by direct numerical solution of the
           Fokker–Planck equation." The Physics of fluids 29.4 (1986):
           1029-1041.
    .. [5] Ji, Jeong-Young, and Eric D. Held. "Closure and transport theory for
           high-collisionality electron-ion plasmas." Physics of Plasmas 20.4
           (2013): 042114.
    """

    @utils.check_quantity({"T_e": {"units": u.K, "can_be_negative": False},
                           "n_e": {"units": u.m**-3},
                           "T_i": {"units": u.K, "can_be_negative": False},
                           "n_i": {"units": u.m**-3},
                           })
    def __init__(self,
                 T_e,
                 n_e,
                 T_i,
                 n_i,
                 ion_particle,
                 m_i=None,
                 Z=None,
                 B=0.0 * u.T,
                 model='Braginskii',
                 field_orientation='parallel',
                 coulomb_log_ei=None,
                 V_ei=None,
                 coulomb_log_ii=None,
                 V_ii=None,
                 hall_e=None,
                 hall_i=None,
                 mu=None,
                 theta=None):
        # check the model
        self.model = model.lower()  # string inputs should be case insensitive
        valid_models = ['braginskii',
                        'spitzer',
                        'spitzer-harm',
                        'ji-held',
                        ]
        is_valid_model = self.model in valid_models
        if not is_valid_model:
            raise ValueError(f"Unknown transport model '{self.model}'")

        # check the field orientation
        self.field_orientation = field_orientation.lower()
        valid_fields = ['parallel', 'par',
                        'perpendicular', 'perp',
                        'cross',
                        'all',
                        ]
        is_valid_field = self.field_orientation in valid_fields
        if not is_valid_field:
            raise ValueError(f"Unknown field orientation "
                             f"'{self.field_orientation}'")

        # density and temperature units have already been checked by decorator
        # so just convert
        self.T_e = T_e.to(u.K, equivalencies=u.temperature_energy())
        self.T_i = T_i.to(u.K, equivalencies=u.temperature_energy())
        self.n_e = n_e.to(u.m**-3)
        self.n_i = n_i.to(u.m**-3)

        # get ion mass and charge state
        if m_i is None:
            try:
                self.m_i = atomic.ion_mass(ion_particle)
            except Exception:
                raise ValueError(f"Unable to find mass of particle: "
                                 f"{ion_particle} in classical_transport")
        else:
            self.m_i = m_i.to(u.kg)
        if Z is None:
            try:
                self.Z = atomic.integer_charge(ion_particle)
            except Exception:
                raise ValueError(f"Unable to find charge of particle: "
                                 f"{ion_particle} in classical_transport.")
        else:
            # red alert: the user has input a Z
            self.Z = Z
            # TODO: we need to make sure Lemmatum's z_mean is implemented here
            # if it's not available for a particular model, they'll complain
            # later
            if Z < 0:
                raise ValueError("Z not allowed to be negative")

        # decide on the particle string for the electrons
        self.e_particle = 'e'
        self.ion_particle = ion_particle

        # save other arguments
        self.B = B
        self.V_ei = V_ei
        self.V_ii = V_ii

        # calculate Coulomb logs if not forced in input
        if coulomb_log_ei:
            self.coulomb_log_ei = coulomb_log_ei
        else:
            self.coulomb_log_ei = Coulomb_logarithm(
                T_e, n_e, [self.e_particle, ion_particle], V_ei)
        if self.coulomb_log_ei < 4:
            warnings.warn(
                f"coulomb_log_ei is {self.coulomb_log_ei}, you might "
                "have strong coupling effects",
                utils.PhysicsWarning)
        if self.coulomb_log_ei < 1:
            raise utils.PhysicsError(
                f"coulomb_log_ei is {self.coulomb_log_ei}, less than 1")

        if coulomb_log_ii is not None:
            self.coulomb_log_ii = coulomb_log_ii
        else:
            # TODO make comment below more clear?
            self.coulomb_log_ii = Coulomb_logarithm(T_i,
                                                    n_e,  # this is not a typo!
                                                    [ion_particle,
                                                     ion_particle],
                                                    V_ii)

        if self.coulomb_log_ii < 4:
            warnings.warn(
                f"coulomb_log_ii is {self.coulomb_log_ii}, you might "
                "have strong coupling effects",
                utils.PhysicsWarning)
        if self.coulomb_log_ii < 1:
            raise utils.PhysicsError(
                f"coulomb_log_ii is {self.coulomb_log_ii}, less than 1")

        # calculate Hall parameters if not forced in input
        if hall_e is not None:
            self.hall_e = hall_e
        else:
            self.hall_e = Hall_parameter(
                n_e,
                T_e,
                B,
                self.e_particle,
                ion_particle,
                coulomb_log_ei,
                V_ei)
        if hall_i is not None:
            self.hall_i = hall_i
        else:
            self.hall_i = Hall_parameter(
                n_i, T_i, B, ion_particle, ion_particle, coulomb_log_ii, V_ii)
        # set up the ion non-dimensional coefficients for the Ji-Held model
        if mu is not None:
            self.mu = mu
        else:
            self.mu = 0  # disable the JH special features by default
#            self.mu = m_e / self.m_i  # enable the JH special features
        if theta is not None:
            self.theta = theta
        else:
            self.theta = self.T_e / self.T_i

    def resistivity(self):
        """
        Calculate the resistivity.

        Returns
        -------
        astropy.units.quantity.Quantity
        """
        alpha_hat = _nondim_resistivity(self.hall_e,
                                        self.Z,
                                        self.e_particle,
                                        self.model,
                                        self.field_orientation)
        tau_e = 1 / collision_rate_electron_ion(self.T_e,
                                                self.n_e,
                                                self.ion_particle,
                                                self.coulomb_log_ei,
                                                self.V_ei)

        alpha = alpha_hat / (self.n_e * e**2 * tau_e / m_e)
        return alpha.to(u.ohm * u.m)

    def thermoelectric_conductivity(self):
        """
        Calculate the thermoelectric conductivity.

        Returns
        -------
        astropy.units.quantity.Quantity
        """
        beta_hat = _nondim_te_conductivity(self.hall_e,
                                           self.Z,
                                           self.e_particle,
                                           self.model,
                                           self.field_orientation)
        beta = beta_hat * u.s / u.s  # yay! already dimensionless
        return beta

    def ion_thermal_conductivity(self):
        """
        Calculate the thermal conductivity for ions.

        Returns
        -------
        astropy.units.quantity.Quantity
        """
        kappa_hat = _nondim_thermal_conductivity(self.hall_i,
                                                 self.Z,
                                                 self.ion_particle,
                                                 self.model,
                                                 self.field_orientation,
                                                 self.mu,
                                                 self.theta)
        tau_i = 1 / collision_rate_ion_ion(self.T_i,
                                           self.n_i,
                                           self.ion_particle,
                                           self.coulomb_log_ii,
                                           self.V_ii)
        kappa = kappa_hat * (self.n_i * k_B**2 * self.T_i * tau_i / self.m_i)
        return kappa.to(u.W / u.m / u.K)

    def electron_thermal_conductivity(self):
        """
        Calculate the thermal conductivity for electrons.

        Returns
        -------
        astropy.units.quantity.Quantity
        """
        kappa_hat = _nondim_thermal_conductivity(self.hall_e,
                                                 self.Z,
                                                 self.e_particle,
                                                 self.model,
                                                 self.field_orientation,
                                                 self.mu,
                                                 self.theta)
        tau_e = 1 / collision_rate_electron_ion(self.T_e,
                                                self.n_e,
                                                self.ion_particle,
                                                self.coulomb_log_ei,
                                                self.V_ei)
        kappa = kappa_hat * (self.n_e * k_B**2 * self.T_e * tau_e / m_e)
        return kappa.to(u.W / u.m / u.K)

    def ion_viscosity(self):
        """
        Calculate the ion viscosity.

        Returns
        -------
        astropy.units.quantity.Quantity
        """
        eta_hat = _nondim_viscosity(self.hall_i,
                                    self.Z,
                                    self.ion_particle,
                                    self.model,
                                    self.field_orientation,
                                    self.mu,
                                    self.theta)
        tau_i = 1 / collision_rate_ion_ion(self.T_i,
                                           self.n_i,
                                           self.ion_particle,
                                           self.coulomb_log_ii,
                                           self.V_ii)
        common_factor = self.n_i * k_B * self.T_i * tau_i
        if np.isclose(self.hall_i, 0, rtol=1e-8):
            eta1 = (eta_hat[0] * common_factor,
                    eta_hat[1] * common_factor,
                    eta_hat[2] * common_factor,
                    eta_hat[3] * common_factor,
                    eta_hat[4] * common_factor)
        else:
            eta1 = (eta_hat[0] * common_factor,
                    eta_hat[1] * common_factor / self.hall_i**2,
                    eta_hat[2] * common_factor / self.hall_i**2,
                    eta_hat[3] * common_factor / self.hall_i,
                    eta_hat[4] * common_factor / self.hall_i)
        if eta1[0].unit == eta1[2].unit == eta1[4].unit:
            unit_val = eta1[0].unit
            eta = (np.array((eta1[0].value,
                             eta1[1].value,
                             eta1[2].value,
                             eta1[3].value,
                             eta1[4].value)) * unit_val).to(u.Pa * u.s)
        return eta

    def electron_viscosity(self):
        """
        Calculate the electron viscosity.

        Returns
        -------
        astropy.units.quantity.Quantity
        """
        eta_hat = _nondim_viscosity(self.hall_e,
                                    self.Z,
                                    self.e_particle,
                                    self.model,
                                    self.field_orientation,
                                    self.mu,
                                    self.theta)
        tau_e = 1 / collision_rate_electron_ion(self.T_e,
                                                self.n_e,
                                                self.ion_particle,
                                                self.coulomb_log_ei,
                                                self.V_ei)
        common_factor = (self.n_e * k_B * self.T_e * tau_e)
        if np.isclose(self.hall_e, 0, rtol=1e-8):
            eta1 = (eta_hat[0] * common_factor,
                    eta_hat[1] * common_factor,
                    eta_hat[2] * common_factor,
                    eta_hat[3] * common_factor,
                    eta_hat[4] * common_factor)
        else:
            eta1 = (eta_hat[0] * common_factor,
                    eta_hat[1] * common_factor / self.hall_e**2,
                    eta_hat[2] * common_factor / self.hall_e**2,
                    eta_hat[3] * common_factor / self.hall_e,
                    eta_hat[4] * common_factor / self.hall_e)
        if eta1[0].unit == eta1[2].unit and eta1[2].unit == eta1[4].unit:
            unit_val = eta1[0].unit
            eta = (np.array((eta1[0].value,
                             eta1[1].value,
                             eta1[2].value,
                             eta1[3].value,
                             eta1[4].value)) * unit_val).to(u.Pa * u.s)
        return eta

    def all_variables(self) -> dict:
        """
        Return all transport variables as a dictionary.

        Returns
        -------
        dict
        """
        d = {'resistivity': self.resistivity(),
             'thermoelectric_conductivity': self.thermoelectric_conductivity(),
             'electron_thermal_conductivity':
                 self.electron_thermal_conductivity(),
             'electron_viscosity': self.electron_viscosity()}
        if self.model != "spitzer":
            d = dict(d, **{'ion_thermal_conductivity':
                           self.ion_thermal_conductivity(),
                           'ion_viscosity': self.ion_viscosity()})
        return d


def _nondim_thermal_conductivity(hall, Z,
                                 particle,
                                 model,
                                 field_orientation,
                                 mu=None,
                                 theta=None):
    """calculate dimensionless classical thermal conductivity coefficients

    This function is a switchboard / wrapper that calls the appropriate
    model-specific functions depending on which model is specified and which
    type of particle (electron or ion) is input. Non-electrons are assumed to
    be ions.
    """
    if _is_electron(particle):
        if model == 'spitzer-harm' or model == 'spitzer':
            kappa_hat = _nondim_tc_e_spitzer(Z)
        elif model == 'braginskii':
            kappa_hat = _nondim_tc_e_braginskii(hall, Z, field_orientation)
        elif model == 'ji-held':
            kappa_hat = _nondim_tc_e_ji_held(hall, Z, field_orientation)
        else:
            raise ValueError(f"Unrecognized model '{model}' in "
                             "_nondim_thermal_conductivity")
    else:
        if model == 'braginskii':
            kappa_hat = _nondim_tc_i_braginskii(hall, field_orientation)
        elif model == 'ji-held':
            kappa_hat = _nondim_tc_i_ji_held(hall, Z, mu, theta,
                                             field_orientation)
        elif model == 'spitzer-harm' or model == 'spitzer':
            raise NotImplementedError("Ion thermal conductivity is not "
                                      "implemented in the Spitzer model.")
        else:
            raise ValueError(f"Unrecognized model '{model}' in "
                             "_nondim_thermal_conductivity")
    return kappa_hat


def _nondim_viscosity(hall,
                      Z,
                      particle,
                      model,
                      field_orientation,
                      mu=None,
                      theta=None):
    """calculate dimensionless classical viscosity coefficients

    This function is a switchboard / wrapper that calls the appropriate
    model-specific functions depending on which model is specified and which
    type of particle (electron or ion) is input. Non-electrons are assumed to
    be ions.
    """

    if _is_electron(particle):
        if model == 'braginskii':
            eta_hat = _nondim_visc_e_braginskii(hall, Z)
        elif model == 'ji-held':
            eta_hat = _nondim_visc_e_ji_held(hall, Z)
        else:
            raise ValueError(f"Unrecognized model '{model}' in "
                             "_nondim_viscosity")
    else:
        if model == 'braginskii':
            eta_hat = _nondim_visc_i_braginskii(hall)
        elif model == 'ji-held':
            eta_hat = _nondim_visc_i_ji_held(hall, Z, mu, theta)
        elif model == 'spitzer-harm' or model == 'spitzer':
            raise NotImplementedError("Ion viscosity is not "
                                      "implemented in the Spitzer model.")
        else:
            raise ValueError(f"Unrecognized model '{model}' in "
                             "_nondim_viscosity")
    return eta_hat


def _nondim_resistivity(hall, Z, particle, model, field_orientation):
    """calculate dimensionless classical resistivity coefficients

    This function is a switchboard / wrapper that calls the appropriate
    model-specific functions depending on which model is specified.
    """

    if model == 'spitzer-harm' or model == 'spitzer':
        alpha_hat = _nondim_resist_spitzer(Z, field_orientation)
    elif model == 'braginskii':
        alpha_hat = _nondim_resist_braginskii(hall, Z, field_orientation)
    elif model == 'ji-held':
        alpha_hat = _nondim_resist_ji_held(hall, Z, field_orientation)
    else:
        raise ValueError(f"Unrecognized model '{model}' in "
                         "_nondim_resistivity")
    return alpha_hat


def _nondim_te_conductivity(hall, Z, particle, model, field_orientation):
    """calculate dimensionless classical thermoelectric coefficients

    This function is a switchboard / wrapper that calls the appropriate
    model-specific functions depending on which model is specified.
    """

    if model == 'spitzer-harm' or model == 'spitzer':
        beta_hat = _nondim_tec_spitzer(Z)
    elif model == 'braginskii':
        beta_hat = _nondim_tec_braginskii(hall, Z, field_orientation)
    elif model == 'ji-held':
        beta_hat = _nondim_tec_ji_held(hall, Z, field_orientation)
    else:
        raise ValueError(f"Unrecognized model '{model}' in "
                         "_nondim_te_conductivity")
    return beta_hat


def _check_Z(allowed_Z, Z):
    """determine if the input Z value is okay given the list of allowed_Z"""
    # first, determine if arbitrary Z values are allowed in the theory
    arbitrary_Z_allowed = False
    the_arbitrary_idx = np.nan
    for idx, allowed_Z_val in enumerate(allowed_Z):
        if allowed_Z_val == 'arbitrary':
            arbitrary_Z_allowed = True
            the_arbitrary_idx = idx
    # next, search the allowed_Z for a match to the current Z
    Z_idx = np.nan
    for idx, allowed_Z_val in enumerate(allowed_Z):
        if Z == allowed_Z_val:
            Z_idx = idx
    # at this point we have looped through allowed_Z and either found a match
    # or not. If we haven't found a match and arbitrary Z aren't allowed, break
    if np.isnan(Z_idx) and not arbitrary_Z_allowed:
        raise utils.PhysicsError(f"{Z} is not an allowed Z value")
    elif np.isnan(Z_idx):  # allowed arbitrary Z
        # return a Z_idx pointing to the 'arbitrary'
        Z_idx = the_arbitrary_idx
    else:  # allowed Z
        pass
    # we have got the Z_idx we want. return
    return Z_idx


def _get_spitzer_harm_coeffs(Z):
    """return numerical coefficients from Spitzer-Harm '53

    Table III, Spitzer and Harm, Phys. Rev. Vol 89, 5, 1953
    """
    allowed_Z = [1, 2, 4, 16, np.inf]
    Z_idx = _check_Z(allowed_Z, Z)
    gamma_E = [0.5816, 0.6833, 0.7849, 0.9225, 1.0000]
    gamma_T = [0.2727, 0.4137, 0.5714, 0.8279, 1.0000]
    delta_E = [0.4652, 0.5787, 0.7043, 0.8870, 1.0000]
    delta_T = [0.2252, 0.3563, 0.5133, 0.7907, 1.0000]
    return gamma_E[Z_idx], gamma_T[Z_idx], delta_E[Z_idx], delta_T[Z_idx]


def _nondim_tc_e_spitzer(Z):
    """dimensionless electron thermal conductivity - Spitzer

    This result is for parallel field or unmagnetized plasma only.
    """
    (gamma_E, gamma_T, delta_E, delta_T) = _get_spitzer_harm_coeffs(Z)
    kappa = (64 / np.pi) * delta_T * \
            (5 / 3 - (gamma_T * delta_E) / (delta_T * gamma_E))
    return kappa


def _nondim_resist_spitzer(Z, field_orientation):
    """dimensionless resistivity - Spitzer

    These are results for both parallel-field / unmagnetized plasmas as well
    as perpendicular-field / strongly magnetized plasma. Summary description
    in Physics of Fully Ionized Gases, Spitzer"""

    alpha_perp = 1
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return alpha_perp

    (gamma_E, gamma_T, delta_E, delta_T) = _get_spitzer_harm_coeffs(Z)
    alpha_par = (3 * np.pi / 32) * (1 / gamma_E)
    if field_orientation == 'parallel' or field_orientation == 'par':
        return alpha_par
#        alpha_par = 0.5064 # Z = 1

    if field_orientation == 'all':
        return alpha_par, alpha_perp


def _nondim_tec_spitzer(Z):
    """dimensionless thermoelectric conductivity - Spitzer

    This result is for parallel field or unmagnetized plasma only.
    """
    (gamma_E, gamma_T, delta_E, delta_T) = _get_spitzer_harm_coeffs(Z)
    beta = 5 / 2 * (8 / 5 * (delta_E / gamma_E) - 1)
#    beta = 0.703
    return beta


def _nondim_tc_e_braginskii(hall, Z, field_orientation):
    """dimensionless electron thermal conductivity - Braginskii

    Braginskii, S. I. "Transport processes in a plasma." Reviews of plasma
    physics 1 (1965): 205.
    """

    allowed_Z = [1, 2, 3, 4, np.inf]
    Z_idx = _check_Z(allowed_Z, Z)

    delta_0 = [3.7703, 1.0465, 0.5814, 0.4106, 0.0961]
    delta_1 = [14.79, 10.80, 9.618, 9.055, 7.482]
    gamma_1_prime = [4.664, 3.957, 3.721, 3.604, 3.25]
    gamma_0_prime = [11.92, 5.118, 3.525, 2.841, 1.20]
    gamma_1_doubleprime = [2.500, 2.500, 2.500, 2.500, 2.500]
    gamma_0_doubleprime = [21.67, 15.37, 13.53, 12.65, 10.23]

    gamma_0 = gamma_0_prime[Z_idx] / delta_0[Z_idx]
    Delta = hall ** 4 + delta_1[Z_idx] * hall ** 2 + delta_0[Z_idx]
    kappa_par = gamma_0
    if field_orientation == 'parallel' or field_orientation == 'par':
        return kappa_par

    kappa_perp = (gamma_1_prime[Z_idx] * hall **
                  2 + gamma_0_prime[Z_idx]) / Delta
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return kappa_perp

    kappa_cross = (gamma_1_doubleprime[Z_idx] * hall ** 3 +
                   gamma_0_doubleprime[Z_idx] * hall) / Delta
    if field_orientation == 'cross':
        return kappa_cross

    if field_orientation == 'all':
        return np.array((kappa_par, kappa_perp, kappa_cross))


def _nondim_tc_i_braginskii(hall, field_orientation):
    """dimensionless ion thermal conductivity - Braginskii

    Braginskii, S. I. "Transport processes in a plasma." Reviews of plasma
    physics 1 (1965): 205.
    """

    kappa_par_coeff_0 = 3.906
    kappa_par = kappa_par_coeff_0
    if field_orientation == 'parallel' or field_orientation == 'par':
        return kappa_par

    kappa_perp_coeff_2 = 2.0
    kappa_perp_coeff_0 = 2.645
    delta_1 = 2.70
    delta_0 = 0.677
    Delta = hall ** 4 + delta_1 * hall ** 2 + delta_0
    kappa_perp = (kappa_perp_coeff_2 * hall ** 2 +
                  kappa_perp_coeff_0) / Delta
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return kappa_perp

    kappa_cross_coeff_3 = 2.5
    kappa_cross_coeff_1 = 4.65
    kappa_cross = (kappa_cross_coeff_3 * hall ** 3 +
                   kappa_cross_coeff_1 * hall) / Delta
    if field_orientation == 'cross':
        return kappa_cross

    if field_orientation == 'all':
        return np.array((kappa_par, kappa_perp, kappa_cross))


def _nondim_visc_e_braginskii(hall, Z):
    """dimensionless electron viscosity - Braginskii

    Braginskii, S. I. "Transport processes in a plasma." Reviews of plasma
    physics 1 (1965): 205.
    """
    allowed_Z = [1]
    _check_Z(allowed_Z, Z)
    eta_prime_0 = 0.733
    eta_doubleprime_2 = 2.05
    eta_doubleprime_0 = 8.50
    eta_tripleprime_2 = 1.0
    eta_tripleprime_0 = 7.91
    delta_1 = 13.8
    delta_0 = 11.6
    eta_0_e = eta_prime_0

    def f_eta_2(hall):
        Delta = hall ** 4 + delta_1 * hall ** 2 + delta_0
        return (eta_doubleprime_2 * hall ** 2 + eta_doubleprime_0) / Delta
    eta_2_e = f_eta_2(hall)
    eta_1_e = f_eta_2(2 * hall)

    def f_eta_4(hall):
        Delta = hall ** 4 + delta_1 * hall ** 2 + delta_0
        return (eta_tripleprime_2 * hall ** 3 +
                eta_tripleprime_0 * hall) / Delta
    eta_4_e = f_eta_4(hall)
    eta_3_e = f_eta_4(2 * hall)
    return np.array((eta_0_e, eta_1_e, eta_2_e, eta_3_e, eta_4_e))


def _nondim_visc_i_braginskii(hall):
    """dimensionless ion viscosity - Braginskii

    Braginskii, S. I. "Transport processes in a plasma." Reviews of plasma
    physics 1 (1965): 205.
    """
    eta_prime_0 = 0.96
    eta_doubleprime_2 = 6 / 5
    eta_doubleprime_0 = 2.23
    eta_tripleprime_2 = 1.0
    eta_tripleprime_0 = 2.38
    delta_1 = 4.03
    delta_0 = 2.33
    eta_0_i = eta_prime_0

    def f_eta_2(hall):
        Delta = hall ** 4 + delta_1 * hall ** 2 + delta_0
        return (eta_doubleprime_2 * hall ** 2 + eta_doubleprime_0) / Delta
    eta_2_i = f_eta_2(hall)
    eta_1_i = f_eta_2(2 * hall)

    def f_eta_4(hall):
        Delta = hall**4 + delta_1 * hall ** 2 + delta_0
        return (eta_tripleprime_2 * hall ** 3 +
                eta_tripleprime_0 * hall) / Delta
    eta_4_i = f_eta_4(hall)
    eta_3_i = f_eta_4(2 * hall)
    return np.array((eta_0_i, eta_1_i, eta_2_i, eta_3_i, eta_4_i))


def _nondim_resist_braginskii(hall, Z, field_orientation):
    """dimensionless resistivity - Braginskii

    Braginskii, S. I. "Transport processes in a plasma." Reviews of plasma
    physics 1 (1965): 205.
    """

    allowed_Z = [1, 2, 3, 4, np.inf]
    Z_idx = _check_Z(allowed_Z, Z)

#    alpha_0 = 0.5129
    delta_0 = [3.7703, 1.0465, 0.5814, 0.4106, 0.0961]
    delta_1 = [14.79, 10.80, 9.618, 9.055, 7.482]
    alpha_1_prime = [6.416, 5.523, 5.226, 5.077, 4.63]
    alpha_0_prime = [1.837, 0.5956, 0.3515, 0.2566, 0.0678]
    alpha_1_doubleprime = [1.704, 1.704, 1.704, 1.704, 1.704]
    alpha_0_doubleprime = [0.7796, 0.3439, 0.2400, 0.1957, 0.0940]

    alpha_0 = 1 - alpha_0_prime[Z_idx] / delta_0[Z_idx]
    Delta = hall ** 4 + delta_1[Z_idx] * hall ** 2 + delta_0[Z_idx]
    alpha_par = alpha_0
    if field_orientation == 'parallel' or field_orientation == 'par':
        return alpha_par

    alpha_perp = (1 - (alpha_1_prime[Z_idx] * hall ** 2 +
                       alpha_0_prime[Z_idx]) / Delta)
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return alpha_perp

    alpha_cross = (alpha_1_doubleprime[Z_idx] * hall ** 3 +
                   alpha_0_doubleprime[Z_idx] * hall) / Delta
    if field_orientation == 'cross':
        return alpha_cross

    if field_orientation == 'all':
        return np.array((alpha_par, alpha_perp, alpha_cross))


def _nondim_tec_braginskii(hall, Z, field_orientation):
    """dimensionless thermoelectric conductivity - Braginskii

    Braginskii, S. I. "Transport processes in a plasma." Reviews of plasma
    physics 1 (1965): 205.
    """

    allowed_Z = [1, 2, 3, 4, np.inf]
    Z_idx = _check_Z(allowed_Z, Z)

    delta_0 = [3.7703, 1.0465, 0.5814, 0.4106, 0.0961]
    delta_1 = [14.79, 10.80, 9.618, 9.055, 7.482]
    beta_1_prime = [5.101, 4.450, 4.233, 4.124, 3.798]
    beta_0_prime = [2.681, 0.9473, 0.5905, 0.4478, 0.1461]
    beta_1_doubleprime = [1.5, 1.5, 1.5, 1.5, 1.5]
    beta_0_doubleprime = [3.053, 1.784, 1.442, 1.285, 0.877]

    Delta = hall**4 + delta_1[Z_idx] * hall**2 + delta_0[Z_idx]
    beta_0 = beta_0_prime[Z_idx] / delta_0[Z_idx]
#    beta_0 = 0.7110

    beta_par = beta_0
    if field_orientation == 'parallel' or field_orientation == 'par':
        return beta_par

    beta_perp = (beta_1_prime[Z_idx] * hall ** 2 +
                 beta_0_prime[Z_idx]) / Delta
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return beta_perp

    beta_cross = (beta_1_doubleprime[Z_idx] * hall ** 3 +
                  beta_0_doubleprime[Z_idx] * hall) / Delta
    if field_orientation == 'cross':
        return beta_cross

    if field_orientation == 'all':
        return np.array((beta_par, beta_perp, beta_cross))


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#               Abandon all hope, ye who enter here
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#


def _nondim_tc_e_ji_held(hall, Z, field_orientation):
    """dimensionless electron thermal conductivity - Ji-Held

    Ji, Jeong-Young, and Eric D. Held. "Closure and transport theory for
    high-collisionality electron-ion plasmas." Physics of Plasmas 20.4 (2013):
    042114.
    """

    allowed_Z = [1, 2, 'arbitrary']
    Z_idx = _check_Z(allowed_Z, Z)
    r = np.abs(Z * hall)

    def f_kappa_par_e(Z):
        numerator = 13.5 * Z ** 2 + 54.4 * Z + 25.2
        denominator = Z ** 3 + 8.35 * Z ** 2 + 15.2 * Z + 4.51
        return numerator / denominator

    def f_kappa_0(Z):
        numerator = 9.91 * Z ** 3 + 75.3 * Z ** 2 + 518 * Z + 333
        denominator = 1000
        return numerator / denominator

    def f_kappa_1(Z):
        numerator = 0.211 * Z ** 3 + 12.7 * Z ** 2 + 48.4 * Z + 6.45
        denominator = Z + 57.1
        return numerator / denominator

    def f_kappa_2(Z):
        numerator = 0.932 * Z ** (7 / 3) + 0.135 * Z ** 2 + 12.3 * Z + 8.77
        denominator = Z + 4.84
        return numerator / denominator

    def f_kappa_3(Z):
        numerator = 0.246 * Z ** 3 + 2.65 * Z ** 2 - 92.8 * Z - 1.96
        denominator = Z ** 2 + 19.9 * Z + 35.3
        return numerator / denominator

    def f_kappa_4(Z):
        numerator = 2.76 * Z ** (5 / 3) - 0.836 * Z ** (2 / 3) - 0.0611
        denominator = Z - 0.214
        return numerator / denominator

    def f_k_0(Z):
        numerator = 0.0396 * Z ** 3 + 46.3 * Z + 176
        denominator = 1000
        return numerator / denominator

    def f_k_1(Z):
        numerator = 15.4 * Z ** 3 + 188 * Z ** 2 + 240 * Z + 35.3
        denominator = 1000 * Z + 397
        return numerator / denominator

    def f_k_2(Z):
        numerator = -0.159 * Z ** 2 - 12.5 * Z + 34.1
        denominator = Z ** (2 / 3) + 0.741 * Z ** (1 / 3) + 31.0
        return numerator / denominator

    def f_k_3(Z):
        numerator = 0.431 * Z ** 2 + 3.69 * Z + 0.0314
        denominator = Z + 3.62
        return numerator / denominator

    def f_k_4(Z):
        numerator = 0.0258 * Z ** 2 - 1.63 * Z + 0.711
        denominator = Z ** (4 / 3) + 4.36 * Z ** (2 / 3) + 2.75
        return numerator / denominator

    def f_k_5(Z):
        numerator = Z ** 3 + 11.9 * Z ** 2 + 28.8 * Z + 9.07
        denominator = 173 * Z + 133
        return numerator / denominator

    kappa_par_e = [3.204, 2.464, f_kappa_par_e(Z)]
    kappa_0 = [0.936, 1.749, f_kappa_0(Z)]
    kappa_1 = [1.166, 2.635, f_kappa_1(Z)]
    kappa_2 = [3.791, 5.644, f_kappa_2(Z)]
    kappa_3 = [-1.635, -2.212, f_kappa_3(Z)]
    kappa_4 = [2.370, 4.129, f_kappa_4(Z)]
    k_0 = [0.222, 0.269, f_k_0(Z)]
    k_1 = [0.343, 0.580, f_k_1(Z)]
    k_2 = [0.655, 0.252, f_k_2(Z)]
    k_3 = [0.899, 1.626, f_k_3(Z)]
    k_4 = [-0.110, -0.201, f_k_4(Z)]
    k_5 = [0.166, 0.255, f_k_5(Z)]

    kappa_par = kappa_par_e[Z_idx]
    if field_orientation == 'parallel' or field_orientation == 'par':
        return Z * kappa_par

    def f_kappa_perp(Z_idx):
        numerator = ((13 / 4 * Z + np.sqrt(2)) * r +
                     kappa_0[Z_idx] * kappa_par_e[Z_idx])
        denominator = (r ** 3 +
                       kappa_4[Z_idx] * r ** (7 / 3) +
                       kappa_3[Z_idx] * r ** 2 +
                       kappa_2[Z_idx] * r ** (5 / 3) +
                       kappa_1[Z_idx] * r +
                       kappa_0[Z_idx])
        return numerator / denominator

    kappa_perp = f_kappa_perp(Z_idx)
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return Z * kappa_perp

    def f_kappa_cross(Z_idx):
        numerator = r * (5 / 2 * r + k_0[Z_idx] / k_5[Z_idx])
        denominator = (r ** 3 +
                       k_4[Z_idx] * r ** (7 / 3) +
                       k_3[Z_idx] * r ** 2 +
                       k_2[Z_idx] * r ** (5 / 3) +
                       k_1[Z_idx] * r +
                       k_0[Z_idx])
        return numerator / denominator

    kappa_cross = f_kappa_cross(Z_idx)
    if field_orientation == 'cross':
        return Z * kappa_cross

    if field_orientation == 'all':
        return np.array((Z * kappa_par, Z * kappa_perp, Z * kappa_cross))


def _nondim_resist_ji_held(hall, Z, field_orientation):
    """dimensionless resistivity - Ji-Held

    Ji, Jeong-Young, and Eric D. Held. "Closure and transport theory for
    high-collisionality electron-ion plasmas." Physics of Plasmas 20.4 (2013):
    042114.
    """

    allowed_Z = [1, 2, 'arbitrary']
    Z_idx = _check_Z(allowed_Z, Z)
    r = np.abs(Z * hall)

    def f_alpha_par_e(Z):
        numerator = Z ** (2 / 3)
        denominator = 1.46 * Z ** (2 / 3) - 0.330 * Z ** (1 / 3) + 0.888
        return 1 - numerator / denominator

    def f_alpha_0(Z):
        return 0.623 * Z ** (5 / 3) - 2.61 * Z ** (4 / 3) + 3.56 * Z + 0.557

    def f_alpha_1(Z):
        return 2.24 * Z ** (2 / 3) - 1.11 * Z ** (1 / 3) + 1.84

    def f_alpha_2(Z):
        return -0.0983 * Z ** (1 / 3) + 0.0176

    def f_a_0(Z):
        return 0.0759 * Z **( 8 / 3) + 0.897 * Z ** 2 + 2.06 * Z + 1.06

    def f_a_1(Z):
        return 2.18 * Z ** (5 / 3) + 5.31 * Z + 3.73

    def f_a_2(Z):
        return 7.41 * Z + 1.11 * Z ** (2 / 3) - 1.17

    def f_a_3(Z):
        return 3.89 * Z ** (2 / 3) - 4.51 * Z ** (1 / 3) + 6.76

    def f_a_4(Z):
        return 2.26 * Z ** (1 / 3) + 0.281

    def f_a_5(Z):
        return 1.18 * Z ** (5 / 3) - 1.03 * Z ** (4 / 3) + 3.60 * Z + 1.32

    alpha_par_e = [0.504, 0.431, f_alpha_par_e(Z)]
    alpha_0 = [2.130, 3.078, f_alpha_0(Z)]
    alpha_1 = [2.970, 3.997, f_alpha_1(Z)]
    alpha_2 = [-0.081, -0.106, f_alpha_2(Z)]
    a_0 = [4.093, 9.250, f_a_0(Z)]
    a_1 = [11.22, 21.27, f_a_1(Z)]
    a_2 = [7.350, 15.41, f_a_2(Z)]
    a_3 = [6.140, 7.253, f_a_3(Z)]
    a_4 = [2.541, 3.128, f_a_4(Z)]
    a_5 = [5.070, 9.671, f_a_5(Z)]

    alpha_par = alpha_par_e[Z_idx]
    if field_orientation == 'parallel' or field_orientation == 'par':
        return alpha_par

    def f_alpha_perp(Z_idx):
        numerator = (1.46 * Z**(2 / 3) * r +
                     alpha_0[Z_idx] * (1 - alpha_par_e[Z_idx]))
        denominator = (r ** (5 / 3) +
                       alpha_2[Z_idx] * r ** (4 / 3) +
                       alpha_1[Z_idx] * r +
                       alpha_0[Z_idx])
        return 1 - numerator / denominator

    alpha_perp = f_alpha_perp(Z_idx)
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return alpha_perp

    def f_alpha_cross(Z_idx):
        numerator = Z ** (2 / 3) * r * (2.53 * r + a_0[Z_idx] / a_5[Z_idx])
        denominator = (r ** (8 / 3) +
                       a_4[Z_idx] * r ** (7 / 3) +
                       a_3[Z_idx] * r ** 2 +
                       a_2[Z_idx] * r ** (5 / 3) +
                       a_1[Z_idx] * r +
                       a_0[Z_idx])
        return numerator / denominator

    alpha_cross = f_alpha_cross(Z_idx)
    if field_orientation == 'cross':
        return alpha_cross

    if field_orientation == 'all':
        return np.array((alpha_par, alpha_perp, alpha_cross))


def _nondim_tec_ji_held(hall, Z, field_orientation):
    """dimensionless thermoelectric conductivity - Ji-Held

    Ji, Jeong-Young, and Eric D. Held. "Closure and transport theory for
    high-collisionality electron-ion plasmas." Physics of Plasmas 20.4 (2013):
    042114.
    """

    allowed_Z = [1, 2, 'arbitrary']
    Z_idx = _check_Z(allowed_Z, Z)
    r = np.abs(Z * hall)

    def f_beta_par_e(Z):
        numerator = Z ** (5 / 3)
        denominator = 0.693 * Z ** (5 / 3) - 0.279 * Z ** (4 / 3) + Z + 0.01
        return numerator / denominator

    def f_beta_0(Z):
        return 0.156 * Z ** (8 / 3) + 0.994 * Z ** 2 + 3.21 * Z - 0.84

    def f_beta_1(Z):
        return 3.69 * Z ** (5 / 3) + 3.77 * Z + 0.77

    def f_beta_2(Z):
        return 9.43 * Z + 4.22 * Z ** (2 / 3) - 12.9 * Z ** (1 / 3) + 4.56

    def f_beta_3(Z):
        return 2.70 * Z ** (2 / 3) + 1.46 * Z ** (1 / 3) - 0.17

    def f_beta_4(Z):
        return 2.58 * Z ** (1 / 3) + 0.17

    def f_b_0(Z):
        numerator = 6.87 * Z ** 3 + 78.2 * Z ** 2 + 623 * Z + 366
        denominator = 1000
        return numerator / denominator

    def f_b_1(Z):
        return 0.134 * Z**2 + 0.977 * Z + 0.17

    def f_b_2(Z):
        return (0.689 * Z ** (4 / 3) - 0.377 * Z ** (2 / 3) +
                3.94 * Z**(1 / 3) + 0.644)

    def f_b_3(Z):
        return -0.109 * Z + 1.33 * Z ** (2 / 3) - 3.80 * Z ** (1 / 3) + 0.289

    def f_b_4(Z):
        return 2.46 * Z ** (2 / 3) + 0.522

    def f_b_5(Z):
        return 0.102 * Z ** 2 + 0.746 * Z + 0.072 * Z ** (1 / 3) + 0.211

    beta_par_e = [0.702, 0.905, f_beta_par_e(Z)]
    beta_0 = [3.520, 10.55, f_beta_0(Z)]
    beta_1 = [8.230, 20.03, f_beta_1(Z)]
    beta_2 = [5.310, 13.87, f_beta_2(Z)]
    beta_3 = [3.990, 5.955, f_beta_3(Z)]
    beta_4 = [2.750, 3.421, f_beta_4(Z)]
    b_0 = [1.074, 1.980, f_b_0(Z)]
    b_1 = [1.281, 2.660, f_b_1(Z)]
    b_2 = [4.896, 6.746, f_b_2(Z)]
    b_3 = [-2.290, -2.605, f_b_3(Z)]
    b_4 = [2.982, 4.427, f_b_4(Z)]
    b_5 = [1.131, 2.202, f_b_5(Z)]

    beta_par = beta_par_e[Z_idx]
    if field_orientation == 'parallel' or field_orientation == 'par':
        return beta_par

    def f_beta_perp(Z_idx):
        numerator = 6.33 * Z ** (5 / 3) * r + beta_0[Z_idx] * beta_par_e[Z_idx]
        denominator = (r ** (8 / 3) +
                       beta_4[Z_idx] * r ** (7 / 3) +
                       beta_3[Z_idx] * r ** 2 +
                       beta_2[Z_idx] * r ** (5 / 3) +
                       beta_1[Z_idx] * r +
                       beta_0[Z_idx])
        return numerator / denominator

    beta_perp = f_beta_perp(Z_idx)
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return beta_perp

    def f_beta_cross(Z_idx):
        numerator = Z * r * (3 / 2 * r + b_0[Z_idx] / b_5[Z_idx])
        denominator = (r ** 3 +
                       b_4[Z_idx] * r ** (7 / 3) +
                       b_3[Z_idx] * r ** 2 +
                       b_2[Z_idx] * r ** (5 / 3) +
                       b_1[Z_idx] * r +
                       b_0[Z_idx])
        return numerator / denominator

    beta_cross = f_beta_cross(Z_idx)
    if field_orientation == 'cross':
        return beta_cross

    if field_orientation == 'all':
        return np.array((beta_par, beta_perp, beta_cross))


def _nondim_visc_e_ji_held(hall, Z):
    """dimensionless electron viscosity - Ji-Held

    Ji, Jeong-Young, and Eric D. Held. "Closure and transport theory for
    high-collisionality electron-ion plasmas." Physics of Plasmas 20.4 (2013):
    042114.
    """

    allowed_Z = [1, 2, 'arbitrary']
    Z_idx = _check_Z(allowed_Z, Z)
    r = np.abs(Z * hall)

    def f_eta_0_e(Z):
        return 1 / (0.55 * Z + 0.083 * Z ** (1 / 3) + 0.732)

    def f_hprime_0(Z):
        return 0.0699 * Z ** 3 + 0.558 * Z ** 2 + 1.66 * Z + 1.06

    def f_hprime_1(Z):
        return 0.657 * Z ** 2 + 1.42 * Z + 0.416

    def f_hprime_2(Z):
        return -0.369 * Z ** (4 / 3) + 0.379 * Z + 0.339 * Z ** (1 / 3) + 2.17

    def f_hprime_3(Z):
        return 2.16 * Z - 0.657 * Z ** (1 / 3) + 0.0347

    def f_hprime_4(Z):
        return -0.0703 * Z ** (2 / 3) - 0.224 * Z ** (1 / 3) + 0.333

    def f_h_0(Z):
        return 0.0473 * Z ** 3 + 0.323 * Z ** 2 + 0.951 * Z + 0.407

    def f_h_1(Z):
        return 0.171 * Z ** 2 + 0.523 * Z + 0.336

    def f_h_2(Z):
        return 0.362 * Z ** (4 / 3) + 0.178 * Z + 1.06 * Z ** (1 / 3) + 1.26

    def f_h_3(Z):
        return 0.599 * Z + 0.106 * Z ** (2 / 3) - 0.444 * Z ** (1 / 3) - 0.161

    def f_h_4(Z):
        return -0.16 * Z ** (2 / 3) + 0.06 * Z ** (1 / 3) + 0.232

    def f_h_5(Z):
        return 0.183 * Z ** 2 + 0.714 * Z + 0.0375 * Z ** (1 / 3) + 0.47

    eta_0_e = [0.733, 0.516, f_eta_0_e(Z)]
    hprime_0 = [3.348, 7.171, f_hprime_0(Z)]
    hprime_1 = [2.493, 5.884, f_hprime_1(Z)]
    hprime_2 = [2.519, 2.425, f_hprime_2(Z)]
    hprime_3 = [1.538, 3.527, f_hprime_3(Z)]
    hprime_4 = [0.039, -0.061, f_hprime_4(Z)]
    h_0 = [1.728, 3.979, f_h_0(Z)]
    h_1 = [1.030, 2.066, f_h_1(Z)]
    h_2 = [2.860, 3.864, f_h_2(Z)]
    h_3 = [0.100, 0.646, f_h_3(Z)]
    h_4 = [0.132, 0.054, f_h_4(Z)]
    h_5 = [1.405, 2.677, f_h_5(Z)]

    eta_0 = eta_0_e[Z_idx]

    def f_eta_2(Z_idx, r):
        numerator = ((6 / 5 * Z + 3 / 5 * np.sqrt(2)) * r +
                     hprime_0[Z_idx] * eta_0_e[Z_idx])
        denominator = (r ** 3 +
                       hprime_4[Z_idx] * r ** (7 / 3) +
                       hprime_3[Z_idx] * r ** 2 +
                       hprime_2[Z_idx] * r ** (5 / 3) +
                       hprime_1[Z_idx] * r +
                       hprime_0[Z_idx])
        return numerator / denominator

    eta_2 = f_eta_2(Z_idx, r)

    eta_1 = f_eta_2(Z_idx, 2 * r)

    def f_eta_4(Z_idx, r):
        numerator = r * (r + h_0[Z_idx] / h_5[Z_idx])
        denominator = (r ** 3 +
                       h_4[Z_idx] * r ** (7 / 3) +
                       h_3[Z_idx] * r ** 2 +
                       h_2[Z_idx] * r ** (5 / 3) +
                       h_1[Z_idx] * r +
                       h_0[Z_idx])
        return numerator / denominator

    eta_4 = f_eta_4(Z_idx, r)

    eta_3 = f_eta_4(Z_idx, 2 * r)

    return np.array((eta_0, eta_1, eta_2, eta_3, eta_4))


def _nondim_tc_i_ji_held(hall, Z, mu, theta, field_orientation, K=3):
    """dimensionless ion thermal conductivity - Ji-Held

    Ji, Jeong-Young, and Eric D. Held. "Closure and transport theory for
    high-collisionality electron-ion plasmas." Physics of Plasmas 20.4 (2013):
    042114.
    """

#    mu = m_e / m_i
#    theta = T_e / T_i
    zeta = 1 / Z * np.sqrt(mu / theta)
    r = np.abs(hall / np.sqrt(2))

#    K = 2  # 2x2 moments, equivalent to original Braginskii
#    K = 3  # 3x3 moments

    if K == 3:
        Delta_par_i1 = 1 + 26.90 * zeta + 187.5 * zeta ** 2 + 346.9 * zeta ** 3
        kappa_par_i = (5.586 + 101.7 * zeta + 289.1 * zeta ** 2) / Delta_par_i1
    elif K == 2:
        Delta_par_i1 = 1 + 13.50 * zeta + 36.46 * zeta ** 2
        kappa_par_i = (5.524 + 30.38 * zeta) / Delta_par_i1
    if field_orientation == 'parallel' or field_orientation == 'par':
        return kappa_par_i / np.sqrt(2)

    if K == 3:
        Delta_perp_i1 = (r ** 6 +
                         (3.635 + 29.15 * zeta + 83 * zeta ** 2) * r ** 4 +
                         (1.395 + 35.64 * zeta + 344.9 * zeta ** 2 +
                          1345 * zeta**3 + 1891 * zeta ** 4) * r ** 2 +
                          0.09163 * Delta_par_i1 ** 2)
        kappa_perp_i = ((np.sqrt(2) + 15 / 2 * zeta) * r ** 4 +
                        (3.841 + 57.59 * zeta + 297.8 * zeta ** 2 +
                         555 * zeta ** 3) * r ** 2 +
                        0.09163 * kappa_par_i * Delta_par_i1 ** 2
                        ) / Delta_perp_i1
    elif K == 2:
        Delta_perp_i1 = (r ** 4 +
                         (1.352 + 12.49 * zeta + 34 * zeta ** 2) * r ** 2 +
                         0.1693 * Delta_par_i1 ** 2)
        kappa_perp_i = ((np.sqrt(2) + 15 / 2 * zeta) * r ** 2 +
                        0.1693 * kappa_par_i * Delta_par_i1 ** 2
                        ) / Delta_perp_i1
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return kappa_perp_i / np.sqrt(2)

    if K == 3:
        kappa_cross_i = (r * (5 / 2 * r ** 4 +
                         (7.963 + 64.40 * zeta + 185 * zeta ** 2) * r ** 2 +
                         1.344 + 44.54 * zeta + 511.9 * zeta ** 2 +
                         2155 * zeta ** 3 + 3063 * zeta ** 4
                         ) / Delta_perp_i1)
    elif K == 2:
        kappa_cross_i = r * (5 / 2 * r ** 2 +
                             2.323 + 22.73 * zeta + 62.5 * zeta ** 2
                             ) / Delta_perp_i1
    if field_orientation == 'cross':
        return kappa_cross_i / np.sqrt(2)

    if field_orientation == 'all':
        return np.array((kappa_par_i / np.sqrt(2),
                         kappa_perp_i / np.sqrt(2),
                         kappa_cross_i / np.sqrt(2)))


def _nondim_visc_i_ji_held(hall, Z, mu, theta, K=3):
    """dimensionless ion viscosity - Ji-Held

    Ji, Jeong-Young, and Eric D. Held. "Closure and transport theory for
    high-collisionality electron-ion plasmas." Physics of Plasmas 20.4 (2013):
    042114.
    """

    zeta = 1 / Z * np.sqrt(mu / theta)
    r = np.abs(hall / np.sqrt(2))
    r13 = 2 * r

#    K = 2  # 2x2 moments, equivalent to original Braginskii
#    K = 3  # 3x3 moments

    if K == 3:
        Delta_par_i2 = 1 + 15.79 * zeta + 63.92 * zeta ** 2 + 71.69 * zeta ** 3
        eta_0_i = (1.365 + 16.75 * zeta + 35.84 * zeta ** 2) / Delta_par_i2

        def Delta_perp_i2(r, zeta, Delta_par_i2):
            Delta_perp_i2 = (r ** 6 +
                             (4.391 + 26.69 * zeta + 56 * zeta ** 2) * r ** 4 +
                             (3.191 + 49.62 * zeta + 306.4 * zeta ** 2 +
                              808.1 * zeta ** 3 + 784 * zeta ** 4) * r ** 2 +
                              0.4483 * Delta_par_i2 ** 2)
            return Delta_perp_i2

        Delta_perp_i2_24 = Delta_perp_i2(r, zeta, Delta_par_i2)
        Delta_perp_i2_13 = Delta_perp_i2(r13, zeta, Delta_par_i2)

        def f_eta_2(r, zeta, Delta_perp_i2):
            eta_2_i = (((3 / 5 * np.sqrt(2) + 2 * zeta) * r ** 4 +
                        (2.680 + 25.98 * zeta + 90.71 * zeta ** 2 +
                         104 * zeta ** 3) * r ** 2 +
                         0.4483 * eta_0_i * Delta_par_i2 ** 2
                       ) / Delta_perp_i2)
            return eta_2_i

        eta_2_i = f_eta_2(r, zeta, Delta_perp_i2_24)
        eta_1_i = f_eta_2(r13, zeta, Delta_perp_i2_13)

        def f_eta_4(r, zeta, Delta_perp_i2):
            eta_4_i = r * (r**4 +
                           (3.535 + 23.30 * zeta + 52 * zeta ** 2) * r ** 2 +
                           0.9538 + 21.81 * zeta + 174.2 * zeta ** 2 +
                           538.4 * zeta ** 3 + 576 * zeta ** 4
                           ) / Delta_perp_i2
            return eta_4_i

        eta_4_i = f_eta_4(r, zeta, Delta_perp_i2_24)
        eta_3_i = f_eta_4(r13, zeta, Delta_perp_i2_13)

    elif K == 2:
        Delta_par_i2 = 1 + 7.164 * zeta + 10.49 * zeta**2
        eta_0_i = (1.357 + 5.243 * zeta) / Delta_par_i2

        def Delta_perp_i2(r, zeta, Delta_par_i2):
            Delta_perp_i2 = (r ** 4 +
                             (2.023 + 11.68 * zeta + 20 * zeta ** 2) * r ** 2 +
                             0.5820 * Delta_par_i2 ** 2)
            return Delta_perp_i2

        Delta_perp_i2_24 = Delta_perp_i2(r, zeta, Delta_par_i2)
        Delta_perp_i2_13 = Delta_perp_i2(r13, zeta, Delta_par_i2)

        def f_eta_2(r, zeta, Delta_perp_i2):
            eta_2_i = ((3 / 5 * np.sqrt(2) + 2 * zeta) * r ** 2 +
                       0.5820 * eta_0_i * Delta_par_i2 ** 2
                       ) / Delta_perp_i2
            return eta_2_i

        eta_2_i = f_eta_2(r, zeta, Delta_perp_i2_24)
        eta_1_i = f_eta_2(r13, zeta, Delta_perp_i2_13)

        def f_eta_4(r, zeta, Delta_perp_i2):
            Delta_perp_i2 = (r ** 4 +
                             (2.023 + 11.68 * zeta + 20 * zeta**2) * r ** 2 +
                             0.5820 * Delta_par_i2 ** 2)
            eta_4_i = r * (r ** 2 +
                           1.188 + 8.283 * zeta + 16 * zeta ** 2
                           ) / Delta_perp_i2
            return eta_4_i

        eta_4_i = f_eta_4(r, zeta, Delta_perp_i2_24)
        eta_3_i = f_eta_4(r13, zeta, Delta_perp_i2_13)

    return np.array((eta_0_i / np.sqrt(2), eta_1_i / np.sqrt(2),
                     eta_2_i / np.sqrt(2), eta_3_i / np.sqrt(2),
                     eta_4_i / np.sqrt(2)))

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#                      end of classical_transport
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
