"""Functions to calculate plasma parameters."""

from astropy import units

from ..constants import (m_e, c, mu0, k_B, e, eps0)
import plasmapy.atomic as atomic
# from plasmapy.atomic import ion_mass, charge_state

import numpy as np
# import warnings

# For future: change these into decorators.  _check_quantity does a
# bit more than @quantity_input as it can allow
import plasmapy.utils as utils
from plasmapy.utils.checks import _check_quantity
from plasmapy.utils.exceptions import PhysicsError  # , PhysicsWarning


r"""
Values should be returned as an Astropy Quantity in SI units.

If a quantity has several names, then the function name should be the
one that provides the most physical insight into what the quantity
represents.  For example, 'gyrofrequency' indicates gyration, while
Larmor frequency indicates that this frequency is somehow related to a
human (or perhaps a cat?) named Larmor.  Similarly, using omega_ce as
a function name for this quantity will make the code less readable to
people who are unfamiliar with the notation or use a different symbol.

The docstrings for plasma parameter methods should describe the
physics associated with these quantities in ways that are
understandable to students who are taking their first course in plasma
physics while still being useful to experienced plasma physicists.

In many cases, units are enough to tell what field a quantity
represents.  The following line is an example.

>>> Alfven_speed(5*units.T, 8e-7*units.kg/units.m**3)

However, some plasma parameters depend on more than one quantity with
the same units.  In the following line, it is difficult to discern which
is the electron temperature and which is the ion temperature.

>>> ion_sound_speed(1e6*units.K, 2e6*units.K)

Remembering that "explicit is better than implicit", it is more
readable and less prone to errors to write:

>>> ion_sound_speed(T_i=1e6*units.K, T_e=2e6*units.K)

SI units that were named after a person should be lower case except at
the beginning of a sentence, even if their symbol is capitalized. For
example, kelvin is a unit while Kelvin was a scientist.

Unit conversions involving angles must be treated with care.  Angles
are dimensionless but do have units.  Angular velocity is often
given in units of radians per second, though dimensionally this is
equivalent to inverse seconds.  Astropy will treat radians
dimensionlessly when using the dimensionless_angles() equivalency,
but dimensionless_angles() does not account for multiplicative
factor of 2*pi that is used when converting between frequency (1 /
s) and angular velocity (rad / s).  A clear way to do this
conversion is to set up an equivalency between cycles/s and Hz:

>>> from astropy import units
>>> f_ce = omega_ce.to(units.Hz, equivalencies=[(units.cy/units.s, units.Hz)])

However, dimensionless_angles() does work when dividing a velocity
by an angular frequency to get a length scale:

>>> d_i = (c/omega_pi).to(units.m, equivalencies=units.dimensionless_angles())

"""


@utils.check_relativistic
def Alfven_speed(B, density, ion="p", z_mean=None):
    r"""
    Returns the Alfven speed.

    Parameters
    ----------
    B : ~astropy.units.Quantity
        The magnetic field magnitude in units convertible to tesla.

    density : ~astropy.units.Quantity
        Either the ion number density in units convertible to 1 / m**3,
        or the mass density in units convertible to kg / m**3.

    ion : string, optional
        Representation of the ion species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to protons.  If no charge state information is
        provided, then the ions are assumed to be singly charged.

    z_mean : Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. If this quantity is not
        given then the atomic charge state (integer) of the ion
        is used. This is effectively an average Alfven speed for the
        plasma where multiple charge states are present.

    Returns
    -------
    V_A : Quantity with units of velocity
        The Alfven velocity of the plasma in units of meters per second.

    Raises
    ------
    TypeError
        The magnetic field and density arguments are not Quantities and
        cannot be converted into Quantities.

    UnitConversionError
        If the magnetic field or density is not in appropriate units.

    RelativityError
        If the Alfven velocity is greater than or equal to the speed of light

    ValueError
        If the density is negative, or the ion mass or charge state
        cannot be found.

    UserWarning
        if units are not provided and SI units are assumed.

    Warns
    -----
    RelativityWarning
        If the Alfven velocity exceeds 10% of the speed of light


    Notes
    -----
    The Alfven velocity :math:`V_A` is the typical propagation speed
    of magnetic disturbances in a plasma, and is given by:

    .. math::

        V_A = \frac{B}{\sqrt{\mu_0\rho}}

    where the mass density is :math:`\rho = n_i m_i + n_e m_e`.

    This expression does not account for relativistic effects, and
    loses validity when the resulting speed is a significant fraction
    of the speed of light.

    This function switches B and density when B has units of number
    density or mass density and density has units of magnetic field
    strength.

    Examples
    --------
    >>> from astropy import units as u
    >>> from plasmapy.constants import m_p, m_e
    >>> B = 0.014*u.T
    >>> n = 5e19*u.m**-3
    >>> rho = n*(m_p+m_e)
    >>> ion = 'p'
    >>> Alfven_speed(B, n, ion)
    <Quantity 43173.87029559 m / s>
    >>> Alfven_speed(B, rho, ion)
    <Quantity 43173.87029559 m / s>
    >>> Alfven_speed(B, rho, ion).to(u.cm/u.us)
    <Quantity 4.31738703 cm / us>
    """

    _check_quantity(B, 'B', 'Alfven_speed', units.T)
    _check_quantity(density, 'density', 'Alfven_speed',
                    [units.m**-3, units.kg / units.m**3], can_be_negative=False)

    B = B.to(units.T)
    density = density.si

    if density.unit == units.m**-3:
        try:
            m_i = atomic.ion_mass(ion)
            if z_mean is None:
                # warnings.warn("No z_mean given, defaulting to atomic charge",
                #               PhysicsWarning)
                try:
                    Z = atomic.charge_state(ion)
                except ValueError:
                    Z = 1
            else:
                # using average ionization provided by user
                Z = z_mean
            # ensuring positive value of Z
            Z = np.abs(Z)
        except Exception:
            raise ValueError("Invalid ion in Alfven_speed.")
        rho = density * m_i + Z * density * m_e

    elif density.unit == units.kg / units.m**3:
        rho = density

    try:
        V_A = (np.abs(B) / np.sqrt(mu0 * rho)).to(units.m / units.s)
    except Exception:
        raise ValueError("Unable to find Alfven speed")

    return V_A


@utils.check_relativistic
@utils.check_quantity({
    'T_i': {'units': units.K, 'can_be_negative': False},
    'T_e': {'units': units.K, 'can_be_negative': False}
})
def ion_sound_speed(*ignore,
                    T_e=0 * units.K,
                    T_i=0 * units.K,
                    gamma_e=1,
                    gamma_i=3,
                    ion='p',
                    z_mean=None):
    r"""
    Returns the ion sound speed for an electron-ion plasma.

    Parameters
    ----------
    T_e : ~astropy.units.Quantity, optional
        Electron temperature in units of temperature or energy per
        particle.  If this is not given, then the electron temperature
        is assumed to be zero.  If only one temperature is entered, it
        is assumed to be the electron temperature.

    T_i : ~astropy.units.Quantity, optional
        Ion temperature in units of temperature or energy per
        particle.  If this is not given, then the ion temperature is
        assumed to be zero.

    gamma_e : float or int
        The adiabatic index for electrons, which defaults to 1.  This
        value assumes that the electrons are able to equalize their
        temperature rapidly enough that the electrons are effectively
        isothermal.

    gamma_i : float or int
        The adiabatic index for ions, which defaults to 3.  This value
        assumes that ion motion has only one degree of freedom, namely
        along magnetic field lines.

    ion : string, optional
        Representation of the ion species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to protons.  If no charge state information is
        provided, then the ions are assumed to be singly charged.

    z_mean : Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. If this quantity is not
        given then the atomic charge state (integer) of the ion
        is used. This is effectively an average ion sound speed for the
        plasma where multiple charge states are present.

    Returns
    -------
    V_S : ~astropy.units.Quantity
        The ion sound speed in units of meters per second.

    Raises
    ------
    TypeError
        If any of the arguments are not entered as keyword arguments
        or are of an incorrect type.

    ValueError
        If the ion mass, adiabatic index, or temperature are invalid.

    PhysicsError
        If an adiabatic index is less than one.

    UnitConversionError
        If the temperature is in incorrect units.

    UserWarning
        If the ion sound speed exceeds 10% of the speed of light, or
        if units are not provided and SI units are assumed.

    Notes
    -----
    The ion sound speed :math:`V_S` is approximately given by

    .. math::

        V_S = \sqrt{\frac{\gamma_e Z k_B T_e + \gamma_i k_B T_i}{m_i}}

    where :math:`\gamma_e` and :math:`\gamma_i` are the electron and
    ion adiabatic indices, :math:`k_B` is the Boltzmann constant,
    :math:`T_e` and :math:`T_i` are the electron and ion temperatures,
    :math:`Z` is the charge state of the ion, and :math:`m_i` is the
    ion mass.

    This function assumes that the product of the wavenumber and the
    Debye length is small. In this limit, the ion sound speed is not
    dispersive (e.g., frequency independent).

    When the electron temperature is much greater than the ion
    temperature, the ion sound velocity reduces to
    :math:`\sqrt{\gamma_e k_B T_e / m_i}`.  Ion acoustic waves can
    therefore occur even when the ion temperature is zero.

    Example
    -------
    >>> from astropy import units as u
    >>> ion_sound_speed(T_e=5e6*u.K, T_i=0*u.K, ion='p', gamma_e=1, gamma_i=3)
    <Quantity 203155.0764042 m / s>
    >>> ion_sound_speed(T_e=5e6*u.K)
    <Quantity 203155.0764042 m / s>
    >>> ion_sound_speed(T_e=500*u.eV, T_i=200*u.eV, ion='D+')
    <Quantity 229586.01860212 m / s>

    """

    if ignore:
        raise TypeError("All arguments are required to be keyword arguments "
                        "in ion_sound_speed to prevent mixing up the electron "
                        "and ion temperatures. An example call that uses the "
                        "units subpackage from astropy is: "
                        "ion_sound_speed(T_e=5*units.K, T_i=0*units.K, "
                        "ion='D+')")

    try:
        m_i = atomic.ion_mass(ion)
        if z_mean is None:
            # warnings.warn("No z_mean given, defaulting to atomic charge",
            #               PhysicsWarning)
            try:
                Z = atomic.charge_state(ion)
            except ValueError:
                Z = 1
        else:
            # using average ionization provided by user
            Z = z_mean
    except Exception:
        raise ValueError("Invalid ion in ion_sound_speed.")

    if not isinstance(gamma_e, (float, int)):
        raise TypeError("The adiabatic index for electrons (gamma_e) must be "
                        "a float or int in ion_sound_speed")
    if not isinstance(gamma_i, (float, int)):
        raise TypeError("The adiabatic index for ions (gamma_i) must be "
                        "a float or int in ion_sound_speed")

    if not 1 <= gamma_e <= np.inf:
        raise PhysicsError("The adiabatic index for electrons must be between "
                           "one and infinity")
    if not 1 <= gamma_i <= np.inf:
        raise PhysicsError("The adiabatic index for ions must be between "
                           "one and infinity")

    T_i = T_i.to(units.K, equivalencies=units.temperature_energy())
    T_e = T_e.to(units.K, equivalencies=units.temperature_energy())

    try:
        V_S_squared = (gamma_e * Z * k_B * T_e + gamma_i * k_B * T_i) / m_i
        V_S = np.sqrt(V_S_squared).to(units.m / units.s)
    except Exception:
        raise ValueError("Unable to find ion sound speed.")

    return V_S


@utils.check_relativistic
@utils.check_quantity({
    'T': {'units': units.K, 'can_be_negative': False}
})
def thermal_speed(T, particle="e", method="most_probable"):
    r"""
    Returns the most probable speed for a particle within a Maxwellian
    distribution.

    Parameters
    ----------
    T : ~astropy.units.Quantity
        The particle temperature in either kelvin or energy per particle

    particle : string, optional
        Representation of the particle species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to electrons.  If no charge state information is
        provided, then the particles are assumed to be singly charged.

    method : string, optional
        Method to be used for calculating the thermal speed. Options are
        'most_probable' (default), 'rms', and 'mean_magnitude'.

    Returns
    -------
    V : ~astropy.units.Quantity
        particle thermal speed

    Raises
    ------
    TypeError
        The particle temperature is not a Quantity

    UnitConversionError
        If the particle temperature is not in units of temperature or
        energy per particle

    ValueError
        The particle temperature is invalid or particle cannot be used to
        identify an isotope or particle

    UserWarning
        If the particle thermal speed exceeds 10% of the speed of light, or
        if units are not provided and SI units are assumed.

    Notes
    -----
    The particle thermal speed is given by:

    .. math::
        V_{th,i} = \sqrt{\frac{2 k_B T_i}{m_i}}

    This function yields the most probable speed within a distribution
    function.  However, the definition of thermal velocity varies by
    the square root of two depending on whether or not this velocity
    absorbs that factor in the expression for a Maxwellian
    distribution.  In particular, the expression given in the NRL
    Plasma Formulary [1] is a square root of two smaller than the
    result from this function.

    Examples
    --------
    >>> from astropy import units as u
    >>> thermal_speed(5*u.eV, 'p')
    <Quantity 30949.69018286 m / s>
    >>> thermal_speed(1e6*u.K, particle='p')
    <Quantity 128486.55193256 m / s>
    >>> thermal_speed(5*u.eV)
    <Quantity 1326205.12123959 m / s>
    >>> thermal_speed(1e6*u.K)
    <Quantity 5505693.98842538 m / s>
    >>> thermal_speed(1e6*u.K, method="rms")
    <Quantity 6743070.47577549 m / s>
    >>> thermal_speed(1e6*u.K, method="mean_magnitude")
    <Quantity 6212510.3969422 m / s>

    """

    T = T.to(units.K, equivalencies=units.temperature_energy())

    try:
        m = atomic.ion_mass(particle)
    except Exception:
        raise ValueError("Unable to find {particle} mass in thermal_speed")

    # different methods, as per https://en.wikipedia.org/wiki/Thermal_velocity
    if method == "most_probable":
        V = (np.sqrt(2 * k_B * T / m)).to(units.m / units.s)
    elif method == "rms":
        V = (np.sqrt(3 * k_B * T / m)).to(units.m / units.s)
    elif method == "mean_magnitude":
        V = (np.sqrt(8 * k_B * T / (m * np.pi))).to(units.m / units.s)
    else:
        raise ValueError("Method {method} not supported in thermal_speed")

    return V


@utils.check_relativistic
@utils.check_quantity({
    'T': {'units': units.K, 'can_be_negative': False}
})
def kappa_thermal_speed(T, kappa, particle="e", method="most_probable"):
    r"""Returns the most probable speed for a particle within a Kappa
    distribution.

    Parameters
    ----------
    T : ~astropy.units.Quantity
        The particle temperature in either kelvin or energy per particle

    kappa: ~astropy.units.Quantity
        The kappa parameter is a dimensionless number which sets the slope
        of the energy spectrum of suprathermal particles forming the tail
        of the Kappa velocity distribution function. Kappa must be greater
        than 3/2.

    particle : string, optional
        Representation of the particle species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to electrons.  If no charge state information is
        provided, then the particles are assumed to be singly charged.

    method : string, optional
        Method to be used for calculating the thermal speed. Options are
        'most_probable' (default), 'rms', and 'mean_magnitude'.

    Returns
    -------
    V : ~astropy.units.Quantity
        particle thermal speed

    Raises
    ------
    TypeError
        The particle temperature is not a Quantity

    UnitConversionError
        If the particle temperature is not in units of temperature or
        energy per particle

    ValueError
        The particle temperature is invalid or particle cannot be used to
        identify an isotope or particle

    UserWarning
        If the particle thermal speed exceeds 10% of the speed of light, or
        if units are not provided and SI units are assumed.

    Notes
    -----
    The particle thermal speed is given by:

    .. math::
        V_{th,i} = \sqrt{(2 \kappa - 3)\frac{2 k_B T_i}{\kappa m_i}}

    Examples
    --------
    >>> from astropy import units as u
    >>> kappa_thermal_speed(5*u.eV, 4, 'p')
    <Quantity 24467.87846359 m / s>
    """
    # Checking thermal units
    T = T.to(units.K, equivalencies=units.temperature_energy())
    if kappa <= 3 / 2:
        raise ValueError(f"Must have kappa > 3/2, instead of {kappa}, for "
                         "kappa distribution function to be valid.")
    # different methods, as per https://en.wikipedia.org/wiki/Thermal_velocity
    if method == "most_probable":
        # thermal velocity of Kappa distribution function is just Maxwellian
        # thermal speed modulated by the following factor.
        # This is only true for "most probable" case. RMS and mean
        # magnitude velocities are same as Maxwellian.
        coeff = np.sqrt((kappa - 3 / 2) / kappa)
        vTh = thermal_speed(T=T,
                            particle=particle,
                            method="most_probable")
        return coeff * vTh
    elif method == "rms":
        vTh = thermal_speed(T=T,
                            particle=particle,
                            method="rms")
        return vTh
    elif method == "mean_magnitude":
        vTh = thermal_speed(T=T,
                            particle=particle,
                            method="mean_magnitude")
        return vTh
    else:
        raise ValueError("Method {method} not supported in thermal_speed")


@utils.check_quantity({
    'B': {'units': units.T}
})
def gyrofrequency(B, particle='e', signed=False, z_mean=None):
    r"""Calculate the particle gyrofrequency in units of radians per second.

    Parameters
    ----------
    B : ~astropy.units.Quantity
        The magnetic field magnitude in units convertible to tesla.

    particle : string, optional
        Representation of the particle species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to electrons.  If no charge state information is
        provided, then the particles are assumed to be singly charged.

    z_mean : Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. If this quantity is not
        given then the atomic charge state (integer) of the ion
        is used. This is effectively an average gyrofrequency for the
        plasma where multiple charge states are present, and should
        not be interpreted as the gyrofrequency for any single particle.

    signed : boolean, optional
        The gyrofrequency can be defined as signed (negative for electron,
        positive for ion). Default is False (unsigned, ie. always positive).

    Returns
    -------
    omega_c : ~astropy.units.Quantity
        The particle gyrofrequency in units of radians per second

    Raises
    ------
    TypeError
        If the magnetic field is not a Quantity or particle is not of an
        appropriate type

    ValueError
        If the magnetic field contains invalid values or particle cannot be
        used to identify an particle or isotope

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----
    The particle gyrofrequency is the angular frequency of particle gyration
    around magnetic field lines and is given by:

    .. math::
        \omega_{ci} = \frac{Z e B}{m_i}

    The particle gyrofrequency is also known as the particle cyclotron
    frequency or the particle Larmor frequency.

    The recommended way to convert from angular frequency to frequency
    is to use an equivalency between cycles per second and Hertz, as
    Astropy's dimensionles_angles() equivalency does not account for
    the factor of 2*pi needed during this conversion.  The
    dimensionless_angles() equivalency is appropriate when dividing a
    velocity by an angular frequency to get a length scale.

    Examples
    --------
    >>> from numpy import pi
    >>> from astropy import units as u
    >>> gyrofrequency(0.1*u.T)
    <Quantity 1.75882002e+10 rad / s>
    >>> gyrofrequency(0.1*u.T, signed=True)
    <Quantity -1.75882002e+10 rad / s>
    >>> gyrofrequency(0.01*u.T, 'p')
    <Quantity 957883.32241481 rad / s>
    >>> gyrofrequency(0.01*u.T, 'p')
    <Quantity 957883.32241481 rad / s>
    >>> gyrofrequency(0.01*u.T, particle='T+')
    <Quantity 319964.54975911 rad / s>
    >>> omega_ce = gyrofrequency(0.1*u.T)
    >>> print(omega_ce)
    17588200236.02124 rad / s
    >>> f_ce = omega_ce.to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])
    >>> print(f_ce)
    2799249007.6528206 Hz

    """
    try:
        m_i = atomic.ion_mass(particle)
        if z_mean is None:
            # warnings.warn("No z_mean given, defaulting to atomic charge",
            #               PhysicsWarning)
            try:
                Z = atomic.charge_state(particle)
            except ValueError:
                Z = 1
        else:
            # using user provided average ionization
            Z = z_mean
    except Exception:
        raise ValueError("Invalid particle {} in gyrofrequency"
                         .format(particle))
    if not signed:
        Z = abs(Z)

    omega_ci = units.rad * (Z * e * np.abs(B) / m_i).to(1 / units.s)

    return omega_ci


def gyroradius(B, *args, Vperp=None, T_i=None, particle='e'):
    r"""Returns the particle gyroradius.

    Parameters
    ----------
    B : ~astropy.units.Quantity
        The magnetic field magnitude in units convertible to tesla.

    Vperp : ~astropy.units.Quantity, optional
        The component of particle velocity that is perpendicular to the
        magnetic field in units convertible to meters per second.

    T_i : ~astropy.units.Quantity, optional
        The particle temperature in units convertible to kelvin.

    particle : string, optional
        Representation of the particle species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to electrons.  If no charge state information is
        provided, then the particles are assumed to be singly charged.

    args : ~astropy.units.Quantity
        If the second positional argument is a Quantity with units
        appropriate to Vperp or T_i, then this argument will take the
        place of that keyword argument.

    Returns
    -------
    r_Li : ~astropy.units.Quantity
        The particle gyroradius in units of meters.  This Quantity will be
        based on either the perpendicular component of particle velocity as
        inputted, or the most probable speed for an particle within a
        Maxwellian distribution for the particle temperature.

    Raises
    ------
    TypeError
        The arguments are of an incorrect type

    UnitConversionError
        The arguments do not have appropriate units

    ValueError
        If any argument contains invalid values

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----
    One but not both of Vperp and T_i must be inputted.

    If any of B, Vperp, or T_i is a number rather than a Quantity,
    then SI units will be assumed and a warning will be raised.

    The particle gyroradius is also known as the particle Larmor radius and is
    given by

    .. math::
        r_{Li} = \frac{V_{\perp}}{omega_{ci}}

    where :math:`V_{\perp}` is the component of particle velocity that is
    perpendicular to the magnetic field and :math:`\omega_{ci}` is the
    particle gyrofrequency.  If a temperature is provided, then
    :math:`V_\perp` will be the most probable thermal velocity of an
    particle at that temperature.

    Examples
    --------
    >>> from astropy import units as u
    >>> gyroradius(0.2*u.T, 1e5*u.K, particle='p')
    <Quantity 0.00212087 m>
    >>> gyroradius(0.2*u.T, 1e5*u.K, particle='p')
    <Quantity 0.00212087 m>
    >>> gyroradius(5*u.uG, 1*u.eV, particle='alpha')
    <Quantity 288002.38837768 m>
    >>> gyroradius(400*u.G, 1e7*u.m/u.s, particle='Fe+++')
    <Quantity 48.23129811 m>
    >>> gyroradius(B = 0.01*u.T, T_i = 1e6*u.K)
    <Quantity 0.00313033 m>
    >>> gyroradius(B = 0.01*u.T, Vperp = 1e6*u.m/u.s)
    <Quantity 0.00056856 m>
    >>> gyroradius(0.2*u.T, 1e5*u.K)
    <Quantity 4.94949252e-05 m>
    >>> gyroradius(5*u.uG, 1*u.eV)
    <Quantity 6744.2598183 m>
    >>> gyroradius(400*u.G, 1e7*u.m/u.s)
    <Quantity 0.00142141 m>

    """

    if Vperp is not None and T_i is not None:
        raise ValueError("Cannot have both Vperp and T_i as arguments to "
                         "gyroradius")

    if len(args) == 1 and isinstance(args[0], units.Quantity):
        arg = args[0].si
        if arg.unit == units.T and B.si.unit in [units.J, units.K,
                                                 units.m / units.s]:
            B, arg = arg, B

        if arg.unit == units.m / units.s:
            Vperp = arg
        elif arg.unit in (units.J, units.K):
            T_i = arg.to(units.K, equivalencies=units.temperature_energy())
        else:
            raise units.UnitConversionError("Incorrect units for positional "
                                            "argument in gyroradius")
    elif len(args) > 0:
        raise ValueError("Incorrect inputs to gyroradius")

    _check_quantity(B, 'B', 'gyroradius', units.T)

    if Vperp is not None:
        _check_quantity(Vperp, 'Vperp', 'gyroradius', units.m / units.s)
    elif T_i is not None:
        _check_quantity(T_i, 'T_i', 'gyroradius', units.K)
        Vperp = thermal_speed(T_i, particle=particle)

    omega_ci = gyrofrequency(B, particle)

    r_Li = np.abs(Vperp) / omega_ci

    return r_Li.to(units.m, equivalencies=units.dimensionless_angles())


@utils.check_quantity({
    'n': {'units': units.m**-3, 'can_be_negative': False}
})
def plasma_frequency(n, particle='e', z_mean=None):
    r"""Calculates the particle plasma frequency.

    Parameters
    ----------
    n : ~astropy.units.Quantity
        Particle number density in units convertible to per cubic meter

    particle : string, optional
        Representation of the particle species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to electrons.  If no charge state information is
        provided, then the particles are assumed to be singly charged.

    z_mean : Quantity, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. If this quantity is not
        given then the atomic charge state (integer) of the ion
        is used. This is effectively an average plasma frequency for the
        plasma where multiple charge states are present. 

    Returns
    -------
    omega_p : ~astropy.units.Quantity
        The particle plasma frequency in radians per second.

    Raises
    ------
    TypeError
        If n_i is not a Quantity or particle is not of an appropriate type

    UnitConversionError
        If n_i is not in correct units

    ValueError
        If n_i contains invalid values or particle cannot be used to
        identify an particle or isotope.

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----
    The particle plasma frequency is

    .. math::
        \omega_{pi} = Z e \sqrt{\frac{n_i}{\epsilon_0 m_i}}

    At present, astropy.units does not allow direct conversions from
    radians/second for angular frequency to 1/second or Hz for
    frequency.  The dimensionless_angles equivalency allows that
    conversion, but does not account for the factor of 2*pi. The
    alternatives are to convert to cycle/second or to do the
    conversion manually, as shown in the examples.

    Example
    -------
    >>> from astropy import units as u
    >>> plasma_frequency(1e19*u.m**-3, particle='p')
    <Quantity 4.16329453e+09 rad / s>
    >>> plasma_frequency(1e19*u.m**-3, particle='D+')
    <Quantity 2.94462452e+09 rad / s>
    >>> plasma_frequency(1e19*u.m**-3)
    <Quantity 1.78398636e+11 rad / s>
    """

    try:
        m = atomic.ion_mass(particle)
        if z_mean is None:
            # warnings.warn("No z_mean given, defaulting to atomic charge",
            #               PhysicsWarning)
            try:
                Z = atomic.charge_state(particle)
            except ValueError:
                Z = 1
        else:
            # using user provided average ionization
            Z = z_mean
        Z = np.abs(Z)
    except Exception:
        raise ValueError(f"Invalid particle, {particle}, in "
                         "plasma_frequency.")

    omega_p = units.rad * Z * e * np.sqrt(n / (eps0 * m))

    return omega_p.si


@utils.check_quantity({
    'T_e': {'units': units.K, 'can_be_negative': False},
    'n_e': {'units': units.m**-3, 'can_be_negative': False}
})
def Debye_length(T_e, n_e):
    r"""Calculate the Debye length.

    Parameters
    ----------
    T_e: ~astropy.units.Quantity
        Electron temperature

    n_e: ~astropy.units.Quantity
        Electron number density

    Returns
    -------
    lambda_D : ~astropy.units.Quantity
        The Debye length in meters

    Raises
    ------
    TypeError
        If either argument is not a ~astropy.units.Quantity

    UnitConversionError
        If either argument is in incorrect units

    ValueError
        If either argument contains invalid values

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----

    The Debye length is the exponential scale length for charge
    screening and is given by

    .. math::
        \lambda_D = \sqrt{\frac{\epsilon_0 k_b T_e}{n_e e^2}}

    for an electron plasma with nearly stationary ions.

    The electrical potential will drop by a factor of 1/e every Debye
    length.

    Plasmas will generally be quasineutral on length scales significantly
    larger than the Debye length.

    See also
    --------
    Debye_number

    Example
    -------
    >>> from astropy import units as u
    >>> Debye_length(5e6*u.K, 5e15*u.m**-3)
    <Quantity 0.00218226 m>

    """

    T_e = T_e.to(units.K, equivalencies=units.temperature_energy())

    try:
        lambda_D = ((eps0 * k_B * T_e / (n_e * e ** 2)) ** 0.5).to(units.m)
    except Exception:
        raise ValueError("Unable to find Debye length.")

    return lambda_D


@utils.check_quantity({
    'T_e': {'units': units.K, 'can_be_negative': False},
    'n_e': {'units': units.m**-3, 'can_be_negative': False}
})
def Debye_number(T_e, n_e):
    r"""Returns the Debye number.

    Parameters
    ----------
    T_e : ~astropy.units.Quantity
        Electron temperature

    n_e : ~astropy.units.Quantity
        Electron number density

    Raises
    ------
    TypeError
        If either argument is not a Quantity

    UnitConversionError
        If either argument is in incorrect units

    ValueError
        If either argument contains invalid values

    UserWarning
        If units are not provided and SI units are assumed

    Returns
    -------
    N_D : ~astropy.units.Quantity
        Number of electrons within a sphere with a radius of the Debye length

    Notes
    -----
    The Debye number is the number of electrons contained within a sphere with
    a radius of a Debye length and is given by

    .. math::
        N_D = \frac{4\pi}{3}n_e\lambda_D^3

    The Debye number is also known as the plasma parameter.

    Collective behavior requires a Debye number significantly larger than one.

    See also
    --------
    Debye_length

    Example
    -------
    >>> from astropy import units as u
    >>> Debye_number(5e6*u.K, 5e9*u.cm**-3)
    <Quantity 2.17658302e+08>

    """

    try:
        lambda_D = Debye_length(T_e, n_e)
        N_D = (4 / 3) * np.pi * n_e * lambda_D ** 3
    except Exception:
        raise ValueError("Unable to find Debye number")

    return N_D.to(units.dimensionless_unscaled)


@utils.check_quantity({
    'n': {'units': units.m**-3, 'can_be_negative': False}
})
def inertial_length(n, particle='e'):
    r"""Calculate the particle inertial length,

    Parameters
    ----------
    n_i : ~astropy.units.Quantity
        Particle number density in units convertible to m**-3.

    particle : string, optional
        Representation of the particle species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to electrons.  If no charge state information is
        provided, then the particles are assumed to be singly charged.

    Returns
    -------
    d_i : ~astropy.units.Quantity
        Particles inertial length in meters.

    Raises
    ------
    TypeError
        If n_i not a Quantity or particle is not a string.

    UnitConversionError
        If n_i is not in units of a number density.

    ValueError
        The particle density does not have an appropriate value.

    UserWarning
        If units are not provided and SI units are assumed.

    Notes
    -----
    The particle inertial length is also known as an particle skin depth and is
    given by:

    .. math::
        d_i = \frac{c}{\omega_{pi}}

    Example
    -------
    >>> from astropy import units as u
    >>> inertial_length(5*u.m**-3, particle='He+')
    <Quantity 2.02985802e+08 m>
    >>> inertial_length(5*u.m**-3)
    <Quantity 2376534.75601976 m>

    """

    try:
        Z = atomic.charge_state(particle)
    except Exception:
        raise ValueError(f"Invalid particle {particle} in inertial_length.")
    if Z:
        Z = abs(Z)

    omega_p = plasma_frequency(n, particle=particle)
    d = (c / omega_p).to(units.m, equivalencies=units.dimensionless_angles())

    return d


@utils.check_quantity({
    'B': {'units': units.T}
})
def magnetic_pressure(B):
    r"""
    Calculate the magnetic pressure.

    Parameters
    ----------
    B : ~astropy.units.Quantity
        The magnetic field in units convertible to tesla.

    Returns
    -------
    p_B : ~astropy.units.Quantity
        The magnetic pressure in units in pascals (newtons per square meter).

    Raises
    ------
    TypeError
        If the input is not a Quantity.

    UnitConversionError
        If the input is not in units convertible to tesla.

    ValueError
        If the magnetic field strength is not a real number between.
        +/- infinity

    UserWarning
        If units are not provided and SI units are assumed.

    Notes
    -----
    The magnetic pressure is given by:

    .. math::
        p_B = \frac{B^2}{2 \mu_0}

    The motivation behind having two separate functions for magnetic
    pressure and magnetic energy density is that it allows greater
    insight into the physics that are being considered by the user and
    thus more readable code.

    See also
    --------
    magnetic_energy_density : returns an equivalent Quantity, except in
        units of joules per cubic meter.

    Example
    -------
    >>> from astropy import units as u
    >>> magnetic_pressure(0.1*u.T).to(u.Pa)
    <Quantity 3978.8735773 Pa>

    """

    p_B = (B ** 2 / (2 * mu0)).to(units.Pa)

    return p_B


@utils.check_quantity({
    'B': {'units': units.T}
})
def magnetic_energy_density(B: units.T):
    r"""
    Calculate the magnetic energy density.

    Parameters
    ----------
    B : ~astropy.units.Quantity
        The magnetic field in units convertible to tesla.

    Returns
    -------
    E_B : ~astropy.units.Quantity
        The magnetic energy density in units of joules per cubic meter.

    Raises
    ------
    TypeError
        If the input is not a Quantity.

    UnitConversionError
        If the input is not in units convertible to tesla.

    ValueError
        If the magnetic field strength does not have an appropriate.
        value.

    UserWarning
        If units are not provided and SI units are assumed.

    Notes
    -----
    The magnetic energy density is given by:

    .. math::
        E_B = \frac{B^2}{2 \mu_0}

    The motivation behind having two separate functions for magnetic
    pressure and magnetic energy density is that it allows greater
    insight into the physics that are being considered by the user and
    thus more readable code.

    See also
    --------
    magnetic_pressure : Returns an equivalent Quantity, except in units
        of pascals.

    Example
    -------
    >>> from astropy import units as u
    >>> magnetic_energy_density(0.1*u.T)
    <Quantity 3978.8735773 J / m3>

    """

    E_B = (B ** 2 / (2 * mu0)).to(units.J / units.m ** 3)

    return E_B


@utils.check_quantity({
    'B': {'units': units.T},
    'n_e': {'units': units.m**-3, 'can_be_negative': False}
})
def upper_hybrid_frequency(B, n_e):
    r"""
    Returns the upper hybrid frequency.

    Parameters
    ----------
    B : ~astropy.units.Quantity
        The magnetic field magnitude in units convertible to tesla.

    n_e : ~astropy.units.Quantity
        The electron number density.

    Returns
    -------
    omega_uh : ~astropy.units.Quantity
        The upper hybrid frequency in radians per second.

    Raises
    ------
    TypeError
        If either of B or n_e is not a Quantity.

    UnitConversionError
        If either of B or n_e is in incorrect units.

    ValueError
        If either of B or n_e contains invalid values or are of
        incompatible dimensions.

    UserWarning
        If units are not provided and SI units are assumed.

    Notes
    -----
    The upper hybrid frequency is given through the relation

    .. math::
        \omega_{uh}^2 = \omega_{ce}^2 + \omega_{pe}^2

    where :math:`\omega_{ce}` is the electron gyrofrequency and
    :math:`\omega_{pe}` is the electron plasma frequency.

    Example
    -------
    >>> from astropy import units as u
    >>> upper_hybrid_frequency(0.2*u.T, n_e=5e19*u.m**-3)
    <Quantity 4.00459419e+11 rad / s>

    """

    try:
        omega_pe = plasma_frequency(n=n_e)
        omega_ce = gyrofrequency(B)
        omega_uh = (np.sqrt(omega_pe**2 + omega_ce**2)).to(units.rad / units.s)
    except Exception:
        raise ValueError("Unable to find upper hybrid frequency.")

    return omega_uh


@utils.check_quantity({
    'B': {'units': units.T},
    'n_i': {'units': units.m**-3, 'can_be_negative': False}
})
def lower_hybrid_frequency(B, n_i, ion='p'):
    r"""
    Returns the lower hybrid frequency.

    Parameters
    ----------
    B : ~astropy.units.Quantity
        The magnetic field magnitude in units convertible to tesla.

    n_i : ~astropy.units.Quantity
        Ion number density.

    ion : string, optional
        Representation of the ion species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to protons.  If no charge state information is
        provided, then the ions are assumed to be singly charged.

    Returns
    -------
    omega_lh : ~astropy.units.Quantity
        The lower hybrid frequency in radians per second.

    Raises
    ------
    TypeError
        If either of B or n_i is not a Quantity, or ion is of an
        inappropriate type.

    UnitConversionError
        If either of B or n_i is in incorrect units.

    ValueError
        If either of B or n_i contains invalid values or are of
        incompatible dimensions, or ion cannot be used to identify an
        ion or isotope.

    UserWarning
        If units are not provided and SI units are assumed.

    Notes
    -----
    The lower hybrid frequency is given through the relation

    .. math::
        \frac{1}{\omega_{lh}^2} = \frac{1}{\omega_{ci}^2 + \omega_{pi}^2} +
        \frac{1}{\omega_{ci}\omega_{ce}}

    where :math:`\omega_{ci}` is the ion gyrofrequency,
    :math:`\omega_{ce}` is the electron gyrofrequency, and
    :math:`\omega_{pi}` is the ion plasma frequency.

    Example
    -------
    >>> from astropy import units as u
    >>> lower_hybrid_frequency(0.2*u.T, n_i=5e19*u.m**-3, ion='D+')
    <Quantity 5.78372733e+08 rad / s>

    """

    # We do not need a charge state here, so the sole intent is to
    # catch invalid ions.
    try:
        atomic.charge_state(ion)
    except Exception:
        raise ValueError("Invalid ion in lower_hybrid_frequency.")

    try:
        omega_ci = gyrofrequency(B, particle=ion)
        omega_pi = plasma_frequency(n_i, particle=ion)
        omega_ce = gyrofrequency(B)
        omega_lh = 1 / np.sqrt((omega_ci * omega_ce) ** -1 + omega_pi ** -2)
        omega_lh = omega_lh.to(units.rad / units.s)
    except Exception:
        raise ValueError("Unable to find lower hybrid frequency.")

    return omega_lh
