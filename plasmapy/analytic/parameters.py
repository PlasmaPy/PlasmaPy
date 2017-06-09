"""Functions to calculate plasma parameters."""

from astropy import units as u

from astropy.units import (UnitConversionError, UnitsError, quantity_input,
                           Quantity)

from ..constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, ion_mass,
                         charge_state)

import numpy as np


"""Values should be returned as an Astropy Quantity in SI units.

If a quantity has several names, then the function name should be the
one that provides the most physical insight into what the quantity
represents.  For example, 'gyrofrequency' indicates gyration, while
Larmor frequency indicates that this frequency is somehow related to
someone named Larmor.  Similarly, using omega_ce as a function name
for this quantity will make the code less readable to people who are
unfamiliar with the notation.

The docstrings for plasma parameter methods should define these
quantities in ways that are understandable to students who are
taking their first course in plasma physics and yet are useful to
experienced plasma physicists.

Units that were named after a person should be lower case except at
the beginning of a sentence, even if their symbol is capitalized.

Unit conversions involving angles must be treated with care.  Angles
are dimensionless but do have units.  Angular velocity is often
given in units of radians per second, though dimensionally this is
equivalent to inverse seconds.  Astropy will treat radians
dimensionlessly when using the dimensionless_angles() equivalency,
but dimensionless_angles() does not account for multiplicative
factor of 2*pi that is used when converting between frequency (1 /
s) and angular velocity (rad / s).  A clear way to do this
conversion is to set up an equivalency between cycles/s and Hz:

  >>> f_ce = omega_ce.to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])

However, dimensionless_angles() does work when dividing a velocity
by an angular frequency to get a length scale:

  >>> d_i = (c/omega_pi).to(u.m, equivalencies=dimensionless_angles())

"""


def Alfven_speed(B, density, ion="p"):
    """Returns the Alfven speed.

    Parameters
    ----------
    B : Quantity
        The magnetic field magnitude in units convertible to tesla

    density: Quantity
        Either the ion number density in units convertible to 1 / m**3,
        or the mass density in units convertible to kg / m**3

    ion : string, optional
        Representation of the ion species.  If not given, then the ions
        are assumed to be protons.

    Returns
    -------
    V_A : Quantity with units of velocity
        The Alfven velocity of the plasma in units of meters per second.

    Raises
    ------
    TypeError
        The magnetic field and density arguments are not Quantities.

    UnitConversionError
        If the magnetic field or density is not in appropriate units.

    UserWarning
        If the Alfven velocity exceeds 10% of the speed of light.

    ValueError
        If the density is negative, or the ion mass or charge state
        cannot be found.

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

    Examples
    --------
    >>> import plasmapy
    >>> from astropy import units as u
    >>> B = 0.014*u.T
    >>> n = 5e19*u.m**-3
    >>> rho = n*(plasmapy.constants.m_p + plasmapy.constants.m_e)
    >>> ion = 'p'
    >>> Alfven_speed(B, n, ion)
    <Quantity 43173.871857213584 m / s>
    >>> Alfven_speed(B, rho, ion)
    <Quantity 43173.871857213584 m / s>
    >>> Alfven_speed(B, rho).to(u.cm/u.us)
    <Quantity 4.317387185721358 cm / us>

    """

    if not isinstance(B, Quantity) or not isinstance(density, Quantity):
        raise TypeError("The inputs to Alfven_speed must be Quantities.")

    if B.si.unit in ['1 / m3', 'kg / m3'] and density.si.unit in ['T']:
        B, density = density, B

    try:
        B = B.to(u.T)
    except Exception:
        raise UnitConversionError("The magnetic field in Alfven_speed cannot "
                                  "be converted to tesla.")

    density = density.si
    if np.any(density.value < 0):
        raise ValueError("The number or mass density in Alfven_speed cannot "
                         "be negative.")

    if density.unit == '1 / m3':

        try:
            m_i = ion_mass(ion)
            Z = charge_state(ion)
            if Z is None:
                Z = 1
        except Exception:
            raise ValueError("Invalid ion in Alfven_speed.")

        rho = density*m_i + Z*density*m_e

    elif density.unit == 'kg / m3':
        rho = density
    else:
        raise UnitsError("One input of Alfven_speed must have units of either "
                         "a number density or mass density.")

    if not np.all(np.isfinite(rho.value)) or not np.all(np.isreal(rho.value)) \
            or np.any(rho.value <= 0):
        raise ValueError("Invalid density in Alfven_speed.")

    if (not np.all(np.isfinite(B.value)) or not np.all(np.isreal(B.value))):
        raise ValueError("Invalid magnetic field strength in Alfven_speed")

    V_A = (np.abs(B)/np.sqrt(mu0*rho)).to(u.m/u.s)

    if np.any(V_A > c):
        raise UserWarning("Alfven speed is greater than the speed of light.")
    elif np.any(V_A > 0.1*c):
        raise UserWarning("Alfven speed is greater than 10% of the speed of "
                          "light.")

    beta = (V_A/c).value

    if np.any(beta > 1):
        raise UserWarning("The Alfven speed is roughly " +
                          str(round(beta, 2)) + " times the speed of light.")
    elif np.any(beta > 0.1):
        raise UserWarning("The Alfven speed is roughly " +
                          str(round(beta*100, 2)) + "% of the speed of light.")

    return V_A


def ion_sound_speed(T_i, ion='p', gamma=5/3):
    """Returns the ion sound speed.

    Parameters
    ----------
    T_i : Quantity
        Ion temperature

    ion : string, optional
        Representation of the ion species.  If not given, then the ions
        are assumed to be protons.

    gamma : float, optional
        The adiabatic index, which defaults to 5/3 as appropriate for a
        monatomic gas.  This quantity is the ratio of the heat capacity at
        constant pressure to the heat capacity at constant volume.

    Returns
    -------
    V_S : Quantity
        The ion sound speed in units of meters per second.

    Raises
    ------
    ValueError
        If the ion mass, adiabatic index, or temperature are invalid.

    UnitConversionError
        If the temperature is in incorrect units.

    UserWarning
        If the ion sound speed exceeds 10% of the speed of light

    Notes
    -----
    The ion sound speed is given by

    .. math::
    V_S = \sqrt{\frac{\gamma k_B T_i}{m_i}}

    Example
    -------
    >>> from astropy import units as u
    >>> ion_sound_speed(5e6*u.K)
    <Quantity 203155.10435273056 m / s>
    >>> ion_sound_speed(1*u.keV)
    <Quantity 309496.9076337826 m / s>

    """

    if gamma < 1:
        raise ValueError("The ratio of specific heats is less than one in "
                         "ion_sound_speed.")
    elif gamma is np.inf:
        return np.inf*u.m/u.s

    try:
        m_i = ion_mass(ion)
    except Exception:
        raise ValueError("Unable to find ion mass in ion_sound_speed.")

    try:
        T_i = T_i.to(u.K, equivalencies=u.temperature_energy())
    except Exception:
        raise u.UnitConversionError("The ion temperature in ion_sound_speed "
                                    "cannot be converted to kelvin.")

    if not np.all(np.isfinite(T_i.value)) or not np.all(np.isreal(T_i.value)) \
            or np.any(T_i.value <= 0):
        raise ValueError("Invalid temperature value in ion_sound_speed.")

    try:
        V_S = (np.sqrt(gamma*k_B*T_i/m_i)).to(u.m/u.s)
    except Exception:
        raise ValueError("Unable to find ion sound speed.")

    beta = (V_S/c).value

    if np.any(beta > 1):
        raise UserWarning("The ion sound speed is roughly " +
                          str(round(beta, 2)) + " times the speed of light.")
    elif np.any(beta > 0.1):
        raise UserWarning("The ion sound speed is roughly " +
                          str(round(beta*100, 2)) + "% of the speed of light.")

    return V_S


def electron_thermal_speed(T_e):
    """Returns the most probable speed for an electron within a
    Maxwellian distribution.

    Parameters
    ----------
    T_e : Quantity
        The electron temperature in either kelvin or electron-volts

    Returns
    -------
    V_Te : Quantity
        Electron thermal speed

    Notes
    -----
    The electron thermal speed is given by:

    .. math::
    V_{th,e} = \sqrt{\frac{2 k_B T_e}{m_e}}

    This function yields the most probable speed within a distribution
    function.  However, the definition of thermal velocity varies by
    the square root of two depending on whether or not this velocity
    absorbs that factor in the expression for a Maxwellian
    distribution.  In particular, the expression given in the NRL
    Plasma Formulary [1] is a square root of two smaller than the
    result from this function.

    Examples
    --------
    >>> electron_thermal_speed(5*u.eV)
    <Quantity 1326205.1454609886 m / s>
    >>> electron_thermal_speed(1e6*u.K)
    <Quantity 5505694.743141063 m / s>

    """

    try:
        T_e = T_e.to(u.K, equivalencies=u.temperature_energy())
    except Exception:
        raise u.UnitConversionError("The electron temperature in "
                                    "electron_thermal_speed "
                                    "cannot be converted to kelvin.")

    if np.any(not 0 <= T_e.value <= np.inf):
        raise ValueError("Invalid temperature in electron_thermal_speed.")

    try:
        V_Te = (np.sqrt(2*k_B*T_e/m_e)).to(u.m/u.s)
    except Exception:
        raise ValueError("Cannot find electron thermal speed")

    beta = (V_Te/c).value

    if np.any(beta > 1):
        raise UserWarning("The electron thermal speed is roughly " +
                          str(round(beta, 2)) + " times the speed of light.")
    elif np.any(beta > 0.1):
        raise UserWarning("The electron thermal is roughly " +
                          str(round(beta*100, 2)) + "% of the speed of light.")

    return V_Te


def ion_thermal_speed(T_i, ion='p'):
    """Returns the most probable speed for an ion within a Maxwellian
    distribution.

    Parameters
    ----------
    T_i : Quantity
        The ion temperature in either kelvin or electron-volts

    ion : string
        Symbol representing the ion species, defaulting to protons

    Returns
    -------
    V_Ti : Quantity
        Ion thermal speed

    Notes
    -----
    The electron thermal speed is given by:

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
    >>> ion_thermal_speed(5*u.eV)
    <Quantity 30949.690763378258 m / s>
    >>> ion_thermal_speed(1e6*u.K, ion='p')
    <Quantity 128486.56960876317 m / s>

    """

    try:
        T_i = T_i.to(u.K, equivalencies=u.temperature_energy())
        m_i = ion_mass(ion)
        V_Ti = (np.sqrt(2*k_B*T_i/m_i)).to(u.m/u.s)
    except Exception:
        raise ValueError("Unable to find ion thermal speed")

    beta = (V_Ti/c).value

    if np.any(beta > 1):
        raise UserWarning("The electron thermal speed is roughly " +
                          str(round(beta, 2)) + " times the speed of light.")
    elif np.any(beta > 0.1):
        raise UserWarning("The electron thermal is roughly " +
                          str(round(beta*100, 2)) + "% of the speed of light.")

    return V_Ti


@u.quantity_input
def electron_gyrofrequency(B: u.T):
    """Calculate the electron gyrofrequency in units of radians per second.

    Parameters
    ----------
    B: Quantity
        Magnetic field strength

    Returns
    -------
    omega_ce: Quantity
        Electron gyrofrequency in radians per second.

    Notes
    -----
    The electron gyrofrequency is the angular frequency of electrons
    gyration around magnetic field lines and is given by:

    .. math::
    omega_{ce} = \frac{e B}{m_e}

    The electron gyrofrequency is also known as the electron cyclotron
    frequency or the electron Larmor frequency.

    The recommended way to convert from angular frequency to frequency
    is to use an equivalency between cycles per second and Hertz, as
    Astropy's dimensionles_angles() equivalency does not account for
    the factor of 2*pi needed during this conversion.  The
    dimensionless_angles() equivalency is appropriate when dividing a
    velocity by an angular frequency to get a length scale.

    Examples
    --------
    >>> from astropy import units as u
    >>> from numpy import pi
    >>> omega_ce = electron_gyrofrequency(0.1*u.T)
    >>> print(omega_ce)
    17588200878.472023 rad / s
    >>> f_ce = omega_ce.to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])
    >>> print(f_ce)
    2799249109.9020386 Hz

    """

    try:
        omega_ce = u.rad*(e*B/m_e).to(1/u.s)
    except Exception:
        raise ValueError("Unable to find electron gyrofrequency")

    return omega_ce


def ion_gyrofrequency(B, ion='p', Z=None):
    """Calculate the ion gyrofrequency in units of radians per second.

    Parameters
    ----------
    B: Quantity
        The magnetic field magnitude in units convertible to tesla

    ion : string, optional
        Representation of the ion species.  If not given, then the ions
        are assumed to be protons.

    Z : integer, optional
        The charge state.  If not given, then the ion is assumed to be
        singly charged.

    Returns
    -------
    omega_ci: Quantity
        The ion gyrofrequency in units of radians per second

    Notes
    -----
    The ion gyrofrequency is the angular frequency of ion gyration
    around magnetic field lines and is given by:

    .. math::
    omega_{ci} = \frac{e B}{m_i}

    The ion gyrofrequency is also known as the ion cyclotron frequency
    or the ion Larmor frequency.

    The recommended way to convert from angular frequency to frequency
    is to use an equivalency between cycles per second and Hertz, as
    Astropy's dimensionles_angles() equivalency does not account for
    the factor of 2*pi needed during this conversion.  The
    dimensionless_angles() equivalency is appropriate when dividing a
    velocity by an angular frequency to get a length scale.

    Example
    -------
    >>> from astropy import units as u
    >>> ion_gyrofrequency(0.01*u.T)
    <Quantity 152451.87138666757 Hz>
    >>> ion_gyrofrequency(0.01*u.T, 'p')
    <Quantity 152451.87138666757 Hz>
    >>> ion_gyrofrequency(0.01*u.T, ion='T', Z=1)


    """

    try:
        B = np.abs(B.to(u.T))
    except Exception:
        raise UnitsError("The first argument of ion_gyrofrequency must "
                         "be a quantity with units convertible to tesla.")

    if Z is None:
        try:
            Z = charge_state(ion)
            if Z is None:
                Z = 1
        except:
            raise ValueError("Unable to find charge state in "
                             "ion_gyrofrequency")
    elif not isinstance(Z, int):
        raise TypeError("In ion_gyrofrequency, the charge state (Z) must be "
                        "an integer.")

    try:
        m_i = ion_mass(ion, Z=Z)
    except Exception:
        raise ValueError("Unable to get ion mass in ion_gyrofrequency")

    try:
        omega_ci = u.rad * (Z*e*B/m_i).to(1/u.s)
    except Exception:
        raise ValueError("Unable to find ion gyrofrequency")

    return omega_ci


def electron_gyroradius(B, Vperp_or_Te):
    """Returns the radius of gyration for an electron in a uniform
    magnetic field.

    Parameters
    ----------
    B : Quantity
        The magnetic field magnitude in units convertible to tesla

    Vperp_or_Te : Quantity
        Either the component of electron velocity that is
        perpendicular to the magnetic field in units convertible to
        meters per second, or the electron temperature in units
        convertible to kelvin.

    Returns
    -------
    r_L : Quantity

        The electron gyroradius in units of meters.  This quantity
        corresponds to either the perpendicular component of electron
        velocity, or the most probable speed for an electron within a
        Maxwellian distribution for the electron temperature.

    Notes
    -----
    The electron gyroradius is also known as the Larmor radius for
    electrons and is given by:

    .. math::
    r_L = \frac{V_{perp}}{omega_{ce}}

    """

    if not isinstance(B, Quantity):
        raise TypeError("The first input to electron_gyroradius is not a "
                        "quantity.")

    if not isinstance(Vperp_or_Te, Quantity):
        raise TypeError("The second input to electron_gyroradius is not a "
                        "quantity.")

    Vperp_or_Te = Vperp_or_Te.si

    if Vperp_or_Te.si.unit == 'T':
        B, Vperp_or_Te = Vperp_or_Te, B

    try:
        B = B.to(u.T)
    except:
        raise UnitConversionError("electron_gyroradius requires a quantity "
                                  "with units convertible to tesla as an "
                                  "argument.")

    if Vperp_or_Te.unit in ['J', 'K']:
        T_e = Vperp_or_Te.to(u.K, equivalencies=u.temperature_energy())
        if T_e < 0*u.K:
            raise ValueError("An argument to electron_gyroradius corresponds "
                             "to a negative temperature.")
        Vperp = electron_thermal_speed(T_e)
    elif Vperp_or_Te.unit == 'm / s':
        Vperp = Vperp_or_Te
    else:
        raise UnitConversionError("electron_gyroradius requires a quantity "
                                  "with units convertible to kelvin or to "
                                  "meters per second as an argument.")

    omega_ce = electron_gyrofrequency(B)
    r_L = (Vperp/omega_ce).to(u.m, equivalencies=u.dimensionless_angles())

    return r_L


def ion_gyroradius(B, Vperp_or_Ti, ion='p'):
    """Returns the ion gyroradius.

    Parameters
    ----------
    B: Quantity
        The magnetic field magnitude in units convertible to tesla

    Vperp_or_Ti: Quantity
        The component of ion velocity that is perpendicular to the
        magnetic field in units convertible to meters per second, or the
        ion temperature in units convertible to kelvin.

    ion : string, optional
        Representation of the ion species.  If not given, then the ions
        are assumed to be protons.

    Notes
    -----

    The ion gyroradius is also known as the ion Larmor radius and is
    given by

    .. math::
    r_L = \frac{V_{perp}}{omega_{ci}}

    """

    if not isinstance(B, Quantity):
        raise TypeError("The first input to ion_gyroradius is not a "
                        "quantity.")

    if not isinstance(Vperp_or_Ti, Quantity):
        raise TypeError("The second input to ion_gyroradius is not a "
                        "quantity.")

    Vperp_or_Ti = Vperp_or_Ti.si

    if Vperp_or_Ti.si.unit == 'T':
        B, Vperp_or_Ti = Vperp_or_Ti, B

    try:
        B = B.to(u.T)
    except:
        raise UnitConversionError("ion_gyroradius requires a quantity "
                                  "with units convertible to tesla as an "
                                  "argument.")

    if Vperp_or_Ti.unit in ['J', 'K']:
        T_i = Vperp_or_Ti.to(u.K, equivalencies=u.temperature_energy())
        if T_i < 0*u.K:
            raise ValueError("An argument to ion_gyroradius corresponds "
                             "to a negative temperature.")
        Vperp = ion_thermal_speed(T_i, ion)
    elif Vperp_or_Ti.unit == 'm / s':
        Vperp = Vperp_or_Ti
    else:
        raise UnitConversionError("ion_gyroradius requires a quantity "
                                  "with units convertible to kelvin or to "
                                  "meters per second as an argument.")

    omega_ci = ion_gyrofrequency(B, ion)
    r_L = (Vperp/omega_ci).to(u.m, equivalencies=u.dimensionless_angles())

    return r_L


@quantity_input
def electron_plasma_frequency(n_e: u.m**-3):
    """Calculates the electron plasma frequency.

    Parameters
    ----------
    n_e: Quantity
        Electron number density

    Returns
    -------
    omega_pe: Quantity
        Electron plasma frequency in radians per second

    Notes
    -----
    The electron plasma frequency is

    .. math::
    \omega_{pe} = e \sqrt{\frac{n_e}{\epsilon_0 m_e}}

    At present, astropy.units does not allow direct conversions from
    radians/second for angular frequency to 1/second or Hz for
    frequency.  The dimensionless_angles equivalency allows that
    conversion, but does not account for the factor of 2*pi.  The
    alternatives are to convert to cycle/second or to do the
    conversion manually, as shown in the examples.

    Example
    -------
    >>>

    """
    try:
        omega_pe = (u.rad*e*np.sqrt(n_e/(eps0*m_e))).to(u.rad/u.s)
    except Exception:
        raise ValueError("Unable to find electron plasma frequency.")

    return omega_pe


def ion_plasma_frequency(n_i, Z=None, ion='p'):
    """Calculates the ion plasma frequency.

    Parameters
    ----------
    n_i : Quantity
        Ion number density
    Z : integer
        The ionization state (e.g., Z=1 for singly ionized)
    ion : string, optional
        Representation of the ion species.  If not given, then the ions
        are assumed to be protons.

    Returns
    -------
    omega_pi : Quantity
        The ion plasma frequency in radians per second.

    Notes
    -----
    The ion plasma frequency is

    .. math::
    \omega_{pi} = Z e \sqrt{\frac{n_i}{\epsilon_0 m_i}}

    At present, astropy.units does not allow direct conversions from
    radians/second for angular frequency to 1/second or Hz for
    frequency.  The dimensionless_angles equivalency allows that
    conversion, but does not account for the factor of 2*pi.  The
    alternatives are to convert to cycle/second or to do the
    conversion manually, as shown in the examples.

    Example
    -------
    >>>

    """

    try:
        m_i = ion_mass(ion, Z=Z)
    except Exception:
        raise ValueError("Unable to get mass of ion in ion_gyrofrequency")

    if Z is None:
        try:
            Z = charge_state(ion)
        except Exception:
            raise ValueError("Unable to get charge state to calculate ion "
                             "plasma frequency.")

    try:
        omega_pi = u.rad*Z*e*np.sqrt(n_i/(eps0*m_e))
    except Exception:
        raise ValueError("Unable to find ion plasma frequency.")

    return omega_pi.si


@u.quantity_input
def Debye_length(T_e: u.K, n_e: u.m**-3):  # Add equivalency related to T in eV
    """Calculate the Debye length.

    Parameters
    ----------
    T_e: Quantity
        Electron temperature

    n_e: Quantity
        Electron number density

    Returns
    -------
    lambda_D: Quantity
        The Debye length

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

    References
    ----------
    [1] Declan Diver plasma formulary.......

    """

    try:
        lambda_D = ((eps0*k_B*T_e / (n_e * e**2))**0.5).to(u.m)
    except Exception:
        raise ValueError("Unable to find Debye length.")

    return lambda_D


@u.quantity_input
def Debye_number(T_e: u.K, n_e: u.m**-3):
    """Returns the Debye number.

    Parameters
    ----------
    T_e: Quantity
        Electron temperature

    n_e: Quantity
        Electron number density

    Returns
    -------
    N_D : Quantity
        Number of electrons within a sphere with a radius of the Debye length

    Notes
    -----
    The Debye number is the number of electrons contained within a sphere with
    a radius of a Debye length and is given by

    .. math::
    N_D = \frac{4}{3}n_e\lambda_D^3

    The Debye number is also known as the plasma parameter.

    Collective behavior requires a Debye number significantly larger than one.

    See also
    --------
    Debye_length

    Example
    -------
    >>> Debye_number(5e6*u.K, 5e9*u.cm**-3)
    <Quantity 69282817.49483868>

    """

    try:
        lambda_D = Debye_length(T_e, n_e)
        N_D = (4/3)*n_e*lambda_D**3
    except Exception:
        raise ValueError("Unable to find Debye number.")

    return N_D.to(u.dimensionless_unscaled)


def ion_inertial_length(n_i, ion='p', Z=None):
    """Calculate the ion inertial length,

    Parameters
    ----------
    n_i : Quantity
        Ion number density
    ion : string, optional
        Representation of the ion species.  If not given, then the ions
        are assumed to be protons.

    Returns
    -------


    Notes
    -----
    The ion inertial length is also known as an ion skin depth and is
    given by:

    .. math::
    d_i = \frac{c}{\omega_{pi}}

    Example
    -------
    >>> from astropy import units as u
    >>> ion_plasma_frequency(5e19*u.m**-3, 1, 'p')
    <Quantity 398911478582.3019 rad / s>
    >>> ion_plasma_frequency(5e19*u.m**-3)  # assumes Z=1 and ion='p'
    <Quantity 398911478582.3019 rad / s>
    >>> from plasmapy.constants import c
    >>> omega_pi = ion_plasma_frequency(1e15*u.m**-3)
    >>> c/omega_pi

    """

    try:
        Z = charge_state(ion)
    except Exception:
        raise ValueError("Unable to find charge state in ion_inertial_length.")

    try:
        omega_pi = ion_plasma_frequency(n_i, Z=1, ion=ion)
        d_i = (c/omega_pi).to(u.m, equivalencies=u.dimensionless_angles())
    except Exception:
        raise ValueError("Unable to find ion inertial length.")

    return d_i


@u.quantity_input
def electron_inertial_length(n_e: u.m**-3):
    """Returns the electron inertial length.

    Parameters
    ----------
    n_e : Quantity
        Electron number density

    Returns
    -------
    d_e : Quantity
        Electron inertial length

    Notes
    -----
    The electron inertial length is also known as an electron skin depth and
    is given by:

    .. math::
    d_e = \frac{c}{\omega_{pe}}

    Example
    -------
    >>>

    """

    try:
        omega_pe = electron_plasma_frequency(n_e)
        d_e = (c/omega_pe).to(u.m, equivalencies=u.dimensionless_angles())
    except Exception:
        raise ValueError("Unable to find electron inertial length.")

    return d_e


@u.quantity_input
def magnetic_pressure(B: u.T):
    """Calculate the magnetic pressure.

    Parameters
    ----------
    B: Quantity
        The magnetic field in units convertible to telsa

    Returns
    -------
    p_B: Quantity
        The magnetic pressure in units in pascals (newtons per square meter)

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
    magnetic_energy_density : returns an equivalent quantity, except in
        units of joules per cubic meter.

    Example
    -------
    >>> from astropy import units as u
    >>> magnetic_pressure(0.1*u.T).to(u.Pa)
    <Quantity 3978.873577297384 Pa>
    """

    try:
        p_B = (B**2/(2*mu0)).to(u.N/u.m**2)
    except Exception:
        raise ValueError("Unable to find magnetic pressure.")

    return p_B


@u.quantity_input
def magnetic_energy_density(B: u.T):
    """Calculate the magnetic energy density.

    Parameters
    ----------
    B: Quantity
        The magnetic field in units convertible to tesla

    Returns
    -------
    E_B: Quantity
        The magnetic energy density in units of joules per cubic meter

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
    magnetic_pressure : returns an equivalent quantity, except in units
        of pascals.

    Example
    -------
    >>> from astropy import units as u
    >>> magnetic_energy_density(0.1*u.T)
    <Quantity 3978.873577297384 J / m3>

    """

    try:
        E_B = (B**2/(2*mu0)).to(u.J/u.m**3)
    except Exception:
        raise ValueError("Unable to find magnetic energy density.")

    return E_B
