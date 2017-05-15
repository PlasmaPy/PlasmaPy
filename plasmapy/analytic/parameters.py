"""Functions to calculate plasma parameters."""

from astropy import units as u

from astropy.units import UnitConversionError, quantity_input

from ..constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, ion_mass)

import numpy as np


"""
Values should be returned as an Astropy Quantity in SI units.

If a quantity has several names, then the function name should be
the one that provides the most physical insight into what the
quantity represents.  For example, 'gyrofrequency' indicates
gyration, while Larmor frequency indicates that this frequency was
discovered by Larmor.

The docstrings for plasma parameter methods should define these
quantities in ways that are understandable to students who are
taking their first course in plasma physics and yet are useful to
experienced plasma physicists.

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

def charge_state(ion):
    """Temporary: gets Z from a string representing an ion."""
    print("Need to update charge_state")
    return 1


def Alfven_speed(B, density, ion="p"):
    """Returns the Alfven speed.

    Parameters
    ----------
    B : Quantity
        The magnetic field in units such as Tesla or Gauss
    density: Quantity
        Either the ion number density in units such as 1 / m**3, or the
        mass density in units such as kg / m**3
    ion : string
        Symbol representing the ion species, defaulting to protons

    Returns
    -------
    V_A : Quantity with units of velocity
        The Alfven velocity of the plasma in units of meters per second.

    Raises
    ------
    UserWarning
        If the Alfven velocity exceeds the speed of light

    Notes
    -----
    The Alfven velocity :math:`V_A` is the typical propagation speed
    of magnetic disturbances in a plasma, and is given by:

    .. math::
    V_A = \frac{B}{\sqrt{\mu_0\rho}}

    where the mass density is :math:`\rho \approx n_i m_i`.

    Examples
    --------

    >>> from astropy import units as u
    >>> import plasmapy
    >>> B = 140 * u.G
    >>> n = 5e19 * u.m**-3
    >>> rho = n * plasmapy.constants.m_i
    >>> element = 'H'
    >>> V_A1 = Alfven_speed(B, n, element)
    >>> print(V_A1)
    43185.62686969444 m / s
    >>> V_A2 = Alfven_speed(B, rho, element)
    >>> print(V_A2)
    43185.62686969444 m / s
    >>> print(V_A2.cgs)
    4318562.686969444 cm / s
    >>> print(V_A2.to(u.km/u.s))
    43.18562686969444 km / s

    """

    try:
        if B.si.unit in ['1 / m3', 'kg / m3'] and density.si.unit in ['T']:
            B, density = density, B
    except:
        pass

    try:
        B = B.to(u.T)
    except:
        raise UnitConversionError("The magnetic field in Alfven_speed cannot "
                                  "be converted to Tesla")

    try:
        density = density.si
    except:
        raise UnitConversionError('The density in Alfven_speed')

    if density.unit not in ['1 / m3', 'kg / m3']:
        raise UnitsError("One input of Alfven_speed must have units of either "
                         "a number density or mass density.")

    if density.unit == '1 / m3':
        m_i = ion_mass(ion)
        Z = charge_state(ion)
        rho = density*m_i + Z*density*m_e
    elif density.unit == 'kg / m3':
        rho = density

    V_A = (np.abs(B)/np.sqrt(mu0*rho)).to(u.m/u.s)

    if V_A > c:
        raise UserWarning("Alfven speed is greater than speed of light")

    return V_A


def ion_sound_speed(T_i, ion='p', gamma=5/3):
    """Returns the ion sound speed.

    Parameters
    ----------
    T_i : Quantity
        Ion temperature.
    ion : string, optional
        Representation of the ion species (assumed to be protons if not given)
    gamma : float, optional
        Heat capacity ratio (also known as the adiabatic index, ratio of
        specific heats, Poisson constant, or isentropic expansion parameter).
        Defaults to 5/3 which is appropriate for a monatomic gas.

    Returns
    -------
    V_S : Quantity
        The ion sound speed.

    Example
    -------
    >>> from astropy import units as u
    >>> ion_sound_speed(5e6*u.K)
    <Quantity 203155.10435273056 m / s>
    >>> ion_sound_speed(1*u.keV)
    <Quantity 309496.9076337826 m / s>

    """

    m_i = ion_mass(ion)
    T_i = T_i.to(u.K, equivalencies=u.temperature_energy())
    V_S = (np.sqrt(gamma*k_B*T_i/m_i)).to(u.m/u.s)
    return V_S


#@u.quantity_input
def electron_thermal_speed(T_e):
    """Returns the most probable speed for an electron within a
    Maxwellian distribution.

    Parameters
    ----------
    T_e : Quantity
        The electron temperature.

    Returns
    -------
    V_Te : Quantity
        Electron thermal speed

    Notes
    -----
    The electron thermal speed is given by:
    
    MATH MATH MATH MATH MATH

    This function yields the most probable speed within a distribution
    function.  However, the definition of thermal velocity varies by
    the square root of two depending on whether or not this velocity
    absorbs that factor in the expression for a Maxwellian
    distribution.  In particular, the expression given in the NRL
    Plasma Formulary [1] is a square root of two smaller than the
    result from this function.

    """

    try:
        T_e = T_e.to(u.K, equivalencies=u.temperature_energy())
    except:
        raise UnitConversionError("Cannot convert input of "
                                  "electron_thermal_velocity to Kelvin")

    V_Te = (np.sqrt(2*k_B*T_e/m_e)).to(u.m/u.s)

    return V_Te


def ion_thermal_speed(T_i, ion='p'):
    """Returns the most probable speed for an ion within a Maxwellian
    distribution.

    Parameters
    ----------
    T_i : Quantity
        The ion temperature

    ion : string
        A representation of the ion species.

    Returns
    -------
    V_Ti : Quantity
        Ion thermal speed

    Notes
    -----
    The electron thermal speed is given by:
    
    MATH MATH MATH MATH MATH

    This function yields the most probable speed within a distribution
    function.  However, the definition of thermal velocity varies by
    the square root of two depending on whether or not this velocity
    absorbs that factor in the expression for a Maxwellian
    distribution.  In particular, the expression given in the NRL
    Plasma Formulary [1] is a square root of two smaller than the
    result from this function.

    Examples
    --------

    >>>

    """
    try:
        T_i = T_i.to(u.K, equivalencies=u.temperature_energy())
    except:
        raise UnitConversionError("Cannot convert input of "
                                  "ion_thermal_velocity to Kelvin")

    V_Ti = (np.sqrt(2*k_B*T_i/m_i)).to(u.m/u.s)

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
    gyration around magnetic field lines.

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

    omega_ce = u.rad*(e*B/m_e).to(1/u.s)

    return omega_ce


@u.quantity_input
def ion_gyrofrequency(B: u.T, ion='p'):
    """Calculate the ion gyrofrequency in units of radians per second.

    Parameters
    ----------
    B: Quantity
        Magnetic field strength
    ion: string
        Representation of the ion

    Returns
    -------
    omega_ci: Quantity
        The ion gyrofrequency in units of radians per second

    Notes
    -----
    The ion gyrofrequency is the angular frequency of ion gyration
    around magnetic field lines.

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
    >>>

    """

    m_i = ion_mass(ion)

    try:
        B = B.to(u.T)
    except:
        raise UnitsError("Argument 'B' to function 'ion_gyrofrequency' must "
                         "be in units convertible to 'T'")

    omega_ci = u.rad * (e*B/m_i).to(1/u.s)

    return omega_ci


#@u.quantity_input
def electron_gyroradius(B, Te_or_Vperp):
    """Returns the radius of gyration for an electron in a uniform
    magnetic field.

    Parameters
    ----------
    B : Quantity
        Magnetic field strength
    Te_or_Vperp : Quantity
        Either the electron temperature or the perpendicular velocity.

    Returns
    -------
    r_L : Quantity
        Electron gyroradius

    Notes
    -----
    The electron gyroradius is also known as the Larmor radius for
    electrons.

    """
    
    try:
        if Te_or_Vperp.si.unit == 'T' and (B.si.unit in ['m / s', 'K'] or 
                                           B.unit[-2:] == 'eV'):
            B, Te_or_Vperp = Te_or_Vperp, B
    except:
        pass

    B = B.to(u.T)
    
    try:
        T_e = Te_or_Vperp.to(u.K, equivalencies=u.temperature_energy())
        Vperp = electron_thermal_speed(T_e).to(u.m/u.s)
    except:
        try:
            Vperp = Te_or_Vperp.to(u.m/u.s)
        except:
            raise UnitConversionError("Incorrect inputs for electron_gyroradius")

#    if Te_Te_or_Vperp.unit[-2:] == 'eV':
#        T_e = Te_or_Vperp.to(u.K, equivalencies=u.temperature_energy())
#        V_perp = electron_thermal_speed(T_e)
#    elif 

    omega_ce = electron_gyrofrequency(B)
    r_L = (Vperp/omega_ce).to(u.m, equivalencies=u.dimensionless_angles())
    
    return r_L

def ion_gyroradius(B=None, V=None, T_i=None, ion='p'):
    """Returns the ion gyroradius.

    Parameters
    ----------
    B: Quantity
        Magnetic field strength
    V: Quantity
        Ion velocity
    ion : str, optional



    Notes
    -----


    The ion gyroradius is also known as the Larmor radius for ions.

    """
    return 0


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

    omega_pe = (u.rad*e*np.sqrt(n_e/(eps0*m_e))).to(u.rad/u.s)

    return omega_pe


def ion_plasma_frequency(n_i, Z=None, ion='p'):
    """Calculates the ion plasma frequency.

    Parameters
    ----------
    n_i : Quantity
        Ion number density
    ion : string
        A representation of the ion species
    Z : integer
        The ionization state (e.g., Z=1 for singly ionized)

    Returns
    -------
    omega_pi : Quantity
        The ion plasma frequency in radians per second.

    Notes
    -----
    The ion plasma frequency is

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

    m_i = ion_mass(ion)
    omega_pi = u.rad*Z*e*np.sqrt(n_i/(eps0*m_e))
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

    MATH MATH MATH MATH MATH

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

    lambda_D = ((eps0*k_B*T_e / (n_e * e**2))**0.5).to(u.m)
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
    a radius of a Debye length.

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

    lambda_D = Debye_length(T_e, n_e)
    N_D = (4/3)*n_e*lambda_D**3
    return N_D.to(u.dimensionless_unscaled)


def ion_inertial_length(n_i, ion='p', Z=None):
    """Calculate the ion inertial length,

    Parameters
    ----------
    n_i : Quantity
        Ion number density

    Returns
    -------


    Notes
    -----
    The ion inertial length is also known as an ion skin depth.

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

    if Z is not None:  # need to think about this more...
        Z = charge_state(ion)

    omega_pi = ion_plasma_frequency(n_i, Z=1, ion=ion)
    d_i = (c/omega_pi).to(u.m, equivalencies=u.dimensionless_angles())
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
    The electron inertial length is also known as an electron skin depth.


    Example
    -------
    >>>

    """

    omega_pe = electron_plasma_frequency(n_e)
    d_e = (c/omega_pe).to(u.m, equivalencies=u.dimensionless_angles())
    return d_e



@u.quantity_input
def magnetic_pressure(B: u.T):
    """Calculate the magnetic pressure.

    Parameters
    ----------
    B: Quantity
        The magnetic field in units such as Tesla or Gauss

    Returns
    -------
    p_B: Quantity
        The magnetic pressure in units Pascals (Newtons per square meter)

    Notes
    -----

    See also
    --------
    magnetic_energy_density : returns an equivalent quantity, except with
        units converted to Joules per cubic meter.

    Example
    -------
    >>> magnetic_pressure(0.1*u.T).to(u.Pa)
    <Quantity 3978.873577297384 Pa>
    """

    p_B = (B**2/(2*mu0)).to(u.N/u.m**2)

    return p_B


@u.quantity_input
def magnetic_energy_density(B: u.T):
    """Calculate the magnetic energy density.

    Parameters
    ----------
    B: Quantity
        The magnetic field in units such as Tesla or Gauss

    Returns
    -------
    E_B: Quantity
        The magnetic energy density in units of Joules per cubic meter

    Notes
    -----
    The expressions for magnetic pressure
    This function is very similar to magnetic_pressure, except for the

    See also
    --------
    magnetic_pressure : returns an equivalent quantity, except with units
        converted to Pascals (Newtons per square meter).

    Example
    -------
    >>> magnetic_energy_density(0.1*u.T)
    <Quantity 3978.873577297384 J / m3>

    """

    E_B = (B**2/(2*mu0)).to(u.J/u.m**3)

    return E_B
