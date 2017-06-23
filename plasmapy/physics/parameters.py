"""Functions to calculate plasma parameters."""

from astropy import units as u

from astropy.units import (UnitConversionError, UnitsError, quantity_input,
                           Quantity)

from ..constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, ion_mass,
                         charge_state)

import numpy as np

from ..utils import _check_quantity, _check_relativistic

r"""Values should be returned as an Astropy Quantity in SI units.

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
    r"""Returns the Alfven speed.

    Parameters
    ----------
    B : Quantity
        The magnetic field magnitude in units convertible to tesla.

    density: Quantity
        Either the ion number density in units convertible to 1 / m**3,
        or the mass density in units convertible to kg / m**3.

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

    if B.si.unit in [u.m**-3, u.kg/u.m**3] and density.si.unit in [u.T]:
        B, density = density, B

    _check_quantity(B, 'B', 'Alfven_speed', u.T)
    _check_quantity(density, 'density', 'Alfven_speed', [u.m**-3, u.kg/u.m**3],
                    can_be_negative=False)

    B = B.to(u.T)
    density = density.si

    if density.unit == u.m**-3:
        try:
            m_i = ion_mass(ion)
            Z = charge_state(ion)
            if Z is None:
                Z = 1
        except Exception:
            raise ValueError("Invalid ion in Alfven_speed.")
        rho = density*m_i + Z*density*m_e

    elif density.unit == u.kg/u.m**3:
        rho = density
    else:
        raise UnitsError("One input of Alfven_speed must have units of either "
                         "a number density or mass density.")

    V_A = (np.abs(B)/np.sqrt(mu0*rho)).to(u.m/u.s)

    _check_relativistic(V_A, 'Alfven_speed')

    return V_A


def ion_sound_speed(T_i=None, T_e=None, ion='p', gamma_e=1, gamma_i=3):
    r"""Returns the ion sound speed for an electron-ion plasma.

    Parameters
    ----------
    T_i : Quantity, optional if T_e is given
        Ion temperature in units of temperature or energy per
        particle.  If this is not given, then the ion temperature is
        assumed to be equal to the electron temperature.

    T_e : Quantity, optional if T_i is given
        Electron temperature in units of temperature or energy per
        particle.  If this is not given, then the electron temperature
        is assumed to be equal to the ion temperature.

    ion : string, optional
        Representation of the ion species.  If not given, then the ions
        are assumed to be protons.  If the charge state is not given, then
        ions are assumed to be singly charged.

    gamma_i : float or int
        The adiabatic index for ions, which defaults to 3.  This value
        assumes that ion motion has only one degree of freedom, namely
        along magnetic field lines.

    gamma_e : float or int
        The adiabatic index for electrons, which defaults to 1.  This
        value assumes that the electrons are able to equalize their
        temperature rapidly enough that the electrons are effectively
        isothermal.

    Returns
    -------
    V_S : Quantity
        The ion sound speed in units of meters per second.

    Raises
    ------
    TypeError
        If T_i or T_e are not quantities.

    ValueError
        If the ion mass, adiabatic index, or temperature are invalid.

    UnitConversionError
        If the temperature is in incorrect units.

    UserWarning
        If the ion sound speed exceeds 10% of the speed of light.

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
    dispersive (frequency independent).

    When the electron temperature is much greater than the ion
    temperature, the ion sound velocity reduces to
    :math:`\sqrt{\gamma_e k_B T_e / m_i}`.  Ion acoustic waves can
    therefore occur even when the ion temperature is zero.

    Example
    -------
    >>> from astropy import units as u
    >>> ion_sound_speed(5e6*u.K, ion='p', gamma=5/3)
    <Quantity 262272.11195171846 m / s>
    >>> ion_sound_speed(1*u.keV, ion='D+', gamma=5/3)
    <Quantity 282600.8484562855 m / s>

    """

    try:
        m_i = ion_mass(ion)
    except Exception:
        raise ValueError("Unable to find ion mass in ion_sound_speed.")

    try:
        Z = charge_state(ion)
        if Z is None:
            Z = 1
    except Exception:
        raise ValueError("Unable to find charge state in ion_sound_speed")

    if not isinstance(gamma_e, (float, int)):
        raise TypeError("The adiabatic index for electrons (gamma_e) must be "
                        "a float or int in ion_sound_speed")
    if not isinstance(gamma_i, (float, int)):
        raise TypeError("The adiabatic index for ions (gamma_i) must be "
                        "a float or int in ion_sound_speed")

    if not 1 <= gamma_e <= np.inf:
        raise ValueError("The adiabatic index for electrons must be between "
                         "one and infinity")
    if not 1 <= gamma_i <= np.inf:
        raise ValueError("The adiabatic index for ions must be between "
                         "one and infinity")

    if (T_i is None) ^ (T_e is None):
        if T_i is None:
            T_i = T_e
        else:
            T_e = T_i

    _check_quantity(T_i, 'T_i', 'ion_sound_speed', u.K, can_be_negative=False)
    _check_quantity(T_e, 'T_e', 'ion_sound_speed', u.K, can_be_negative=False)

    T_i = T_i.to(u.K, equivalencies=u.temperature_energy())
    T_e = T_e.to(u.K, equivalencies=u.temperature_energy())

    V_S = np.sqrt((gamma_e*Z*k_B*T_e + gamma_i*k_B*T_i)/m_i).to(u.m/u.s)

    _check_relativistic(V_S, 'ion_sound_speed')

    return V_S


def electron_thermal_speed(T_e):
    r"""Returns the most probable speed for an electron within a
    Maxwellian distribution.

    Parameters
    ----------
    T_e : Quantity
        The electron temperature in either kelvin or energy per particle

    Returns
    -------
    V_Te : Quantity
        Electron thermal speed

    Raises
    ------
    TypeError
        The electron temperature is not a Quantity

    UnitConversionError
        If the electron temperature is not in appropriate units

    ValueError
        The electron temperature is invalid

    UserWarning
        If the electron thermal speed exceeds 10% of the speed of light.

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
    >>> from astropy import units as u
    >>> electron_thermal_speed(5*u.eV)
    <Quantity 1326205.1454609886 m / s>
    >>> electron_thermal_speed(1e6*u.K)
    <Quantity 5505694.743141063 m / s>

    """

    if not isinstance(T_e, Quantity):
        raise TypeError("The input to electron_thermal_speed must "
                        "be a Quantity.")

    try:
        T_e = T_e.to(u.K, equivalencies=u.temperature_energy())
    except Exception:
        raise u.UnitConversionError("The electron temperature in "
                                    "electron_thermal_speed "
                                    "cannot be converted to kelvin.")

    if np.any(not 0 <= T_e.value <= np.inf):
        raise ValueError("Invalid temperature in electron_thermal_speed.")

    V_Te = (np.sqrt(2*k_B*T_e/m_e)).to(u.m/u.s)

    beta = (V_Te/c).value

    if np.any(beta > 1):
        raise UserWarning("The electron thermal speed is roughly " +
                          str(round(beta, 2)) + " times the speed of light.")
    elif np.any(beta > 0.1):
        raise UserWarning("The electron thermal is roughly " +
                          str(round(beta*100, 2)) + "% of the speed of light.")

    return V_Te


def ion_thermal_speed(T_i, ion='p'):
    r"""Returns the most probable speed for an ion within a Maxwellian
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
    >>> from astropy import units as u
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


def electron_gyrofrequency(B):
    r"""Calculate the electron gyrofrequency in units of radians per second.

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

    _check_quantity(B, 'B', 'electron_gyrofrequency', u.T)

    omega_ce = u.rad*(e*np.abs(B)/m_e).to(1/u.s)

    return omega_ce


def ion_gyrofrequency(B, ion='p'):
    r"""Calculate the ion gyrofrequency in units of radians per second.

    Parameters
    ----------
    B: Quantity
        The magnetic field magnitude in units convertible to tesla

    ion : string, optional
        Representation of the ion species.  If not given, then the ions
        are assumed to be protons.

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

    Examples
    --------
    >>> from astropy import units as u
    >>> ion_gyrofrequency(0.01*u.T)
    <Quantity 152451.87138666757 Hz>
    >>> ion_gyrofrequency(0.01*u.T, 'p')
    <Quantity 152451.87138666757 Hz>
    >>> ion_gyrofrequency(0.01*u.T, ion='T')

    """

    _check_quantity(B, 'B', 'ion_gyrofrequency', u.T)

    try:
        Z = charge_state(ion)
        if Z is None:
            Z = 1
    except Exception:
        raise ValueError("Unable to find charge state in "
                         "ion_gyrofrequency")

    try:
        m_i = ion_mass(ion)
    except Exception:
        raise ValueError("Unable to get ion mass in ion_gyrofrequency")

    omega_ci = u.rad * (Z*e*np.abs(B)/m_i).to(1/u.s)

    return omega_ci


def electron_gyroradius(B, Vperp_or_Te):
    r"""Returns the radius of gyration for an electron in a uniform
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

    Examples
    --------
    >>> from astropy import units as u
    >>> electron_gyroradius(0.2*u.T, 1e5*u.K)
    <Quantity 4.949493018143766e-05 m>
    >>> electron_gyroradius(5*u.uG, 1*u.eV)
    <Quantity 6744.259695124416 m>
    >>> electron_gyroradius(400*u.G, 1e7*u.m/u.s)
    <Quantity 0.00142140746360249 m>

    """

    if Vperp_or_Te.si.unit == 'T':
        B, Vperp_or_Te = Vperp_or_Te, B

    _check_quantity(B, 'B', 'electron_gyroradius', u.T)
    _check_quantity(Vperp_or_Te, 'Vperp_or_Te', 'electron_gyroradius',
                    [u.m/u.s, u.K])

    if Vperp_or_Te.unit in ['J', 'K']:
        T_e = Vperp_or_Te.to(u.K, equivalencies=u.temperature_energy())
        if T_e < 0*u.K:
            raise ValueError("An argument to electron_gyroradius corresponds "
                             "to a negative temperature.")
        Vperp = electron_thermal_speed(T_e)
    elif Vperp_or_Te.unit == 'm / s':
        Vperp = Vperp_or_Te

    omega_ce = electron_gyrofrequency(B)
    r_L = (Vperp/omega_ce).to(u.m, equivalencies=u.dimensionless_angles())

    return r_L


def ion_gyroradius(B, Vperp_or_Ti, ion='p'):
    r"""Returns the ion gyroradius.

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

    Examples
    --------
    >>> from astropy import units as u
    >>> ion_gyroradius(0.2*u.T, 1e5*u.K)
    <Quantity 0.0021208751836230026 m>
    >>> ion_gyroradius(0.2*u.T, 1e5*u.K, ion='p')
    <Quantity 0.0021208751836230026 m>
    >>> ion_gyroradius(5*u.uG, 1*u.eV, ion='alpha')
    <Quantity 288002.3964615791 m>
    >>> ion_gyroradius(400*u.G, 1e7*u.m/u.s, ion='Fe+++')
    <Quantity 48.23129633674924 m>

    """

    if Vperp_or_Ti.unit.si == 'T':
        B, Vperp_or_Ti = Vperp_or_Ti, B

    _check_quantity(B, 'B', 'ion_gyroradius', u.T)
    _check_quantity(Vperp_or_Ti, 'Vperp_or_Ti', 'ion_gyroradius',
                    [u.m/u.s, u.K])

    if Vperp_or_Ti.unit in ['J', 'K']:
        T_i = Vperp_or_Ti.to(u.K, equivalencies=u.temperature_energy())
        if T_i < 0*u.K:
            raise ValueError("An argument to ion_gyroradius corresponds "
                             "to a negative temperature.")
        Vperp = ion_thermal_speed(T_i, ion)
    elif Vperp_or_Ti.unit == 'm / s':
        Vperp = Vperp_or_Ti

    omega_ci = ion_gyrofrequency(B, ion)
    r_L = (Vperp/omega_ci).to(u.m, equivalencies=u.dimensionless_angles())

    return r_L


def electron_plasma_frequency(n_e):
    r"""Calculates the electron plasma frequency.

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
    >>> from astropy import units as u
    >>> from plasmapy import electron_plasma_frequency
    >>> electron_plasma_frequency(1e19*u.m**-3)
    <Quantity 178398636622.99567 rad / s>

    """

    _check_quantity(n_e, 'n_e', 'electron_plasma_frequency', u.m**-3,
                    can_be_negative=False)

    omega_pe = (u.rad*e*np.sqrt(n_e/(eps0*m_e))).to(u.rad/u.s)

    return omega_pe


def ion_plasma_frequency(n_i, ion='p'):
    r"""Calculates the ion plasma frequency.

    Parameters
    ----------
    n_i : Quantity
        Ion number density

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
    >>> from astropy import units as u
    >>> ion_plasma_frequency(1e19*u.m**-3)
    <Quantity 178398636622.99567 rad / s>
    >>> ion_plasma_frequency(1e19*u.m**-3, ion='p')
    <Quantity 178398636622.99567 rad / s>

    """

    _check_quantity(n_i, 'n_i', 'ion_plasma_frequency', u.m**-3,
                    can_be_negative=False)

    try:
        m_i = ion_mass(ion)
    except Exception:
        raise ValueError("Unable to get mass of ion in ion_gyrofrequency")

    try:
        Z = charge_state(ion)
    except Exception:
        raise ValueError("Unable to get charge state to calculate ion "
                         "plasma frequency.")

    omega_pi = u.rad*Z*e*np.sqrt(n_i/(eps0*m_e))

    return omega_pi.si


def Debye_length(T_e, n_e):
    r"""Calculate the Debye length.

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
    <Quantity 0.0021822555159125854 m>

    """

    _check_quantity(T_e, 'T_e', 'Debye_length', u.K, can_be_negative=False)

    lambda_D = ((eps0*k_B*T_e / (n_e * e**2))**0.5).to(u.m)

    return lambda_D


def Debye_number(T_e, n_e):
    r"""Returns the Debye number.

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
    >>> from astropy import units as u
    >>> Debye_number(5e6*u.K, 5e9*u.cm**-3)
    <Quantity 69282817.49483868>

    """

    _check_quantity(T_e, 'T_e', 'Debye_number', u.K, can_be_negative=False)
    _check_quantity(n_e, 'n_e', 'Debye_number', u.m**-3, can_be_negative=False)

    lambda_D = Debye_length(T_e, n_e)
    N_D = (4/3)*n_e*lambda_D**3

    return N_D.to(u.dimensionless_unscaled)


def ion_inertial_length(n_i, ion='p'):
    r"""Calculate the ion inertial length,

    Parameters
    ----------
    n_i : Quantity
        Ion number density in units convertible to

    ion : string, optional
        Representation of the ion species.  If not given, then the ions
        are assumed to be protons.

    Returns
    -------
    d_i : Quantity
        Ion inertial length in meters

    Notes
    -----
    The ion inertial length is also known as an ion skin depth and is
    given by:

    .. math::
    d_i = \frac{c}{\omega_{pi}}

    Example
    -------
    >>> from astropy import units as u
    >>> ion_inertial_length(5*u.m**-3, ion='He+')
    <Quantity 2376534.754 m>

    """

    try:
        Z = charge_state(ion)
    except Exception:
        raise ValueError("Unable to find charge state in ion_inertial_length.")

    _check_quantity(n_i, 'n_i', 'ion_inertial_length', u.m**-3,
                    can_be_negative=False)

    omega_pi = ion_plasma_frequency(n_i, ion=ion)
    d_i = (c/omega_pi).to(u.m, equivalencies=u.dimensionless_angles())

    return d_i


def electron_inertial_length(n_e):
    r"""Returns the electron inertial length.

    Parameters
    ----------
    n_e : Quantity
        Electron number density

    Returns
    -------
    d_e : Quantity
        Electron inertial length in meters

    Notes
    -----
    The electron inertial length is also known as an electron skin depth and
    is given by:

    .. math::
    d_e = \frac{c}{\omega_{pe}}

    Example
    -------
    >>> from astropy import units as u
    >>> electron_inertial_length(5*u.m**-3)
    <Quantity 2376534.754 m>

    """

    _check_quantity(n_e, 'n_e', 'electron_inertial_length', u.m**-3,
                    can_be_negative=False)

    omega_pe = electron_plasma_frequency(n_e)
    d_e = (c/omega_pe).to(u.m, equivalencies=u.dimensionless_angles())

    return d_e


def magnetic_pressure(B):
    r"""Calculate the magnetic pressure.

    Parameters
    ----------
    B: Quantity
        The magnetic field in units convertible to telsa

    Returns
    -------
    p_B: Quantity
        The magnetic pressure in units in pascals (newtons per square meter)

    Raises
    ------
    TypeError
        If the input is not a Quantity

    UnitsError
        If the input is not in units convertible to tesla

    ValueError
        If the magnetic field strength is not a real number between
        +/- infinity

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
    <Quantity 3978.873577297384 Pa>

    """

    _check_quantity(B, 'B', 'magnetic_pressure', u.T)

    p_B = (B**2/(2*mu0)).to(u.N/u.m**2)

    return p_B


def magnetic_energy_density(B: u.T):
    r"""Calculate the magnetic energy density.

    Parameters
    ----------
    B: Quantity
        The magnetic field in units convertible to tesla

    Returns
    -------
    E_B: Quantity
        The magnetic energy density in units of joules per cubic meter

    Raises
    ------
    TypeError
        If the input is not a Quantity

    UnitsError
        If the input is not in units convertible to tesla

    ValueError
        If the magnetic field strength is not a real number between
        +/- infinity

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
    magnetic_pressure : returns an equivalent Quantity, except in units
        of pascals.

    Example
    -------
    >>> from astropy import units as u
    >>> magnetic_energy_density(0.1*u.T)
    <Quantity 3978.873577297384 J / m3>

    """

    _check_quantity(B, 'B', 'magnetic_energy_density', u.T)

    E_B = (B**2/(2*mu0)).to(u.J/u.m**3)

    return E_B
