"""Functions to calculate plasma parameters."""

from astropy import units

from astropy.units import (UnitConversionError, UnitsError, quantity_input,
                           Quantity)

from ..constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi)
from ..atomic import (ion_mass, charge_state)

import numpy as np

# For future: change these into decorators.  _check_quantity does a
# bit more than @quantity_input as it can allow

from ..utils import _check_quantity, check_relativistic, check_quantity


r"""Values should be returned as an Astropy Quantity in SI units.

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


@check_relativistic
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
        Representation of the ion species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to protons.  If no charge state information is
        provided, then the ions are assumed to be singly charged.

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

    UserWarning
        If the Alfven velocity exceeds 10% of the speed of light, or
        if units are not provided and SI units are assumed.

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
    <Quantity 43173.87029559353 m / s>
    >>> Alfven_speed(B, rho, ion)
    <Quantity 43173.87029559353 m / s>
    >>> Alfven_speed(B, rho, ion).to(u.cm/u.us)
    <Quantity 4.317387029559353 cm / us>
    """

    _check_quantity(B, 'B', 'Alfven_speed', units.T)
    _check_quantity(density, 'density', 'Alfven_speed',
                    [units.m**-3, units.kg/units.m**3], can_be_negative=False)

    B = B.to(units.T)
    density = density.si

    if density.unit == units.m**-3:
        try:
            m_i = ion_mass(ion)
            Z = charge_state(ion)
            if Z is None:
                Z = 1
        except Exception:
            raise ValueError("Invalid ion in Alfven_speed.")
        rho = density*m_i + Z*density*m_e

    elif density.unit == units.kg/units.m**3:
        rho = density

    try:
        V_A = (np.abs(B)/np.sqrt(mu0*rho)).to(units.m/units.s)
    except Exception:
        raise ValueError("Unable to find Alfven speed")

    return V_A


@check_relativistic
@check_quantity({
    'T_i': {'units': units.K, 'can_be_negative': False},
    'T_e': {'units': units.K, 'can_be_negative': False}
})
def ion_sound_speed(*ignore, T_e=0*units.K, T_i=0*units.K,
                    gamma_e=1, gamma_i=3, ion='p'):
    r"""Returns the ion sound speed for an electron-ion plasma.

    Parameters
    ----------
    T_e : Quantity, optional
        Electron temperature in units of temperature or energy per
        particle.  If this is not given, then the electron temperature
        is assumed to be zero.  If only one temperature is entered, it
        is assumed to be the electron temperature.

    T_i : Quantity, optional
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

    Returns
    -------
    V_S : Quantity
        The ion sound speed in units of meters per second.

    Raises
    ------
    TypeError
        If any of the arguments are not entered as keyword arguments
        or are of an incorrect type.

    ValueError
        If the ion mass, adiabatic index, or temperature are invalid.

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
    <Quantity 203155.07640420322 m / s>
    >>> ion_sound_speed(T_e=5e6*u.K)
    <Quantity 203155.07640420322 m / s>
    >>> ion_sound_speed(T_e=500*u.eV, T_i=200*u.eV, ion='D+')
    <Quantity 229586.01860212447 m / s>

    """

    if ignore:
        raise TypeError("All arguments are required to be keyword arguments "
                        "in ion_sound_speed to prevent mixing up the electron "
                        "and ion temperatures. An example call that uses the "
                        "units subpackage from astropy is: "
                        "ion_sound_speed(T_e=5*units.K, T_i=0*units.K, "
                        "ion='D+')")

    try:
        m_i = ion_mass(ion)
        Z = charge_state(ion)
        if Z is None:
            Z = 1
    except Exception:
        raise ValueError("Invalid ion in ion_sound_speed.")

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

    T_i = T_i.to(units.K, equivalencies=units.temperature_energy())
    T_e = T_e.to(units.K, equivalencies=units.temperature_energy())

    try:
        V_S_squared = (gamma_e*Z*k_B*T_e + gamma_i*k_B*T_i)/m_i
        V_S = np.sqrt(V_S_squared).to(units.m/units.s)
    except Exception:
        raise ValueError("Unable to find ion sound speed.")

    return V_S


@check_relativistic
@check_quantity({
    'T_e': {'units': units.K, 'can_be_negative': False}
})
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
        If the electron temperature is not in units of temperature or
        energy per particle

    ValueError
        The electron temperature is invalid

    UserWarning
        If the electron thermal speed exceeds 10% of the speed of
        light, or if units are not provided and SI units are assumed.

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
    <Quantity 1326205.1212395933 m / s>
    >>> electron_thermal_speed(1e6*u.K)
    <Quantity 5505693.988425379 m / s>

    """

    T_e = T_e.to(units.K, equivalencies=units.temperature_energy())
    V_Te = (np.sqrt(2*k_B*T_e/m_e)).to(units.m/units.s)

    return V_Te


@check_relativistic
@check_quantity({
    'T_i': {'units': units.K, 'can_be_negative': False}
})
def ion_thermal_speed(T_i, ion='p'):
    r"""Returns the most probable speed for an ion within a Maxwellian
    distribution.

    Parameters
    ----------
    T_i : Quantity
        The ion temperature in either kelvin or energy per particle

    ion : string, optional
        Representation of the ion species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to protons.  If no charge state information is
        provided, then the ions are assumed to be singly charged.

    Returns
    -------
    V_Ti : Quantity
        Ion thermal speed

    Raises
    ------
    TypeError
        The ion temperature is not a Quantity

    UnitConversionError
        If the ion temperature is not in units of temperature or
        energy per particle

    ValueError
        The ion temperature is invalid or ion cannot be used to
        identify an isotope or ion

    UserWarning
        If the ion thermal speed exceeds 10% of the speed of light, or
        if units are not provided and SI units are assumed.

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
    <Quantity 30949.690182856546 m / s>
    >>> ion_thermal_speed(1e6*u.K, ion='p')
    <Quantity 128486.55193256242 m / s>

    """

    T_i = T_i.to(units.K, equivalencies=units.temperature_energy())

    try:
        m_i = ion_mass(ion)
    except Exception:
        raise ValueError("Unable to find ion mass in ion_thermal_speed")

    V_Ti = (np.sqrt(2*k_B*T_i/m_i)).to(units.m/units.s)

    return V_Ti


@check_quantity({
    'B': {'units': units.T}
})
def electron_gyrofrequency(B):
    r"""Calculate the electron gyrofrequency in units of radians per second.

    Parameters
    ----------
    B: Quantity
        The magnetic field magnitude in units convertible to tesla.

    Returns
    -------
    omega_ce: Quantity
        Electron gyrofrequency in radians per second.

    Raises
    ------
    TypeError
        The magnetic field is not a Quantity

    UnitConversionError
        If the magnetic field is in incorrect units

    ValueError
        If the magnetic field has an invalid value

    UserWarning
        If units are not provided and SI units are assumed

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
    17588200236.02124 rad / s
    >>> f_ce = omega_ce.to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])
    >>> print(f_ce)
    2799249007.6528206 Hz

    """

    omega_ce = units.rad*(e*np.abs(B)/m_e).to(1/units.s)

    return omega_ce


@check_quantity({
    'B': {'units': units.T}
})
def ion_gyrofrequency(B, ion='p'):
    r"""Calculate the ion gyrofrequency in units of radians per second.

    Parameters
    ----------
    B: Quantity
        The magnetic field magnitude in units convertible to tesla.

    ion : string, optional
        Representation of the ion species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to protons.  If no charge state information is
        provided, then the ions are assumed to be singly charged.

    Returns
    -------
    omega_ci: Quantity
        The ion gyrofrequency in units of radians per second

    Raises
    ------
    TypeError
        If the magnetic field is not a Quantity or ion is not of an
        appropriate type

    ValueError
        If the magnetic field contains invalid values or ion cannot be
        used to identify an ion or isotope

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----
    The ion gyrofrequency is the angular frequency of ion gyration
    around magnetic field lines and is given by:

    .. math::
    omega_{ci} = \frac{Z e B}{m_i}

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
    <Quantity 957883.3224148067 rad / s>
    >>> ion_gyrofrequency(0.01*u.T, 'p')
    <Quantity 957883.3224148067 rad / s>
    >>> ion_gyrofrequency(0.01*u.T, ion='T')
    <Quantity 319964.54975910933 rad / s>
    """

    try:
        m_i = ion_mass(ion)
        Z = charge_state(ion)
        if Z is None:
            Z = 1
    except Exception:
        raise ValueError("Invalid ion in ion_gyrofrequency")

    omega_ci = units.rad * (Z*e*np.abs(B)/m_i).to(1/units.s)

    return omega_ci


def electron_gyroradius(B, *args, Vperp=None, T_e=None):
    r"""Returns the radius of gyration for an electron in a uniform
    magnetic field.

    Parameters
    ----------
    B : Quantity
        The magnetic field magnitude in units convertible to tesla.
        If no units are given, a UserWarning will be raised and units
        of tesla will be assumed.

    Vperp : Quantity, optional
        The component of electron velocity that is perpendicular to
        the magnetic field in units convertible to meters per second.

    T_e : Quantity, optional
        The electron temperature in units convertible to kelvin.

    args : Quantity
        If the second positional argument is a Quantity with units
        appropriate to Vperp or T_e, then this argument will take the
        place of that keyword argument.

    Returns
    -------
    r_Le : Quantity
        The electron gyroradius in units of meters.  This Quantity
        will be based on either the perpendicular component of
        electron velocity as inputted, or the most probable speed for
        an electron within a Maxwellian distribution for the electron
        temperature.

    Raises
    ------
    TypeError
        If either of the inputs is not a Quantity

    UnitConversionError
        If either argument is in incorrect units

    ValueError
        If either argument contains invalid values

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----
    One but not both of Vperp and T_i must be inputted.

    If any of B, Vperp, or T_e is a number rather than a Quantity,
    then SI units will be assumed and a UserWarning will be raised.

    Formula
    -----
    The electron gyroradius is also known as the Larmor radius for
    electrons and is given by:

    .. math::
    r_{Le} = \frac{V_{perp}}{omega_{ce}}

    where :math:`V_{\perp}` is the component of electron velocity that
    is perpendicular to the magnetic field and :math:`\omega_{ce}` is
    the electron gyrofrequency.  If a temperature is provided, then
    :math:`V_\perp` will be the most probable thermal velocity of an
    electron at that temperature.

    Examples
    --------
    >>> from astropy import units as u
    >>> electron_gyroradius(B = 0.01*u.T, T_e = 1e6*u.K)
    <Quantity 0.0031303339253265536 m>
    >>> electron_gyroradius(B = 0.01*u.T, Vperp = 1e6*u.m/u.s)
    <Quantity 0.0005685630062091092 m>
    >>> electron_gyroradius(0.2*u.T, 1e5*u.K)
    <Quantity 4.9494925204636764e-05 m>
    >>> electron_gyroradius(5*u.uG, 1*u.eV)
    <Quantity 6744.259818299466 m>
    >>> electron_gyroradius(400*u.G, 1e7*u.m/u.s)
    <Quantity 0.0014214075155227729 m>

    """

    if Vperp is not None and T_e is not None:
        raise ValueError("Cannot have both Vperp and T_e as arguments to "
                         "electron_gyroradius")

    if len(args) == 1 and isinstance(args[0], units.Quantity):
        arg = args[0].si
        if arg.unit == units.T and B.si.unit in [units.J, units.K,
                                                 units.m/units.s]:
            B, arg = arg, B

        if arg.unit == units.m/units.s:
            Vperp = arg
        elif arg.unit in (units.J, units.K):
            T_e = arg.to(units.K, equivalencies=units.temperature_energy())
        else:
            raise units.UnitConversionError("Incorrect units for positional "
                                            "argument in electron_gyroradius")
    elif len(args) > 0:
        raise ValueError("Incorrect inputs to electron_gyroradius")

    _check_quantity(B, 'B', 'electron_gyroradius', units.T)

    if Vperp is not None:
        _check_quantity(Vperp, 'Vperp', 'electron_gyroradius', units.m/units.s)
    elif T_e is not None:
        _check_quantity(T_e, 'T_e', 'electron_gyroradius', units.K)
        Vperp = electron_thermal_speed(T_e)

    omega_ce = electron_gyrofrequency(B)
    r_L = np.abs(Vperp)/omega_ce

    return r_L.to(units.m, equivalencies=units.dimensionless_angles())


def ion_gyroradius(B, *args, Vperp=None, T_i=None, ion='p'):
    r"""Returns the ion gyroradius.

    Parameters
    ----------
    B: Quantity
        The magnetic field magnitude in units convertible to tesla.

    Vperp: Quantity, optional
        The component of ion velocity that is perpendicular to the
        magnetic field in units convertible to meters per second.

    T_i: Quantity, optional
        The ion temperature in units convertible to kelvin.

    ion : string, optional
        Representation of the ion species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to protons.  If no charge state information is
        provided, then the ions are assumed to be singly charged.

    args : Quantity
        If the second positional argument is a Quantity with units
        appropriate to Vperp or T_i, then this argument will take the
        place of that keyword argument.

    Returns
    -------
    r_Li : Quantity
        The ion gyroradius in units of meters.  This Quantity will be
        based on either the perpendicular component of ion velocity as
        inputted, or the most probable speed for an ion within a
        Maxwellian distribution for the ion temperature.

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

    Formula
    -------
    The ion gyroradius is also known as the ion Larmor radius and is
    given by

    .. math::
    r_{Li} = \frac{V_{\perp}}{omega_{ci}}

    where :math:`V_{\perp}` is the component of ion velocity that is
    perpendicular to the magnetic field and :math:`\omega_{ci}` is the
    ion gyrofrequency.  If a temperature is provided, then
    :math:`V_\perp` will be the most probable thermal velocity of an
    ion at that temperature.

    Examples
    --------
    >>> from astropy import units as u
    >>> ion_gyroradius(0.2*u.T, 1e5*u.K)
    <Quantity 0.002120874971411475 m>
    >>> ion_gyroradius(0.2*u.T, 1e5*u.K, ion='p')
    <Quantity 0.002120874971411475 m>
    >>> ion_gyroradius(5*u.uG, 1*u.eV, ion='alpha')
    <Quantity 288002.38837768475 m>
    >>> ion_gyroradius(400*u.G, 1e7*u.m/u.s, ion='Fe+++')
    <Quantity 48.23129811339086 m>

    """

    if Vperp is not None and T_i is not None:
        raise ValueError("Cannot have both Vperp and T_i as arguments to "
                         "ion_gyroradius")

    if len(args) == 1 and isinstance(args[0], units.Quantity):
        arg = args[0].si
        if arg.unit == units.T and B.si.unit in [units.J, units.K,
                                                 units.m/units.s]:
            B, arg = arg, B

        if arg.unit == units.m/units.s:
            Vperp = arg
        elif arg.unit in (units.J, units.K):
            T_i = arg.to(units.K, equivalencies=units.temperature_energy())
        else:
            raise units.UnitConversionError("Incorrect units for positional "
                                            "argument in ion_gyroradius")
    elif len(args) > 0:
        raise ValueError("Incorrect inputs to ion_gyroradius")

    _check_quantity(B, 'B', 'ion_gyroradius', units.T)

    if Vperp is not None:
        _check_quantity(Vperp, 'Vperp', 'ion_gyroradius', units.m/units.s)
    elif T_i is not None:
        _check_quantity(T_i, 'T_i', 'ion_gyroradius', units.K)
        Vperp = ion_thermal_speed(T_i, ion=ion)

    omega_ci = ion_gyrofrequency(B, ion)

    r_Li = np.abs(Vperp)/omega_ci

    return r_Li.to(units.m, equivalencies=units.dimensionless_angles())


@check_quantity({
    'n_e': {'units': units.m**-3, 'can_be_negative': False}
})
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

    Raises
    ------
    TypeError
        If n_e is not a Quantity

    UnitConversionError
        If n_e is in incorrect units

    ValueError
        If n_e contains invalid values

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----
    In a simple one dimensional model, the separation of charge
    within a plasma creates an electric field proportional to
    separation distance with magnitude:

    .. math::
    E = \frac{e}{\epsilon_0} n_e x

    Where x is the separation distance.

    The electrons will move under the action of this field with force
    magnitude:

    .. math::
    F = e E = \frac{e^2}{\epsilon_0} n x = m_e \frac{d^2 x}{d t^2}

    This is a simple harmonic oscillator. Computing its eigenfrequency
    yields the electron plasma frequency, which is given by:

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
    >>> electron_plasma_frequency(1e19*u.m**-3)
    <Quantity 178398636471.3789 rad / s>

    """

    omega_pe = (units.rad*e*np.sqrt(n_e/(eps0*m_e))).to(units.rad/units.s)

    return omega_pe


@check_quantity({
    'n_i': {'units': units.m**-3, 'can_be_negative': False}
})
def ion_plasma_frequency(n_i, ion='p'):
    r"""Calculates the ion plasma frequency.

    Parameters
    ----------
    n_i : Quantity
        Ion number density in units convertible to per cubic meter

    ion : string, optional
        Representation of the ion species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to protons.  If no charge state information is
        provided, then the ions are assumed to be singly charged.

    Returns
    -------
    omega_pi : Quantity
        The ion plasma frequency in radians per second.

    Raises
    ------
    TypeError
        If n_i is not a Quantity or ion is not of an appropriate type

    UnitConversionError
        If n_i is not in correct units

    ValueError
        If n_i contains invalid values or ion cannot be used to
        identify an ion or isotope.

    UserWarning
        If units are not provided and SI units are assumed

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
    <Quantity 4163294530.6925354 rad / s>
    >>> ion_plasma_frequency(1e19*u.m**-3, ion='p')
    <Quantity 4163294530.6925354 rad / s>

    """

    try:
        m_i = ion_mass(ion)
        Z = charge_state(ion)
        if Z is None:
            Z = 1
    except Exception:
        raise ValueError("Invalid ion in ion_gyrofrequency")

    omega_pi = units.rad*Z*e*np.sqrt(n_i/(eps0*m_i))

    return omega_pi.si


@check_quantity({
    'T_e': {'units': units.K, 'can_be_negative': False},
    'n_e': {'units': units.m**-3, 'can_be_negative': False}
})
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
        The Debye length in meters

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
    <Quantity 0.002182255218625608 m>

    """

    T_e = T_e.to(units.K, equivalencies=units.temperature_energy())

    try:
        lambda_D = ((eps0*k_B*T_e/(n_e*e**2))**0.5).to(units.m)
    except Exception:
        raise ValueError("Unable to find Debye length.")

    return lambda_D


@check_quantity({
    'T_e': {'units': units.K, 'can_be_negative': False},
    'n_e': {'units': units.m**-3, 'can_be_negative': False}
})
def Debye_number(T_e, n_e):
    r"""Returns the Debye number.

    Parameters
    ----------
    T_e: Quantity
        Electron temperature

    n_e: Quantity
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
    <Quantity 217658301.50749832>

    """

    try:
        lambda_D = Debye_length(T_e, n_e)
        N_D = (4/3)*np.pi*n_e*lambda_D**3
    except Exception:
        raise ValueError("Unable to find Debye number")

    return N_D.to(units.dimensionless_unscaled)


@check_quantity({
    'n_i': {'units': units.m**-3, 'can_be_negative': False}
})
def ion_inertial_length(n_i, ion='p'):
    r"""Calculate the ion inertial length,

    Parameters
    ----------
    n_i : Quantity
        Ion number density in units convertible to m**-3

    ion : string, optional
        Representation of the ion species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to protons.  If no charge state information is
        provided, then the ions are assumed to be singly charged.

    Returns
    -------
    d_i : Quantity
        Ion inertial length in meters

    Raises
    ------
    TypeError
        If n_i not a Quantity or ion is not a string

    UnitConversionError
        If n_i is not in units of a number density

    ValueError
        The ion density does not have an appropriate value.

    UserWarning
        If units are not provided and SI units are assumed

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
    <Quantity 202985801.8507889 m>

    """

    try:
        Z = charge_state(ion)
    except Exception:
        raise ValueError("Invalid ion in ion_inertial_length.")

    omega_pi = ion_plasma_frequency(n_i, ion=ion)
    d_i = (c/omega_pi).to(units.m, equivalencies=units.dimensionless_angles())

    return d_i


@check_quantity({
    'n_e': {'units': units.m**-3, 'can_be_negative': False}
})
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

    Raises
    ------
    TypeError
        If n_e is not a Quantity

    UnitConversionError
        If n_e is not in units of per cubic meter

    ValueError
        If n_e contains invalid values

    UserWarning
        If units are not provided and SI units are assumed

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
    <Quantity 2376534.756019761 m>

    """

    omega_pe = electron_plasma_frequency(n_e)
    d_e = (c/omega_pe).to(units.m, equivalencies=units.dimensionless_angles())

    return d_e


@check_quantity({
    'B': {'units': units.T}
})
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

    UnitConversionError
        If the input is not in units convertible to tesla

    ValueError
        If the magnetic field strength is not a real number between
        +/- infinity

    UserWarning
        If units are not provided and SI units are assumed

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

    p_B = (B**2/(2*mu0)).to(units.Pa)

    return p_B


@check_quantity({
    'B': {'units': units.T}
})
def magnetic_energy_density(B: units.T):
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

    UnitConversionError
        If the input is not in units convertible to tesla

    ValueError
        If the magnetic field strength does not have an appropriate
        value.

    UserWarning
        If units are not provided and SI units are assumed

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

    E_B = (B**2/(2*mu0)).to(units.J/units.m**3)

    return E_B


@check_quantity({
    'B': {'units': units.T},
    'n_e': {'units': units.m**-3, 'can_be_negative': False}
})
def upper_hybrid_frequency(B, n_e):
    r"""Returns the upper hybrid frequency.

    Parameters
    ----------
    B : Quantity
        The magnetic field magnitude in units convertible to tesla.

    n_e : Quantity
        The electron number density

    Returns
    -------
    omega_uh : Quantity
        The upper hybrid frequency in radians per second

    Raises
    ------
    TypeError
        If either of B or n_e is not a Quantity

    UnitConversionError
        If either of B or n_e is in incorrect units

    ValueError
        If either of B or n_e contains invalid values or are of
        incompatible dimensions

    UserWarning
        If units are not provided and SI units are assumed

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
    <Quantity 400459419447.72076 rad / s>

    """

    try:
        omega_pe = electron_plasma_frequency(n_e=n_e)
        omega_ce = electron_gyrofrequency(B)
        omega_uh = (np.sqrt(omega_pe**2 + omega_ce**2)).to(units.rad/units.s)
    except Exception:
        raise ValueError("Unable to find upper hybrid frequency.")

    return omega_uh


@check_quantity({
    'B': {'units': units.T},
    'n_i': {'units': units.m**-3, 'can_be_negative': False}
})
def lower_hybrid_frequency(B, n_i, ion='p'):
    r"""Returns the lower hybrid frequency.

    Parameters
    ----------
    B : Quantity
        The magnetic field magnitude in units convertible to tesla.

    n_i : Quantity
        Ion number density

    ion : string, optional
        Representation of the ion species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to protons.  If no charge state information is
        provided, then the ions are assumed to be singly charged.

    Returns
    -------
    omega_lh : Quantity
        The lower hybrid frequency in radians per second

    Raises
    ------
    TypeError
        If either of B or n_i is not a Quantity, or ion is of an
        inappropriate type

    UnitConversionError
        If either of B or n_i is in incorrect units

    ValueError
        If either of B or n_i contains invalid values or are of
        incompatible dimensions, or ion cannot be used to identify an
        ion or isotope

    UserWarning
        If units are not provided and SI units are assumed

    Notes
    -----
    The lower hybrid frequency is given through the relation

    .. math::
    \frac{1}{\omega_{lh}^2} = \frac{1}{\omega_{ci}^2 + \omega_{pi}^2}
    + \frac{1}{\omega_{ci}\omega_{ce}}

    where .. math::`\omega_{ci}` is the ion gyrofrequency,
    .. math::`\omega_{ce}` is the electron gyrofrequency, and
    .. math::`\omega_{pi}` is the ion plasma frequency.

    Example
    -------
    >>> from astropy import units as u
    >>> lower_hybrid_frequency(0.2*u.T, n_i=5e19*u.m**-3, ion='D+')
    <Quantity 578372732.5155494 rad / s>

    """

    # We do not need a charge state here, so the sole intent is to
    # catch invalid ions.
    try:
        charge_state(ion)
    except Exception:
        raise ValueError("Invalid ion in lower_hybrid_frequency.")

    try:
        omega_ci = ion_gyrofrequency(B, ion=ion)
        omega_pi = ion_plasma_frequency(n_i, ion=ion)
        omega_ce = electron_gyrofrequency(B)
        omega_lh = 1/np.sqrt((omega_ci*omega_ce)**-1+omega_pi**-2)
        omega_lh = omega_lh.to(units.rad/units.s)
    except Exception:
        raise ValueError("Unable to find lower hybrid frequency.")

    return omega_lh
