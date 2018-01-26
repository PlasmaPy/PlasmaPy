"""Functions that retrieve or are related to elemental or isotopic data."""

import numpy as np
import warnings
from typing import (Union, Optional, List)

from astropy import units as u, constants as const
from astropy.units import Quantity

from .elements import _Elements
from .isotopes import _Isotopes

from ..utils import (
    AtomicWarning,
    InvalidElementError,
    InvalidIsotopeError,
    InvalidIonError,
    AtomicError,
    MissingAtomicDataError,
    MissingAtomicDataWarning,
    InvalidParticleError,
    ChargeError)

from .names import (
    atomic_symbol,
    isotope_symbol,
    _extract_integer_charge,
    _is_proton,
    _is_positron,
    _is_antineutron,
    _is_antiproton,
    _is_electron,
    _is_neutron)


# TODO: refactor the atomic subpackage
#  - Create a Particle class
#  - Create a decorator that:
#     1. Takes the inputs of a function corresponding to a particle
#     2. Checks to make sure the particle is valid
#     3. Checks to make sure that the inputs are of the correct types
#     4. Returns the Particle class which would be the input to the
#        inner function

# TODO: Create a particle_mass function
# TODO: Create lepton_number and baryon_number functions


def atomic_number(argument: str) -> str:
    r"""Returns the number of protons in an atom, isotope, or ion.

    Parameters
    ----------

    argument: string
        A string representing an element, isotope, or ion.

    Returns
    -------

    atomic_number: integer
        An integer representing the atomic number of the element or
        isotope.

    Raises
    ------

    InvalidElementError
        If the argument is a valid particle but not a valid element.

    InvalidParticleError
        If the argument does not correspond to a valid particle.

    TypeError
        If the argument is not a string.

    See also
    --------

    mass_number : returns the mass number (the total number of protons
        and neutrons) of an isotope.

    Examples
    --------
    >>> atomic_number("H")
    1
    >>> atomic_number("tritium")
    1
    >>> atomic_number("alpha")
    2
    >>> atomic_number("oganesson")
    118

    """

    try:
        element = atomic_symbol(argument)
        atomic_numb = _Elements[element]['atomic_number']
    except (InvalidElementError, KeyError):
        raise InvalidElementError(
            f"The argument {argument} to atomic_number does not correspond to "
            "a valid element.") from None
    except InvalidParticleError:
        raise InvalidParticleError(f"The argument {argument} to atomic_number "
                                   "is not a valid particle.") from None
    except TypeError:
        raise TypeError(f"The argument {argument} to atomic number is not a "
                        f"string.") from None

    return atomic_numb


def mass_number(isotope: str) -> int:
    r"""Get the mass number (the number of protons and neutrons) of an
    isotope.

    Parameters
    ----------

    isotope : string
        A string representing an isotope or a neutron.

    Returns
    -------

    mass_number : integer
       The total number of protons plus neutrons in a nuclide.

    Raises
    ------

    InvalidParticleError
        If the argument does not correspond to a valid particle.

    InvalidIsotopeError
        If the argument does not correspond to a valid isotope.

    TypeError
        The argument is not a string.


    See also
    --------

    atomic_number : returns the number of protons in an isotope or
        element

    Examples
    --------

    >>> mass_number("H-1")
    1
    >>> mass_number("Pb-208")
    208
    >>> mass_number("tritium")
    3
    >>> mass_number("n")
    1
    >>> mass_number("alpha")
    4

    """

    try:
        isotope = isotope_symbol(isotope)
        mass_numb = _Isotopes[isotope]["mass_number"]
    except InvalidIsotopeError:
        raise InvalidIsotopeError("Invalid isotope in mass_number.")
    except InvalidParticleError:
        raise InvalidParticleError(
            f"The argument {isotope} to mass_number does not "
            "correspond to a valid isotope.") from None
    except TypeError:
        raise TypeError("The argument to mass_number is not a string.")

    return mass_numb


def standard_atomic_weight(argument: Union[str, int]) -> Quantity:
    r"""Returns the standard (conventional) atomic weight of an element
    based on the relative abundances of isotopes in terrestrial
    environments.

    Parameters
    ----------

    argument: string or integer
        A string representing an element or an integer representing an
        atomic number

    Returns
    -------

    atomic_weight: astropy.units.Quantity with units of u
        The standard atomic weight of an element based on values from
        NIST

    Raises
    ------

    InvalidElementError
        If the argument is a valid particle but not a valid element.

    InvalidParticleError
        If the argument does not correspond to a valid particle.

    TypeError
        If the argument is not a string or integer.

    See also
    --------

    isotope_mass : returns the atomic mass of an isotope.

    ion_mass : returns the mass of an ion of an element or isotope,
        accounting for reduction of mass from the neutral state due to
        different numbers of electrons.

    Notes
    -----

    Standard atomic weight data are most readily available for the
    terrestrial environment, so this function may not be wholly
    appropriate for space and astrophysical environments.

    The relative abundances of different isotopes of an element
    sometimes vary naturally in different locations within the
    terrestrial environment.  The CIAAW provides ranges for these
    element, which include H, Li, B, C, N, O, Mg, Si, S, Cl, Br, Tl.
    This function provides a single value from the CIAWW 2015 standard
    values when a single value is given, and the lower accuracy
    conventional value given by Meija et al. (2013,
    doi:10.1515/pac-2015-0305) for the elements where a range is
    given.

    Examples
    --------

    >>> from astropy import units as u
    >>> standard_atomic_weight("H")
    <Quantity 1.008 u>
    >>> # the following result accounts for small amount of deuterium
    >>> standard_atomic_weight("H").to(u.kg)
    <Quantity 1.67382335e-27 kg>
    >>> isotope_mass("H-1")
    <Quantity 1.00782503 u>
    >>> standard_atomic_weight(82)
    <Quantity 207.2 u>
    >>> standard_atomic_weight("lead")
    <Quantity 207.2 u>

    """

    try:
        element = atomic_symbol(argument)
    except InvalidParticleError:
        raise InvalidParticleError(f"{argument} is an invalid argument to "
                                   "standard_atomic_weight.")
    except InvalidElementError:
        raise InvalidElementError(f"{argument} is not a valid element in "
                                  "standard_atomic_weight.")

    try:
        integer_charge(argument)
    except ChargeError:
        pass
    else:
        raise AtomicError("Use ion_mass to get masses of ions.")

    try:
        isotope_symbol(argument)
    except InvalidIsotopeError:
        pass
    else:
        raise AtomicError("Use isotope_mass to get masses of isotopes.")

    try:
        atomic_weight = _Elements[element]['atomic_mass']
    except KeyError as e:
        raise MissingAtomicDataError(
            f"No standard atomic weight is available for {element}.") from None

    return atomic_weight


def isotope_mass(argument: Union[str, int],
                 mass_numb: int = None) -> Quantity:
    r"""Return the mass of an isotope.

    Parameters
    ----------

    argument : string or integer
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    mass_numb : integer (optional)
        The mass number of the isotope.

    Returns
    -------

    isotope_mass : Quantity
        The atomic mass of a neutral atom of an isotope.

    Raises
    ------

    InvalidIsotopeError
        If the argument is a valid particle but not a valid isotope.

    InvalidParticleError
        If the argument does not correspond to a valid particle.

    AtomicError
        If the charge of the particle is given, in which case ion_mass
        should be used instead.

    TypeError
        If the argument is not a string.

    See also
    --------

    standard_atomic_weight : returns atomic weight of an element based
        on terrestrial abundances of isotopes

    ion_mass : returns the mass of an ion of an element or isotope,
        accounting for loss of electrons

    Notes
    -----

    The masses of rare isotopes may be unavailable.

    Examples
    --------

    >>> from astropy import units as u
    >>> isotope_mass("H-1")
    <Quantity 1.00782503 u>
    >>> isotope_mass("H-1").to(u.kg)
    <Quantity 1.67353281e-27 kg>
    >>> isotope_mass("He", 4)
    <Quantity 4.00260325 u>
    >>> isotope_mass(2, 4)
    <Quantity 4.00260325 u>

    """

    argument, Z = _extract_integer_charge(argument)

    if Z is not None and Z != 0:
        raise AtomicError("Use ion_mass instead of isotope_mass for masses of "
                          "charged particles.")

    try:
        isotope = isotope_symbol(argument, mass_numb)
        atomic_mass = _Isotopes[isotope]['atomic_mass']
    except InvalidParticleError:
        raise InvalidParticleError("Invalid particle in isotope_mass.")
    except InvalidIsotopeError:
        raise InvalidIsotopeError("Cannot identify isotope in isotope_mass.")
    except TypeError:
        raise TypeError("Invalid input to isotope_mass")

    return atomic_mass


def ion_mass(argument: Union[str, int, Quantity], Z: int = None,
             mass_numb: int = None) -> Quantity:
    r"""Returns the mass of an ion by finding the standard atomic
    weight of an element or the atomic mass of an isotope, and then
    accounting for the change in mass due to loss of electrons from
    ionization.

    Parameters
    ----------

    argument: string, integer, or Quantity
        A string representing an element, isotope, or ion; an integer
        representing an atomic number; or a Quantity in units of mass
        within the mass range that is appropriate for an ion.

    Z: integer (optional)
        The ionization state of the ion (defaulting to a charge of
        Z=1)

    mass_numb: integer (optional)
        The mass number of an isotope.

    Returns
    -------

    m_i: Quantity
        The mass of a single ion of the isotope or element with charge
        state Z.

    Raises
    ------

    TypeError
        The argument is not a string, integer, or Quantity.

    InvalidIonError
        If the argument represents a particle other than an ion, the
        ionization state exceeds the atomic number, or no isotope mass
        or standard atomic weight is available.

    AtomicWarning
        If a mass was inputted and it is outside of the range of known
        isotopes or electrons/positrons.

    UnitConversionError
        If the argument is a Quantity but does not have units of mass.

    See also
    --------

    standard_atomic_mass : returns the conventional atomic mass of an
        element based on terrestrial values and assuming the atom is
        neutral.

    isotope_mass : returns the mass of an isotope (if available)
        assuming the atom is neutral.

    Notes
    -----

    This function in general finds the mass of an isotope (or the
    standard atomic weight based on terrestrial values if a unique
    isotope cannot be identified), and then substracts the mass of Z
    electrons.  If Z is not provided as an input, then this function
    assumes that the ion is singly ionized.

    Specific values are returns for protons, deuterons, tritons, alpha
    particles, and positrons.

    Calling ion_mass('H') does not return the mass of a proton but
    instead uses hydrogen's standard atomic weight based on
    terrestrial values.  To get the mass of a proton, use
    ion_mass('p').

    This function can accept a Quantity in units of mass.  If the
    Quantity is close to the mass of an electron or positron, it will
    return the mass of an electron or positron to full known
    precision.  If the Quantity is within the mass range of known
    isotopes, it will return the mass that is inputted to it.

    Examples
    --------

    >>> print(ion_mass('p').si.value)
    1.672621898e-27
    >>> ion_mass('H+')  # assumes terrestrial abundance of D
    <Quantity 1.67291241e-27 kg>
    >>> ion_mass('H+') == ion_mass('p')
    False
    >>> ion_mass('H-1') == ion_mass('p')
    True
    >>> ion_mass('P+')  # phosphorus
    <Quantity 5.14322301e-26 kg>
    >>> ion_mass('He-4', Z=2)
    <Quantity 6.64465709e-27 kg>
    >>> ion_mass('T+')
    <Quantity 5.00735666e-27 kg>
    >>> ion_mass(26, Z=1, mass_numb=56)
    <Quantity 9.28812345e-26 kg>
    >>> ion_mass('Fe-56 1+')
    <Quantity 9.28812345e-26 kg>
    >>> ion_mass(9.11e-31*u.kg).si.value
    9.10938356e-31
    >>> ion_mass(1.67e-27*u.kg)
    <Quantity 1.67e-27 kg>

    """

    if isinstance(argument, Quantity) and Z is None and mass_numb is None:

        try:
            m_i = argument.to(u.kg)
        except u.UnitConversionError:
            raise u.UnitConversionError("If the ion in given as a Quantity, "
                                        "then it must have units of mass.")

        if np.isclose(m_i.value, const.m_e.value, atol=1e-33):  # positrons
            return const.m_e
        elif 1.66e-27 <= m_i.value < 7e-25:  # mass range of known isotopes
            return m_i
        else:
            warnings.warn(
                "The mass that was inputted to ion_mass and is being returned"
                " from ion_mass is outside of the range of known isotopes or "
                "electrons/ions.", AtomicWarning)
            return m_i

    if _is_electron(argument) or _is_positron(argument):
        return const.m_e
    elif _is_proton(argument, Z, mass_numb) or _is_antiproton(argument):
        return const.m_p
    elif _is_neutron(argument, mass_numb):
        raise InvalidIonError("Use isotope_mass or m_n to get mass of neutron")

    if isinstance(argument, str):
        argument, Z_from_arg = _extract_integer_charge(argument)
    else:
        Z_from_arg = None

    if Z is None and Z_from_arg is None:
        Z = 1
    elif Z is not None and Z_from_arg is not None and Z != Z_from_arg:
        raise InvalidIonError("Inconsistent charge information in"
                              "ion_mass.")
    elif Z is None and Z_from_arg is not None:
        Z = Z_from_arg

    if isinstance(Z, str) and Z.isdigit():
        Z = int(Z)
    if isinstance(mass_numb, str) and mass_numb.isdigit():
        mass_numb = int(mass_numb)

    if not isinstance(Z, int):
        raise TypeError("In ion_mass, Z must be an integer representing the "
                        "ionization state (e.g., Z=1 for singly ionized).")

    if not isinstance(mass_numb, int) and mass_numb is not None:
        raise TypeError("In ion_mass, mass_numb must be an integer "
                        "representing the mass number of an isotope.")

    if atomic_number(argument) < Z:
        raise InvalidParticleError("The ionization state cannot exceed the "
                                   "atomic number in ion_mass")

    try:
        isotope = isotope_symbol(argument, mass_numb)
    except InvalidParticleError as e:
        raise InvalidParticleError("Invalid particle in ion_mass.")
    except InvalidIsotopeError:
        is_isotope = False
    else:
        is_isotope = True

    if is_isotope:

        if isotope == 'H-1' and Z == 1:
            return const.m_p
        elif isotope == 'D' and Z == 1:
            return 3.343583719e-27 * u.kg
        elif isotope == 'T' and Z == 1:
            return 5.007356665e-27 * u.kg

        atomic_mass = isotope_mass(isotope)

    else:

        try:
            atomic_mass = standard_atomic_weight(argument)
        except MissingAtomicDataError:  # coveralls: ignore

            errormessage = ("No isotope mass or standard atomic weight is "
                            f"available to get ion mass for {argument}")

            if isinstance(mass_numb, int):
                errormessage += f" with mass number {mass_numb}"

            raise

    m_i = (atomic_mass - Z * const.m_e).to(u.kg)

    return m_i


def isotopic_abundance(argument: Union[str, int],
                       mass_numb: int = None) -> Quantity:
    r"""Returns the isotopic abundances if known, and otherwise zero.

    Parameters
    ----------

    argument: string or integer
        A string representing an element or isotope, or an integer
        representing the atomic number of an element.

    mass_numb: integer
        The mass number of an isotope, which is required if and only
        if the first argument can only be used

    Returns
    -------

    iso_comp: float
        The relative isotopic abundance in the terrestrial environment

    Raises
    ------

    InvalidIsotopeError
        If the argument is a valid particle but not a valid isotope.

    InvalidParticleError
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    TypeError
        If the argument is not a string or integer.


    Notes
    -----

    Isotopic composition data are most readily available for the
    terrestrial environment, so this function may not be wholly
    appropriate for space and astrophysical applications.

    The data retrieved from this routine are those recommended by NIST
    as of 2017.

    Examples
    --------

    >>> isotopic_abundance('Pb-208')
    0.524
    >>> isotopic_abundance('hydrogen', 1)
    0.999885
    >>> isotopic_abundance(118, 294)  # Og-294
    0.0

    """

    if _is_neutron(argument):
        raise InvalidIsotopeError(
            "Neutrons do not have an isotopic abundance.")

    try:
        isotope = isotope_symbol(argument, mass_numb)
    except InvalidParticleError:
        raise InvalidParticleError("Invalid particle in isotopic_abundance.")
    except InvalidIsotopeError:
        raise InvalidIsotopeError("Invalid isotope in isotopic_abundance.")

    iso_comp = _Isotopes[isotope].get('isotopic_abundance', 0.0)

    return iso_comp


def integer_charge(particle: str) -> int:
    r"""Returns the integer charge of an ion or other particle.

    Parameters
    ----------

    particle : string
        String representing a particle.

    Returns
    -------

    Z : integer
        The integer charge, or None if it is not available.

    Raises
    ------

    InvalidParticleError
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    ChargeError
        If charge information for the particle is not available.

    AtomicWarning
        If the input represents an ion with an integer charge that is
        below -3.

    Notes
    -----

    This function supports two formats for integer charge information.

    The first format is a string that has information for the element
    or isotope at the beginning, a space in between, and the integer
    charge information in the form of an integer followed by a plus or
    minus sign, or a plus or minus sign followed by an integer.

    The second format is a string containing element information at
    the beginning, following by one or more plus or minus signs.

    This function returns -1 for electrons, +1 for positrons, and 0
    for neutrons.

    Examples
    --------

    >>> integer_charge('Fe-56 2+')
    2
    >>> integer_charge('He -2')
    -2
    >>> integer_charge('H+')
    1
    >>> integer_charge('N-14++')
    2

    """

    if _is_electron(particle) or _is_antiproton(particle):
        return -1
    elif _is_positron(particle):
        return 1
    elif _is_neutron(particle) or _is_antineutron(particle):
        return 0

    particle, Z = _extract_integer_charge(particle)

    try:
        atomic_symbol(particle)
    except InvalidParticleError:
        raise InvalidParticleError("Invalid particle in integer_charge")

    if Z is None:
        raise ChargeError(f"Unable to find charge of {particle}.")

    atomic_numb = atomic_number(particle)

    if Z > atomic_numb:
        raise InvalidParticleError("The integer charge cannot be greater than "
                                   "the atomic number in integer_charge.")

    if Z < -atomic_numb - 1 or Z < -3:
        warnings.warn(f"Element {atomic_symbol(particle)} has a charge of {Z}"
                      " which is unlikely to occur in nature.", AtomicWarning)

    return Z


def electric_charge(particle: str) -> Quantity:
    r"""Returns the electric charge (in coulombs) of an ion or other
    particle

    Parameters
    ----------

    particle : string
        String representing an element or isotope followed by integer
        charge information.

    Returns
    -------

    charge: Quantity
        The electric charge in coulombs.

    Raises
    ------

    InvalidParticleError
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    ChargeError
        If charge information for the particle is not available.

    AtomicWarning
        If the input represents an ion with an integer charge that is
        below -3.

    Notes
    -----

    This function supports two formats for integer charge information.

    The first format is a string that has information for the element
    or isotope at the beginning, a space in between, and the integer
    charge information in the form of an integer followed by a plus or
    minus sign, or a plus or minus sign followed by an integer.

    The second format is a string containing element information at
    the beginning, following by one or more plus or minus signs.

    This function returns -1.6021766208e-19 C for electrons and
    1.6021766208e-19 C for positrons.

    Examples
    --------

    >>> electric_charge('p')
    <Quantity 1.60217662e-19 C>
    >>> electric_charge('e')
    <Quantity -1.60217662e-19 C>

    """

    try:
        charge = integer_charge(particle) * const.e.to('C')
    except InvalidParticleError:
        raise InvalidParticleError("Invalid particle in electric_charge.")
    except ChargeError:
        raise ChargeError("Charge information is not being provided to "
                          "electric_charge.")
    except TypeError:
        raise TypeError("Invalid type in electric_charge.")

    return charge


def is_isotope_stable(argument: Union[str, int],
                      mass_numb: int = None) -> bool:
    r"""Returns true for stable isotopes and false otherwise.

    Parameters
    ----------

    argument: integer or string
        A string representing an isotope or an integer representing an
        atomic number

    mass_numb: integer
        The mass number of the isotope.

    Returns
    -------

    is_stable: boolean
        True if the isotope is stable, False if it is unstable.

    Raises
    ------

    InvalidIsotopeError
        If the arguments correspond to a valid particle but not a
        valid isotope.

    InvalidParticleError
        If the arguments do not correspond to a valid particle.

    TypeError
        If the argument is not a string or integer.

    MissingAtomicDataError
        If stability information is not available.

    Examples
    --------

    >>> is_isotope_stable("H-1")
    True
    >>> is_isotope_stable("tritium")
    False

    """

    try:
        isotope = isotope_symbol(argument, mass_numb)
        is_stable = _Isotopes[isotope]['is_stable']
    except InvalidIsotopeError:
        raise InvalidIsotopeError("Invalid isotope in is_isotope_stable")
    except InvalidParticleError:
        raise InvalidParticleError("Invalid particle in is_isotope_stable")
    except KeyError:
        raise MissingAtomicDataError(f"No data on stability of {isotope}.")
    except TypeError:
        raise TypeError("The argument to is_isotope_stable must be a string.")

    return is_stable


def half_life(argument: Union[int, str], mass_numb: int = None) -> Quantity:
    r"""Returns the half-life in seconds for unstable isotopes, and
    numpy.inf for stable isotopes.

    Parameters
    ----------

    argument: integer or string
        A string representing an isotope or an integer representing an
        atomic number

    mass_numb: integer
        The mass number of the isotope.

    Returns
    -------

    half_life_sec: astropy Quantity
        The half-life in units of seconds.

    Raises:
    -------

    InvalidIsotopeError
        If the argument is a valid particle but not a valid isotope.

    InvalidParticleError
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    MissingAtomicDataError
        If no half-life data is available for the isotope.

    TypeError
        The argument is not an integer or string or the mass number is
        not an integer.

    AtomicWarning
        The half-life is unavailable so the routine returns None.

    Notes:
    ------

    At present there is limited half-life data available.

    Examples:
    ---------

    >>> half_life('T')
    <Quantity 3.888e+08 s>
    >>> half_life('n')
    <Quantity 881.5 s>
    >>> half_life('H-1')
    <Quantity inf s>

    """

    try:

        isotope = isotope_symbol(argument, mass_numb)

        if _Isotopes[isotope]['is_stable']:
            half_life_sec = np.inf * u.s
        else:
            half_life_sec = _Isotopes[isotope].get('half_life', None)

    except InvalidParticleError:
        raise InvalidParticleError("Invalid element in isotope_symbol.")
    except InvalidIsotopeError:
        raise InvalidIsotopeError(
            "Cannot determine isotope information from these inputs to "
            f"half_life: {argument}, {mass_numb}")
    except TypeError:
        raise TypeError("Incorrect argument type for half_life")

    if half_life_sec is None:
        warnings.warn(f"The half-life for isotope {isotope} is not"
                      "available; returning None.", MissingAtomicDataWarning)

    return half_life_sec


def known_isotopes(argument: Union[str, int] = None) -> List[str]:
    r"""Returns a list of all known isotopes of an element, or a list
    of all known isotopes of every element if no input is provided.

    Parameters
    ----------

    argument: integer or string, optional
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    Returns
    -------

    isotopes_list: list of strings or empty list
        List of all of the isotopes of an element that have been
        discovered, sorted from lowest mass number to highest mass
        number.  If no argument is provided, then a list of all known
        isotopes of every element will be returned that is sorted by
        atomic number, with entries for each element sorted by mass
        number.

    Raises
    ------

    InvalidElementError
        If the argument is a valid particle but not a valid element.

    InvalidParticleError
        If the argument does not correspond to a valid particle.

    TypeError
        If the argument is not a string or integer.

    Notes
    -----

    This list returns both natural and artifically produced isotopes.

    See also
    --------

    common_isotopes : returns isotopes with non-zero isotopic
        abundances

    stable_isotopes : returns isotopes that are stable against
        radioactive decay

    Examples
    --------
    >>> known_isotopes('H')
    ['H-1', 'D', 'T', 'H-4', 'H-5', 'H-6', 'H-7']
    >>> known_isotopes('helium 1+')
    ['He-3', 'He-4', 'He-5', 'He-6', 'He-7', 'He-8', 'He-9', 'He-10']
    >>> known_isotopes()[0:10]
    ['H-1', 'D', 'T', 'H-4', 'H-5', 'H-6', 'H-7', 'He-3', 'He-4', 'He-5']
    >>> len(known_isotopes())
    3352

    """

    def known_isotopes_for_element(argument):
        element = atomic_symbol(argument)
        isotopes = []
        for isotope in _Isotopes.keys():
            if element + '-' in isotope and isotope[0:len(element)] == element:
                isotopes.append(isotope)
        if element == 'H':
            isotopes.insert(1, 'D')
            isotopes.insert(2, 'T')
        mass_numbers = [mass_number(isotope) for isotope in isotopes]
        sorted_isotopes = [mass_number for (isotope, mass_number) in
                           sorted(zip(mass_numbers, isotopes))]
        return sorted_isotopes

    if argument is not None:
        try:
            element = atomic_symbol(argument)
            isotopes_list = known_isotopes_for_element(element)
        except InvalidElementError:
            raise InvalidElementError("known_isotopes is unable to get "
                                      f"isotopes from an input of: {argument}")
        except InvalidParticleError:
            raise InvalidParticleError("Invalid particle in known_isotopes.")
    elif argument is None:
        isotopes_list = []
        for atomic_numb in range(1, 119):
            isotopes_list += known_isotopes_for_element(atomic_numb)

    return isotopes_list


def common_isotopes(argument: Union[str, int] = None,
                    most_common_only: bool = False) -> List[str]:
    r"""Returns a list of isotopes of an element with an isotopic
    abundances greater than zero, or if no input is provided, a list
    of all such isotopes for every element.

    Parameters
    ----------

    argument: integer or string, optional
        A string or integer representing an atomic number or element,
        or a string represnting an isotope.

    most_common_only: boolean
        If set to True, return only the most common isotope

    Returns
    -------

    isotopes_list: list of strings or empty list
        List of all isotopes of an element with isotopic abundances
        greater than zero, sorted from most abundant to least
        abundant.  If no isotopes have isotopic abundances greater
        than zero, this function will return an empty list.  If no
        arguments are provided, then a list of all common isotopes of
        all elements will be provided that is sorted by atomic number,
        with entries for each element sorted from most abundant to
        least abundant.

    Raises
    ------

    InvalidElementError
        If the argument is a valid particle but not a valid element.

    InvalidParticleError
        If the argument does not correspond to a valid particle.

    TypeError
        If the argument is not a string or integer.

    Notes
    -----

    The isotopic abundances are based on the terrestrial environment
    and may not be wholly appropriate for space and astrophysical
    applications.

    See also
    --------

    known_isotopes : returns a list of isotopes that have been
        discovered

    stable_isotopes : returns isotopes that are stable against
        radioactive decay

    isotopic_abundance : returns the relative isotopic abundance

    Examples
    --------

    >>> common_isotopes('H')
    ['H-1', 'D']
    >>> common_isotopes(44)
    ['Ru-102', 'Ru-104', 'Ru-101', 'Ru-99', 'Ru-100', 'Ru-96', 'Ru-98']
    >>> common_isotopes('beryllium 2+')
    ['Be-9']
    >>> common_isotopes('Fe')
    ['Fe-56', 'Fe-54', 'Fe-57', 'Fe-58']
    >>> common_isotopes('Fe', most_common_only=True)
    ['Fe-56']
    >>> common_isotopes()[0:7]
    ['H-1', 'D', 'He-4', 'He-3', 'Li-7', 'Li-6', 'Be-9']

    """

    def common_isotopes_for_element(argument: Union[str, int],
                                    most_common_only: Optional[bool]) \
            -> List[str]:

        isotopes = known_isotopes(argument)
        CommonIsotopes = [isotope for isotope in isotopes if
                          'isotopic_abundance' in _Isotopes[isotope].keys()]
        isotopic_abundances = [_Isotopes[isotope]['isotopic_abundance']
                               for isotope in CommonIsotopes]
        sorted_isotopes = [iso_comp for (isotope, iso_comp) in
                           sorted(zip(isotopic_abundances, CommonIsotopes))]

        sorted_isotopes.reverse()

        if most_common_only and len(sorted_isotopes) > 1:
            sorted_isotopes = sorted_isotopes[0:1]

        return sorted_isotopes

    if argument is not None:

        try:
            element = atomic_symbol(argument)
            isotopes_list = \
                common_isotopes_for_element(element, most_common_only)
        except InvalidParticleError:
            raise InvalidParticleError("Invalid particle")
        except InvalidElementError:
            raise InvalidElementError(
                "common_isotopes is unable to get isotopes "
                f"from an input of: {argument}")

    elif argument is None:
        isotopes_list = []
        for atomic_numb in range(1, 119):
            isotopes_list += \
                common_isotopes_for_element(atomic_numb, most_common_only)

    return isotopes_list


def stable_isotopes(argument: Union[str, int] = None,
                    unstable: bool = False) -> List[str]:
    r"""Returns a list of all stable isotopes of an element, or if no
    input is provided, a list of all such isotopes for every element.

    Parameters
    ----------

    argument: integer or string
        A string or integer representing an atomic number or element,
        or a string represnting an isotope.

    unstable: boolean
        If set to True, this function will return a list of the
        unstable isotopes instead of the stable isotopes.

    Returns
    -------

    StableIsotopes: list of strings or empty list
        List of all stable isotopes of an element, sorted from lowest
        mass number.  If an element has no stable isotopes, this
        function returns an empty list.

    Raises
    ------

    InvalidElementError
        If the argument is a valid particle but not a valid element.

    InvalidParticleError
        If the argument does not correspond to a valid particle.

    TypeError
        If the argument is not a string or integer.

    Notes
    -----

    There are 254 isotopes for which no radioactive decay has been
    observed.  It is possible that some isotopes will be discovered to
    be unstable but with extremely long half-lives.  For example,
    bismuth-209 was recently discovered to have a half-life of about
    1.9e19 years.  However, such isotopes can be regarded as virtually
    stable for most applications.

    See also
    --------

    known_isotopes : returns a list of isotopes that have been
        discovered

    common_isotopes : returns isotopes with non-zero isotopic
        abundances

    Examples
    --------

    >>> stable_isotopes('H')
    ['H-1', 'D']
    >>> stable_isotopes(44)
    ['Ru-96', 'Ru-98', 'Ru-99', 'Ru-100', 'Ru-101', 'Ru-102', 'Ru-104']
    >>> stable_isotopes('beryllium')
    ['Be-9']
    >>> stable_isotopes('Pb-209')
    ['Pb-204', 'Pb-206', 'Pb-207', 'Pb-208']
    >>> stable_isotopes(118)
    []

    Find unstable isotopes

    >>> stable_isotopes('U', unstable=True)[:5] # only first five
    ['U-217', 'U-218', 'U-219', 'U-220', 'U-221']

    """

    def stable_isotopes_for_element(argument: Union[str, int],
                                    stable_only: Optional[bool]) -> List[str]:
        KnownIsotopes = known_isotopes(argument)
        StableIsotopes = [isotope for isotope in KnownIsotopes if
                          _Isotopes[isotope]['is_stable'] == stable_only]
        return StableIsotopes

    if argument is not None:
        try:
            element = atomic_symbol(argument)
            isotopes_list = \
                stable_isotopes_for_element(element, not unstable)
        except InvalidParticleError:
            raise InvalidParticleError("Invalid particle in stable_isotopes")
        except InvalidElementError:
            raise InvalidElementError(
                "stable_isotopes is unable to get isotopes "
                f"from an input of: {argument}")
    elif argument is None:
        isotopes_list = []
        for atomic_numb in range(1, 119):
            isotopes_list += \
                stable_isotopes_for_element(atomic_numb, not unstable)

    return isotopes_list
