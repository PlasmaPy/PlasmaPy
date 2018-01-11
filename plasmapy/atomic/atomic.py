"""Functions that retrieve or are related to elemental or isotopic data."""

import numpy as np
import re
from astropy import units as u, constants as const
from astropy.units import Quantity
import warnings
from .elements import atomic_symbols, atomic_symbols_dict, Elements
from .isotopes import Isotopes
from .particles import _Particles, _get_standard_symbol
from ..utils import (AtomicWarning, ElementError, IsotopeError, IonError,
                     ChargeError, AtomicError, MissingAtomicDataError)
from typing import (Union, Optional, Any, List, Tuple)

# The code contained within atomic_symbol(), isotope_symbol(), and
# _extract_charge_state() is designed to catch all of the special
# cases for different inputs.  Complexity is concentrated in these
# functions so that the rest of the functions can be simpler.

# TODO: Create an ion_symbol function
# TODO: Create a particle_symbol function
# TODO: Create a particle_mass function
# TODO: Create lepton_number and baryon_number functions
# TODO: Maybe create is_antimatter, is_lepton, is_baryon, is_boson, is_fermion


def atomic_symbol(argument: Union[str, int]) -> str:
    r"""Returns the atomic symbol.

    Parameters
    ----------

    argument: string or integer
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    Returns
    -------

    symbol: string
        The atomic symbol of the element, isotope, or nucleon.

    Raises
    ------

    ElementError
        If the argument does not correspond to a valid element.

    IsotopeError
        If the argument does not correspond to a valid isotope.

    IonError
        If the argument does not correspond to a valid ion.

    ChargeError
        If the argument has invalid charge information.

    TypeError
        If the argument is not a string or integer.

    See also
    --------

    isotope_symbol : returns isotope symbol instead of atomic symbol.

    element_name : returns the name of an element.

    Notes
    -----

    This function returns the symbol of the element rather than the
    symbol of an isotope.  For example, 'deuterium', 'T', or
    'hydrogen-2' will yield 'H'; 'alpha' will yield 'He'; and
    'iron-56' or 'Fe-56' will yield 'Fe'.

    This function is case insensitive when there is no ambiguity
    associated with case.  However, this function will return 'H' for
    hydrogen for lower case 'p' but capital 'P' if the argument is 'P'
    for phosphorus.  This function will return 'N' for nitrogen if the
    argument is capital 'N', but will not accept lower case 'n' for
    neutrons.

    Examples
    --------

    >>> atomic_symbol('helium')
    'He'
    >>> atomic_symbol(42)
    'Mo'
    >>> atomic_symbol('D')
    'H'
    >>> atomic_symbol('C-13')
    'C'
    >>> atomic_symbol('alpha')
    'He'
    >>> atomic_symbol('79')
    'Au'
    >>> atomic_symbol('N')  # Nitrogen
    'N'
    >>> atomic_symbol('P'), atomic_symbol('p')  # Phosphorus, proton
    ('P', 'H')

    """

    if _is_neutron(argument):
        raise ElementError("Neutrons do not have an atomic symbol")
    elif _is_antiproton(argument):
        raise ElementError("Antiprotons do not have an atomic symbol")

    try:
        argument, Z = _extract_charge_state(argument)
    except ChargeError:
        raise ChargeError("Invalid charge in atomic_symbol")

    if not isinstance(argument, (str, int)):
        raise TypeError("The first argument in atomic_symbol must be either "
                        "a string representing an element or isotope, or an "
                        "integer representing the atomic number (or 0 for "
                        "neutrons).")

    if isinstance(argument, str) and argument.isdigit():
        argument = int(argument)

    if isinstance(argument, int):

        try:
            element = atomic_symbols[argument]
        except KeyError:
            raise ElementError(f"{argument} is an invalid atomic number in "
                               "atomic_symbol")

    elif _is_hydrogen(argument):
        element = 'H'
    elif _is_alpha(argument):
        element = 'He'
    elif isinstance(argument, str):

        if argument.count('-') == 1:
            dash_position = argument.find('-')
            mass_numb = argument[dash_position+1:]
            if not mass_numb.isdigit():
                raise IsotopeError("Invalid isotope format in atomic_symbol")
            argument = argument[:dash_position]
        else:
            mass_numb = ''

        if argument.lower() in atomic_symbols_dict.keys():
            element = atomic_symbols_dict[argument.lower()]
        elif argument in atomic_symbols.values():
            element = argument.capitalize()
        else:
            raise ElementError(f"{argument} is an invalid argument for "
                               "atomic_symbol")

        if mass_numb.isdigit():

            isotope = element.capitalize() + '-' + mass_numb

            if isotope not in Isotopes.keys():
                raise IsotopeError("The input in atomic_symbol corresponding "
                                   f"to {isotope} is not a valid isotope.")

    if Z is not None and \
            Z > Elements[element]['atomic_number']:
        raise IonError("Cannot have an ionization state greater than the "
                       "atomic number.")

    return element


def isotope_symbol(argument: Union[str, int], mass_numb: int = None) -> str:
    r"""Returns the symbol representing an isotope.

    Parameters
    ----------

    argument: integer or string
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    mass_numb: integer or string
        An integer or string representing the mass number of the
        isotope.

    Returns
    -------

    symbol: string
        The isotopic symbol. The result will generally be returned as
        something like 'He-4' or 'Au-197', but will return 'D' for
        deuterium and 'T' for tritium.

    Raises
    ------

    ElementError
        If the argument does not correspond to a valid element.

    IsotopeError
        If the argument does not correspond to a valid isotope, or
        there is contradictory isotope information.

    IonError
        If the argument does not correspond to a valid ion.

    ChargeError
        If the argument has invalid charge information.

    TypeError
        If the argument is not a string or integer.

    Warns
    -----

    AtomicWarning
        If redundant isotope information is provided.

    See also
    --------
    atomic_symbol : returns atomic symbol instead of isotopic symbol

    Notes
    -----
    This function returns the symbol of the element rather than the
    symbol of an isotope.  For example, 'deuterium', 'T', or
    'hydrogen-2' will yield 'H'; 'alpha' will yield 'He'; and
    'iron-56' or 'Fe-56' will yield 'Fe'.

    Examples
    --------
    >>> isotope_symbol('He', 4)
    'He-4'
    >>> isotope_symbol(79, 197)
    'Au-197'
    >>> isotope_symbol('hydrogen-2')
    'D'
    >>> isotope_symbol('carbon-13')
    'C-13'
    >>> isotope_symbol('alpha')
    'He-4'

    """

    # If the argument is already in our standard form for an isotope,
    # return the argument.

    if mass_numb is None and argument in Isotopes.keys():
        return argument

    if isinstance(argument, str):
        argument, charge_state = _extract_charge_state(argument)

    if isinstance(argument, str) and argument.isdigit():
        argument = int(argument)

    if isinstance(mass_numb, str) and mass_numb.isdigit():
        mass_numb = int(mass_numb)

    # This routine allows several forms of input, and must be able to handle
    # all of the exceptions that can arise with useful error messages.

    if not isinstance(argument, (str, int)):
        raise TypeError("The first argument in isotope_symbol must be either "
                        "a string representing an element or isotope, or an "
                        "integer representing the atomic number (or 0 for "
                        " neutrons).")

    if not (isinstance(mass_numb, int) or
            mass_numb is None):  # coveralls: ignore
        raise TypeError("The second argument in isotope_symbol must be an "
                        "integer (or a string containing an integer) that "
                        "represents the mass number of an isotope.")

    if isinstance(argument, int):
        if not 0 <= argument <= 118:
            raise ElementError("Invalid atomic number in isotope_symbol")
        if mass_numb is None:
            raise IsotopeError("Insufficient information to determine "
                               "isotope in isotope_symbol.")

    # TODO: Remove this functionality when particle_symbol comes online
    if _is_neutron(argument, mass_numb):
        return 'n'

    try:
        element = atomic_symbol(argument)
    except ElementError:
        raise ElementError(f"The argument {argument} to isotope_symbol does "
                           f"not correspond to a valid element.")
    except IsotopeError:
        raise IsotopeError(f"The argument {argument} to isotope_symbol does "
                           f"not correspond to a valid isotope.")

    # Get mass number from argument, check for redundancies, and take
    # care of special cases.

    if isinstance(argument, str):
        if argument.count('-') == 1:
            dash_position = argument.find('-')
            mass_numb_from_arg = argument[dash_position+1:].strip()
            mass_numb_from_arg = int(mass_numb_from_arg)
        elif argument in ['p', 'p+'] or \
                argument.lower() in ['protium', 'proton']:
            mass_numb_from_arg = 1
        elif argument.lower() in ['d', 'deuterium', 'deuteron']:
            mass_numb_from_arg = 2
        elif argument.lower() in ['t', 'tritium', 'triton']:
            mass_numb_from_arg = 3
        elif _is_alpha(argument):
            mass_numb_from_arg = 4
        else:
            mass_numb_from_arg = None

        if mass_numb is None and mass_numb_from_arg is None:
            raise IsotopeError("Insufficient information to determine the mass"
                               " number from the inputs to isotope_symbol.")

        if mass_numb is not None and mass_numb_from_arg is not None:
            if mass_numb == mass_numb_from_arg:
                warnings.warn(
                    "Redundant mass number information in isotope_symbol "
                    f"from inputs: {argument}, {mass_numb}", AtomicWarning)
            else:  # coveralls: ignore
                raise IsotopeError("Contradictory mass number information in "
                                   "isotope_symbol.")

        if mass_numb_from_arg is not None:
            mass_numb = mass_numb_from_arg

    isotope = f"{element}-{mass_numb}"

    if isotope == 'H-2':
        isotope = 'D'
    elif isotope == 'H-3':
        isotope = 'T'

    if atomic_number(element) > mass_numb:
        raise IsotopeError("The atomic number cannot exceed the mass number "
                           "in isotope_symbol.")

    return isotope


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

    ElementError
        If the argument does not correspond to a valid element.

    IsotopeError
        If the argument does not correspond to a valid isotope.

    IonError
        If the argument does not correspond to a valid ion.

    ChargeError
        If the argument has invalid charge information.

    TypeError
        If the argument is not a string or integer.

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
        atomic_numb = Elements[element]['atomic_number']
    except ElementError:
        raise ElementError("Unable to identify element in atomic_number.")
    except IsotopeError:
        raise IsotopeError("Invalid isotope in atomic_number.")
    except IonError:
        raise IonError("Invalid ion in atomic_number.")
    except ChargeError:
        raise ChargeError("Invalid charge information in atomic_number.")
    except TypeError:
        raise TypeError("Invalid type in atomic_number")

    return atomic_numb


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

    IsotopeError
        If isotope cannot be determined or is not known.

    Examples
    --------

    >>> is_isotope_stable("H-1")
    True
    >>> is_isotope_stable("tritium")
    False

    """

    try:
        isotope = isotope_symbol(argument, mass_numb)
    except ElementError:
        raise ElementError("Invalid element in is_isotope_stable")
    except IsotopeError:
        raise IsotopeError("Invalid isotope in is_isotope_stable")
    except IonError:
        raise IonError("Invalid ion in is_isotope_stable")
    except TypeError:
        raise TypeError("Invalid input to is_isotope_stable")

    try:
        is_stable = Isotopes[isotope]['is_stable']
    except KeyError:  # coveralls: ignore
        raise MissingAtomicDataError("No data on stability of " + isotope)

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

    MissingAtomicDataError
        If no half-life data is available for the isotope.

    TypeError
        The argument is not an integer or string.

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
    except ElementError:
        raise ElementError("Invalid element in isotope_symbol.")
    except IsotopeError:
        raise IsotopeError("Cannot determine isotope information from these " +
                           f"inputs to half_life: {argument}, {mass_numb}")
    except TypeError:
        raise TypeError("Incorrect argument type for half_life")

    try:
        if Isotopes[isotope]['is_stable']:
            half_life_sec = np.inf * u.s
        else:
            half_life_sec = Isotopes[isotope]['half_life']
    except KeyError:
        half_life_sec = None
        warnings.warn(f"The half-life for isotope {isotope} is not available; "
                      "returning None.", AtomicWarning)

    return half_life_sec


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

    ElementError
        If the argument does not correspond to a valid element.

    IsotopeError
        If the argument does not correspond to a valid isotope.

    IonError
        If the argument does not correspond to a valid ion.

    ChargeError
        If the argument has invalid charge information.

    TypeError
        The first argument is not a string.


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
        mass_numb = Isotopes[isotope]["mass_number"]
    except TypeError:
        raise TypeError("Incorrect type for mass_number input.")
    except ElementError:
        raise ElementError("Invalid element in mass_number.")
    except IsotopeError:
        raise IsotopeError("Invalid isotope in mass_number.")
    except IonError:
        raise IonError("Invalid ion in mass_number.")
    except ChargeError:
        raise ChargeError("Invalid charge information in mass_number.")

    return mass_numb


def element_name(argument: Union[str, int]) -> str:
    r"""Returns the name of an element.

    Parameters
    ----------

    argument : string or integer
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    Returns
    -------

    name : string
        The name of the element.

    See also
    --------

    atomic_symbol : returns the atomic symbol

    isotope_symbol : returns the symbol of an isotope

    Examples
    --------

    >>> element_name("H")
    'hydrogen'
    >>> element_name("T")
    'hydrogen'
    >>> element_name("alpha")
    'helium'
    >>> element_name(42)
    'molybdenum'
    >>> element_name("C-12")
    'carbon'

    """

    try:
        element = atomic_symbol(argument)
        name = Elements[element]["name"]
    except ElementError:
        raise ElementError("Invalid element in element_name.")
    except IsotopeError:
        raise IsotopeError("Invalid isotope in element_name.")
    except IonError:
        raise IonError("Invalid ion in element_name.")
    except ChargeError:
        raise ChargeError("Invalid charge information in element_name.")
    except TypeError:
        raise TypeError("Invalid input to element_name.")

    return name


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

    ElementError:
        If the argument cannot be used to identify an element; the
        argument represents an isotope, ion, or neutron; or no
        standard atomic weight is provided for an element.

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

    >>> from astropy import units
    >>> standard_atomic_weight("H")
    <Quantity 1.008 u>
    >>> # the following result accounts for small amount of deuterium
    >>> standard_atomic_weight("H").to(units.kg)
    <Quantity 1.67382335e-27 kg>
    >>> isotope_mass("H-1")
    <Quantity 1.00782503 u>
    >>> standard_atomic_weight(82)
    <Quantity 207.2 u>
    >>> standard_atomic_weight("lead")
    <Quantity 207.2 u>

    """

    argument, charge_state = _extract_charge_state(argument)

    if charge_state is not None and charge_state != 0:
        raise AtomicError("Use ion_mass to get masses of charged particles.")

    try:
        isotope = isotope_symbol(argument)
    except Exception:
        isotope = ''

    if _is_neutron(isotope):
        raise ElementError("Use isotope_mass('n') or plasmapy.constants.m_n "
                           "instead of standard_atomic_weight to get neutron "
                           "mass.")
    elif '-' in isotope or isotope in ['D', 'T']:
        raise AtomicError("Use isotope_mass to get masses of isotopes.")

    try:
        element = atomic_symbol(argument)
        atomic_weight = Elements[element]['atomic_mass']
    except KeyError:
        raise MissingAtomicDataError(
            f"No standard atomic weight is available for {element}.")
    except ElementError:
        raise ElementError("Invalid element in standard_atomic_weight.")

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

    IsotopeError
        Contradictory or insufficient isotope information is provided.

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

    argument, charge_state = _extract_charge_state(argument)

    if charge_state is not None and charge_state != 0:
        raise AtomicError("Use ion_mass instead of isotope_mass for masses of "
                          "charged particles.")

    try:
        isotope = isotope_symbol(argument, mass_numb)
        atomic_mass = Isotopes[isotope]['atomic_mass']
    except ElementError:
        raise ElementError("Invalid element in isotope_mass.")
    except IsotopeError:
        raise IsotopeError("Unable to identify isotope in isotope_mass.")
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

    IonError
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
        except Exception:
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
        raise ElementError("Use isotope_mass or m_n to get mass of neutron")

    if isinstance(argument, str):
        argument, Z_from_arg = _extract_charge_state(argument)
    else:
        Z_from_arg = None

    if Z is None and Z_from_arg is None:
        Z = 1
    elif Z is not None and Z_from_arg is not None and Z != Z_from_arg:
        raise IonError("Inconsistent charge state information in ion_mass")
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
        raise IonError("The ionization state cannot exceed the "
                       "atomic number in ion_mass")

    try:
        isotope = isotope_symbol(argument, mass_numb)
    except Exception:
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
        except Exception:  # coveralls: ignore

            errormessage = ("No isotope mass or standard atomic weight is "
                            f"available to get ion mass for {argument}")

            if isinstance(mass_numb, int):
                errormessage += f" with mass number {mass_numb}"

            raise

    m_i = (atomic_mass - Z * const.m_e).to(u.kg)

    return m_i


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

    ValueError
        The argument cannot be used to determine an element.

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
        for isotope in Isotopes.keys():
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
        except ElementError:
            raise ElementError("known_isotopes is unable to get isotopes from "
                               f"an input of: {argument}")
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

    ValueError
        The argument cannot be used to determine an element.

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
                                    most_common_only: Optional[bool])\
            -> List[str]:

        isotopes = known_isotopes(argument)
        CommonIsotopes = [isotope for isotope in isotopes if
                          'isotopic_abundance' in Isotopes[isotope].keys()]
        isotopic_abundances = [Isotopes[isotope]['isotopic_abundance']
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
        except Exception:
            raise ElementError("common_isotopes is unable to get isotopes "
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
                          Isotopes[isotope]['is_stable'] == stable_only]
        return StableIsotopes

    if argument is not None:
        try:
            element = atomic_symbol(argument)
            isotopes_list = \
                stable_isotopes_for_element(element, not unstable)
        except Exception:
            raise ElementError("stable_isotopes is unable to get isotopes "
                               f"from an input of: {argument}")
    elif argument is None:
        isotopes_list = []
        for atomic_numb in range(1, 119):
            isotopes_list += \
                stable_isotopes_for_element(atomic_numb, not unstable)

    return isotopes_list


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
    ValueError
        Invalid isotope input, or the input corresponded to neutrons.

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

    try:
        element = atomic_symbol(argument)
        isotope = isotope_symbol(argument, mass_numb)
    except ElementError:
        raise ElementError("Invalid element in isotopic_abundance.")
    except IsotopeError:
        raise IsotopeError("Invalid isotope in isotopic_abundance.")
    except Exception:
        raise

    if isotope == 'n':
        raise IsotopeError("Neutrons do not have an isotopic abundance.")

    try:
        iso_comp = Isotopes[isotope]['isotopic_abundance']
    except KeyError:
        iso_comp = 0.0

    return iso_comp


def charge_state(particle: str) -> int:
    r"""Returns the charge state of an ion or other particle.

    Parameters
    ----------

    particle : string
        String representing a particle.

    Returns
    -------

    Z : integer
        The charge state, or None if it is not available.

    Raises
    ------

    ValueError
        If the charge state or isotope information is invalid, or the
        charge state exceeds the atomic number.

    AtomicWarning
        If the input represents an ion with a charge state that is
        below -3.

    Notes
    -----

    This function supports two formats for the charge state
    information.

    The first format is a string that has information for the element
    or isotope at the beginning, a space in between, and the charge
    state information in the form of an integer followed by a plus or
    minus sign, or a plus or minus sign followed by an integer.

    The second format is a string containing element information at
    the beginning, following by one or more plus or minus signs.

    This function returns -1 for electrons, +1 for positrons, and 0
    for neutrons.

    Examples
    --------

    >>> charge_state('Fe-56 2+')
    2
    >>> charge_state('He -2')
    -2
    >>> charge_state('H+')
    1
    >>> charge_state('N-14++')
    2

    """

    if _is_electron(particle) or _is_antiproton(particle):
        return -1
    elif _is_positron(particle):
        return 1
    elif _is_neutron(particle) or _is_antineutron(particle):
        return 0

    particle, Z = _extract_charge_state(particle)

    try:
        atomic_numb = atomic_number(particle)
    except ElementError:
        raise ElementError("Invalid element in charge_state")
    except IsotopeError:
        raise IsotopeError("Invalid isotope in charge_state")
    except IonError:
        raise IonError("Invalid ion in charge_state")

    if Z is not None and Z > atomic_numb:
        raise IonError("The charge state cannot be greater than the atomic "
                       "number.")

    if Z is not None and (Z < -atomic_numb-1 or Z < -3):
        warnings.warn(f"Element {atomic_symbol(particle)} has a charge of {Z}"
                      " which is unlikely to occur in nature.", AtomicWarning)

    if Z is None:
        raise ChargeError(f"Unable to find charge of {particle}")

    return Z


def electric_charge(particle: str) -> Quantity:
    r"""Returns the electric charge (in coulombs) of an ion or other
    particle

    Parameters
    ----------

    particle : string
        String representing an element or isotope followed by charge
        state information.

    Returns
    -------

    charge: Quantity
        The electric charge in coulombs.

    Raises
    ------

    ValueError
        If the charge state or isotope information is invalid, or the
        charge state exceeds the atomic number.

    AtomicWarning
        If the input represents an ion with a charge state that is
        below -3.

    Notes
    -----

    This function supports two formats for the charge state
    information.

    The first format is a string that has information for the element
    or isotope at the beginning, a space in between, and the charge
    state information in the form of an integer followed by a plus or
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
        charge = charge_state(particle) * const.e.to('C')
        return charge
    except Exception:  # coveralls: ignore
        raise


def _extract_charge_state(argument: str) -> Tuple[str, int]:
    r"""Splits strings containing element or isotope and charge state
    information into a string without the charge state information and
    the charge state as an integer (or None if no charge state
    information is given).

    Parameters
    ----------

    argument : string
        String containing information for an element or isotope in any
        of the allowed formats, followed by charge state information.

    Returns
    -------

    argument : string
        The original string with charge state information removed.

    Z : integer
        The charge state of an ion (e.g., this will return 1 if one
        electron has been removed and -1 if one electron has been
        gained)

    Raises
    ------

    ChargeError
        If invalid charge information is included.

    Notes
    -----

    If the argument is not a string, this function will return the
    original argument and None.

    This function does not check that the charge state is valid.

    Examples
    --------

    >>> isotope, Z = _extract_charge_state('Fe-56+++')
    >>> print(isotope)
    Fe-56
    >>> print(Z)
    3
    >>> _extract_charge_state('D +1')
    ('D', 1)

    """

    if not isinstance(argument, str):
        return argument, None

    if argument in ['n', 'neutron', 'n-1']:
        return argument, 0
    elif argument in ['p', 'p+'] or argument.lower() in \
            ['proton', 'deuteron', 'triton']:
        return argument, 1
    elif argument.lower() == 'alpha':
        return argument, 2
    elif argument == 'e+' or argument.lower() == 'positron':
        return argument, 1
    elif _is_antiproton(argument):
        return argument, -1

    if argument.count(' ') == 1:  # For cases like 'Fe +2' and 'Fe-56 2+'

        argument, ion_info = argument.split()

        check1 = (ion_info.endswith(('-', '+')) ^
                  ion_info.startswith(('-', '+')))

        check2 = ((ion_info.count('-') == 1 and ion_info.count('+') == 0) or
                  (ion_info.count('-') == 0 and ion_info.count('+') == 1))

        if '-' in ion_info:
            sign = -1
            charge = ion_info.replace('-', '')
        elif '+' in ion_info:
            sign = +1
            charge = ion_info.replace('+', '')

        try:
            charge_state = sign*int(charge)
            check3 = True
        except Exception:
            check3 = False

        if not (check1 and check2 and check3):
            raise ChargeError(
                "The following input does not have valid charge: "
                f"state information: {argument}{ion_info}")

    elif argument.endswith(('-', '+')):  # For cases like 'Fe++' or 'Si-'

        char = argument[-1]
        match = re.match(r"["+char+"]*", argument[::-1])

        charge_state = match.span()[1]

        if char == '-':
            charge_state = -charge_state

        argument = argument[0:len(argument)-match.span()[1]]

        if argument.endswith(('-', '+')):
            raise ChargeError("Invalid charge state information")

    else:
        charge_state = None

    if charge_state is not None and charge_state < -3:
        warnings.warn(f"Element {atomic_symbol(argument)} has a charge of "
                      f"{charge_state} which is unlikely to occur in nature.",
                      AtomicWarning)

    return argument, charge_state


def _is_neutron(argument: Any, mass_numb: int = None) -> bool:
    r"""Returns True if the argument corresponds to a neutron, and
    False otherwise."""

    if argument == 0 and mass_numb == 1:
        return True
    elif isinstance(argument, str) and mass_numb is None:
        if argument in ('n', 'n-1') or argument.lower() in ['neutron', 'n0']:
            return True
        else:
            return False
    else:
        return False


def _is_hydrogen(argument: Any,
                 can_be_atomic_number: Optional[bool] = False) -> bool:
    r"""Returns True if the argument corresponds to hydrogen, and False
    otherwise."""

    case_sensitive_aliases = ['p', 'p+', 'H', 'D', 'T']

    case_insensitive_aliases = ['proton', 'protium', 'deuterium',
                                'deuteron', 'triton', 'tritium', 'hydrogen']

    for mass_numb in range(1, 8):
        case_sensitive_aliases.append(f'H-{mass_numb}')
        case_insensitive_aliases.append(f'hydrogen-{mass_numb}')

    if isinstance(argument, str):

        argument, Z = _extract_charge_state(argument)

        if argument in case_sensitive_aliases:
            is_hydrogen = True
        else:
            is_hydrogen = argument.lower() in case_insensitive_aliases

        if is_hydrogen and Z is not None and Z > 1:
            raise IonError("Invalid charge state of hydrogen")

    elif argument == 1 and can_be_atomic_number:
        is_hydrogen = True
    else:
        is_hydrogen = False

    return is_hydrogen


def _is_electron(arg: Any) -> bool:
    r"""Returns True if the argument corresponds to an electron, and False
    otherwise."""

    if not isinstance(arg, str):
        return False

    return arg in ['e', 'e-'] or arg.lower() == 'electron'


def _is_positron(arg: Any) -> bool:
    r"""Returns True if the argument corresponds to a positron, and False
    otherwise."""

    if not isinstance(arg, str):
        return False

    return arg == 'e+' or arg.lower() == 'positron'


def _is_antiproton(arg: Any) -> bool:
    r"""Returns True if the argument corresponds to an antiproton, and
    False otherwise."""

    if not isinstance(arg, str):
        return False

    return arg == 'p-' or arg.lower() == 'antiproton'


def _is_antineutron(arg: Any) -> bool:
    r"""Returns True if the argument corresponds to an antineutron, and
    False otherwise."""

    if not isinstance(arg, str):
        return False

    return arg.lower() == 'antineutron'


def _is_proton(arg: Any, Z: int = None, mass_numb: int = None) -> bool:
    r"""Returns True if the argument corresponds to a proton, and
    False otherwise.  This function returns False for 'H-1' if no
    charge state is given."""

    try:

        isotope = isotope_symbol(arg, mass_numb)

        if Z is None:
            Z = charge_state(arg)

        return isotope == 'H-1' and Z == 1

    except Exception:

        return False


def _is_alpha(arg: Any) -> bool:
    r"""Returns True if the argument corresponds to an alpha particle,
    and False otherwise."""

    if not isinstance(arg, str):
        return False

    if arg.lower() == 'alpha':
        return True
    else:
        arg, Z = _extract_charge_state(arg)

        if Z != 2 or arg[-2:] != '-4':
            return False
        else:

            dash_position = arg.find('-')
            arg = arg[:dash_position]

            return arg.lower() == 'helium' or arg == 'He'
