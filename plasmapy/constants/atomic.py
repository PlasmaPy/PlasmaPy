"""Functions that retrieve or are related to elemental or isotopic data."""

import numpy as np
import re
from astropy import units as u, constants as const
from .elements import atomic_symbols_list, atomic_symbols_dict, Elements
from .isotopes import Isotopes


# The code contained within element_symbol(), isotope_symbol(), and
# __extract_charge_state() is designed to catch all of the special
# cases for different inputs.  Complexity is concentrated in these
# functions so that the rest of the functions can be simpler.

def element_symbol(argument):
    """Returns the atomic symbol.

    Parameters
    ----------
    argument: integer or string
        An integer representing the atomic number, or a string containing a
        representation of an element or an isotope.

    Returns
    -------
    symbol: string
        The atomic symbol of the element, isotope, or nucleon.

    Raises
    ------
    TypeError:
        If the argument is not a string or integer.

    ValueError:
        If the argument cannot be used to identify the element.

    See also
    --------
    isotope_symbol : returns isotope symbol instead of atomic symbol.

    element_name : returns the name of an element.

    Notes
    -----
    This function returns the symbol of the element rather than the symbol of
    an isotope.  For example, 'deuterium', 'T', or 'hydrogen-2' will yield 'H';
    'alpha' will yield 'He'; and 'iron-56' or 'Fe-56' will yield 'Fe'.

    This function is case insensitive when there is no ambiguity associated
    with case.  However, this function will return 'H' for hydrogen for lower
    case 'p' but capital 'P' if the argument is 'P' for phosphorus.  This
    function will also return lower case 'n' if the argument is lower case 'n'
    for neutrons, but capital 'N' for nitrogen if the argument is capital 'N'.

    Examples
    --------
    >>> element_symbol('helium')
    'He'
    >>> element_symbol(42)
    'Mo'
    >>> element_symbol('D')
    'H'
    >>> element_symbol('C-13')
    'C'
    >>> element_symbol('alpha')
    'He'
    >>> element_symbol('79')
    'Au'
    >>> element_symbol('N'), element_symbol('n')  # Nitrogen, neutron
    ('N', 'n')
    >>> element_symbol('P'), element_symbol('p')  # Phosphorus, proton
    ('P', 'p')

    """

    argument, charge_state = __extract_charge_state(argument)

    if not isinstance(argument, (str, int)):
        raise TypeError("The first argument in element_symbol must be either "
                        "a string representing an element or isotope, or an "
                        "integer representing the atomic number (or 0 for "
                        " neutrons).")

    if isinstance(argument, str) and argument.isdigit():
        argument = int(argument)

    if isinstance(argument, int):

        if 0 <= argument <= 118:
            symbol = atomic_symbols_list[argument]
        else:
            raise ValueError(str(argument)+" is an invalid atomic number in "
                             "element_symbol")

    elif isinstance(argument, str):

        if argument in ['n', 'neutron', 'n-1']:
            return 'n'
        elif argument in ['p', 'p+'] or argument.lower() in \
                ['d', 't', 'proton', 'protium', 'deuterium', 'deuteron',
                 'triton', 'tritium']:
            return 'H'
        elif argument.lower() == 'alpha':
            return 'He'

        if argument.count('-') == 1:
            dash_position = argument.find('-')
            mass_numb = argument[dash_position+1:]
            if not mass_numb.isdigit():
                raise ValueError("Invalid isotope format in element_symbol")
            argument = argument[:dash_position]
        else:
            mass_numb = ''

        if atomic_symbols_dict.keys().__contains__(argument.lower()):
            symbol = atomic_symbols_dict[argument.lower()]
        elif atomic_symbols_list.__contains__(argument.capitalize()):
            symbol = argument.capitalize()
        else:
            raise ValueError(argument+" is an invalid argument for "
                             "element_symbol")

        if mass_numb.isdigit():
            isotope = symbol.capitalize() + '-' + mass_numb
            if isotope == 'H-2':
                isotope = 'D'
            if isotope == 'H-3':
                isotope = 'T'

            if isotope not in Isotopes.keys():
                raise ValueError("The input in element_symbol corresponding "
                                 "to " + isotope + " is not a valid isotope.")

    if charge_state is not None and \
            charge_state > Elements[symbol]['atomic_number']:

        raise ValueError("Cannot have an ionization state greater than the "
                         "atomic number in element_name.")

    return symbol


def isotope_symbol(argument, mass_numb=None):
    """Returns the symbol representing an isotope.

    Parameters
    ----------
    argument: integer or string
        An integer representing the atomic number, or a string containing a
        representation of an element or an isotope.

    mass_numb: integer or string
        An integer or string representing the mass number of the isotope.

    Returns
    -------
    symbol: string
        The isotopic symbol. The result will generally be returned as something
        like 'He-4' or 'Au-197', but will return 'D' for deuterium and 'T' for
        tritium.

    Raises
    ------
    ValueError:
        If insufficient or contradictory isotope information is provided, the
        element cannot be determined from the first argument, or the mass
        number exceeds the atomic number.

    TypeError:
        If isotope information cannot be found because one or both inputs is
        of an inappropriate type.

    UserWarning:
        If redundant isotope information is provided.

    See also
    --------
    element_symbol : returns atomic symbol instead of isotopic symbol

    Notes
    -----
    This function returns the symbol of the element rather than the symbol of
    an isotope.  For example, 'deuterium', 'T', or 'hydrogen-2' will yield 'H';
    'alpha' will yield 'He'; and 'iron-56' or 'Fe-56' will yield 'Fe'.

    Examples
    --------
    >>> isotope_symbol('He, 4')
    'He-4'
    >>> isotope_symbol(79, 197)
    'Au-197
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
        argument, charge_state = __extract_charge_state(argument)

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

    if not (isinstance(mass_numb, int) or mass_numb is None):
        raise TypeError("The second argument in isotope_symbol must be an "
                        "integer (or a string containing an integer) that "
                        "represents the mass number of an isotope.")

    if isinstance(argument, int) and mass_numb is None:
        raise ValueError("Insufficient information to determine element and "
                         "mass number in isotope_symbol.")

    try:
        element = element_symbol(argument)
    except Exception:
        raise ValueError("The first argument of isotope_symbol (" +
                         str(argument) + ") does not correspond to a valid "
                         "element or isotope.")

    # Get mass number from argument, check for redundancies, and take
    # care of special cases.

    if isinstance(argument, str):
        if argument.count('-') == 1:

            dash_position = argument.find('-')
            mass_numb_from_arg = argument[dash_position+1:].strip()

            try:
                mass_numb_from_arg = int(mass_numb_from_arg)
            except Exception:
                raise ValueError("Unable to extract mass number from the first"
                                 " argument of isotope_symbol, which is: " +
                                 str(argument))

        elif argument == 'n' or argument.lower() == 'neutron':
            mass_numb_from_arg = 1
        elif argument in ['p', 'p+'] or \
                argument.lower() in ['protium', 'proton']:
            mass_numb_from_arg = 1
        elif argument.lower() in ['d', 'deuterium', 'deuteron']:
            mass_numb_from_arg = 2
        elif argument.lower() in ['t', 'tritium', 'triton']:
            mass_numb_from_arg = 3
        elif argument.lower() in ['alpha']:
            mass_numb_from_arg = 4
        else:
            mass_numb_from_arg = None

        if mass_numb is None and mass_numb_from_arg is None:
            raise ValueError("Insufficient information to determine the mass "
                             "number from the inputs to isotope_symbol.")

        if mass_numb is not None and mass_numb_from_arg is not None:
            if mass_numb == mass_numb_from_arg:
                raise UserWarning("Redundant mass number information in " +
                                  "isotope_symbol from inputs: (" +
                                  str(argument)+", " + str(mass_numb) + ")")
            else:
                raise ValueError("Contradictory mass number information in "
                                 "isotope_symbol.")

        if mass_numb_from_arg is not None:
            mass_numb = mass_numb_from_arg

    isotope = element + '-' + str(mass_numb)

    if isotope == 'n-1':
        isotope = 'n'
    elif isotope == 'H-2':
        isotope = 'D'
    elif isotope == 'H-3':
        isotope = 'T'

    if atomic_number(element) > mass_numb:
        raise ValueError("The atomic number cannot exceed the mass number in "
                         "isotope_symbol.")

    if isotope not in Isotopes.keys():
        raise ValueError("The isotope " + isotope + "returned by "
                         "isotope_symbol is unknown and may not exist.")

    return isotope


def atomic_number(argument):
    """Returns the atomic number (the number of protons in an atom)
    from an atomic symbol or name.

    Parameters
    ----------
    argument: string
        A string represnting the symbol or name of an element or isotope.

    Returns
    -------
    atomic_number: integer
        An integer representing the atomic number of the element or isotope.

    See also
    --------
    mass_number : returns the mass number (the total number of protons and
        neutrons) of an isotope.

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

    atomic_symbol = element_symbol(argument)
    atomic_numb = Elements[atomic_symbol]['atomic_number']
    return atomic_numb


def is_isotope_stable(argument, mass_numb=None):
    """Returns true for stable isotopes and false otherwise.

    Parameters
    ----------
    argument: integer or string
        A string or integer representing an atomic number or element, or a
        string represnting an isotope.

    mass_numb: integer
        The mass number of the isotope.

    Returns
    -------
    is_stable: boolean
        True if the isotope is stable, False if it is unstable.

    Raises
    ------
    ValueError:
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
    except (UserWarning, ValueError):
        raise ValueError("Unable to determine isotope in input for "
                         "is_isotope_stable")

    try:
        is_stable = Isotopes[isotope]['is_stable']
    except Exception:
        ValueError("No data on stability of " + isotope)

    return is_stable


def half_life(argument, mass_numb=None):
    """Returns the half-life in seconds for unstable isotopes, and numpy.inf
    for stable isotopes.

    Parameters
    ----------
    argument: integer or string
        A string or integer representing an atomic number or element, or a
        string represnting an isotope.
    mass_numb: integer
        The mass number of the isotope.

    Returns
    -------
    half_life_sec: astropy Quantity
        The half-life in units of seconds.

    Raises:
    -------
    ValueError:
        If no half-life data is available for the isotope.

    Notes:
    ------
    At present there is limited half-life data available.

    Examples:
    ---------
    >>> half_life('T')
    <Quantity 388800000.0 s>
    >>> half_life('n')
    <Quantity 881.5 s>
    >>> half_life('H-1')
    <Quantity inf s>

    """

    try:
        isotope = isotope_symbol(argument, mass_numb)
    except ValueError:
        raise ValueError("Cannot determine isotope information from these " +
                         "inputs to half_life: (" + str(argument) + ", " +
                         str(mass_numb) + ")")
    except TypeError:
        raise TypeError("Incorrect argument type for half_life")

    try:
        if Isotopes[isotope]['is_stable']:
            half_life_sec = np.inf * u.s
        else:
            half_life_sec = Isotopes[isotope]['half_life']
    except Exception:
        half_life_sec = None
        raise UserWarning("The half-life for isotope " + isotope +
                          " is not available; returning None.")

    return half_life_sec


def mass_number(isotope):
    """Get the mass number (the number of protons and neutrals) of an isotope.

    Parameters
    ----------
    isotope : string
        A string representing the symbol or name of an isotope.

    Returns
    -------
    mass_number : integer
       An integer representing the mass number of the isotope.

    Raises
    ------
    ValueError
        If the mass number cannot be found.

    See also
    --------
    atomic_number : returns the number of protons in an isotope or element

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
    >>> mass_number("N")
    7
    >>> mass_number("alpha")
    4

    """

    try:
        symbol = isotope_symbol(isotope)
        mass_numb = Isotopes[symbol]["mass_number"]
    except TypeError:
        raise("Incorrect type for mass_number input.")
    except ValueError:
        raise ValueError("Mass number not able to be found from input " +
                         str(isotope))

    return mass_numb


def element_name(argument):
    """Returns the name of an element.

    Parameters
    ----------
    argument : string or integer
        String representing an element or isotope, or an integer representing
        the atomic number of an element.

    Returns
    -------
    name : string
        The name of the element.

    See also
    --------
    element_symbol : returns the atomic symbol

    isotope_symbol : returns the symbol of an isotope

    Examples
    --------
    >>> element_name("H")
    'hydrogen'
    >>> element_name("T")
    'hydrogen'
    >>> element_name("alpha")
    'helium'
    >>> element_name()
    >>> element_name(42)
    'molybdenum'
    >>> element_name("C-12")
    'carbon'

    """

    atomic_symbol = element_symbol(argument)
    name = Elements[atomic_symbol]["name"]

    return name


def standard_atomic_weight(argument):
    """Returns the standard (conventional) atomic weight of an element based on
    the relative abundances of isotopes in terrestrial environments.

    Parameters
    ----------
    argument: string or integer
        A string representing the symbol or name of an element or isotope, or
        an integer representing an element's atomic number.

    Returns
    -------
    atomic_weight: astropy.units.Quantity with units of astropy.units.u
        The standard atomic weight of an element based on values frmo

    Raises
    ------
    ValueError:
        If the argument cannot be used to identify an element, or if the
        argument represents an isotope.

    UserWarning:
        If the standard atomic weight is not provided for an element (e.g., if
        an element does not naturally occur).

    See also
    --------
    isotope_mass : returns the atomic mass of an isotope.

    ion_mass : returns the mass of an ion of an element or isotope, accounting
        for reduction of mass from the neutral state due to different numbers
        of electrons.

    Notes
    -----
    The relative abundances of different isotopes of an element sometimes vary
    naturally in different locations within the terrestrial environment.  The
    CIAAW provides ranges for these element, which include H, Li, B, C, N, O,
    Mg, Si, S, Cl, Br, Tl.  This function provides a single value from the
    CIAWW 2015 standard values when a single value is given, and the lower
    accuracy conventional value given by Meija et al. (2013,
    doi:10.1515/pac-2015-0305) for the elements where a range is given.

    Examples
    --------
    >>> from astropy import units
    >>> standard_atomic_weight("H")
    <Quantity 1.008 u>
    >>> standard_atomic_weight("H").to(units.kg)
    <Quantity 1.673823232368e-27 kg>  # accounts for small amount of deuterium
    >>> isotope_mass("H-1")
    <Quantity 1.6735326915759943e-27 kg>  # pure hydrogen-1
    >>> standard_atomic_weight(82)
    <Quantity 238.02891 u>
    >>> standard_atomic_weight("lead")
    <Quantity 238.02891 u>
    >>> standard_atomic_weight(118) is None
    UserWarning: No standard atomic weight is available for Og. Returning None.
    True

    """

    argument, charge_state = __extract_charge_state(argument)

    if charge_state is not None and charge_state != 0:
        raise ValueError("Use ion_mass to get masses of charged particles.")

    try:
        isotope = isotope_symbol(argument)
    except Exception:
        isotope = ''

    if isinstance(argument, str) and (argument in ['p', 'p+'] or
                                      argument.lower() in ['proton', 'alpha']):
        raise ValueError("Use ion_mass to get masses of protons and alpha "
                         "particles instead of standard_atomic_weight")
    elif isotope == 'n':
        raise ValueError("Use isotope_mass('n') or plasmapy.constants.m_n "
                         "instead of standard_atomic_weight to get neutron "
                         "mass")
    elif '-' in isotope or isotope in ['D', 'T']:
        raise ValueError("Use isotope_mass to get masses of isotopes")

    atomic_symbol = element_symbol(argument)

    try:
        atomic_weight = Elements[atomic_symbol]['atomic_mass']
    except Exception:
        raise ValueError("No standard atomic weight is available for " +
                         atomic_symbol)

    return atomic_weight


def isotope_mass(argument, mass_numb=None):
    """Return the mass of an isotope.

    Parameters
    ----------
    argument : string or integer
        A string representing an element or isotope or an integer representing
        an atomic number

    mass_numb : integer (optional)
        The mass number of the isotope.

    Returns
    -------
    isotope_mass : Quantity
        The atomic mass of a neutral atom of an isotope.

    Raises
    ------
    UserWarning:
        Redundant (but consistent) isotope information is provided.

    ValueError:
        Contradictory or insufficient isotope information is provided.

    See also
    --------
    standard_atomic_weight : returns atomic weight of an element based on
        terrestrial abundances of isotopes

    ion_mass : returns the mass of an ion of an element or isotope, accounting
        for loss of electrons

    Notes
    -----
    The masses of rare isotopes may be unavailable.

    Examples
    --------
    >>>from astropy import units as u
    >>>isotope_mass("H-1")
    <Quantity 1.00782503223 u>
    >>>isotope_mass("H-1").to(units.kg)
    <Quantity 1.6735326915759943e-27 kg>
    >>> isotope_mass("He", 4)
    <Quantity 4.00260325413 u>
    >>> isotope_mass(2, 4)

    """

    argument, charge_state = __extract_charge_state(argument)

    if charge_state is not None and charge_state != 0:
        raise ValueError("Use ion_mass instead of isotope_mass for masses of "
                         "charged particles")

    if argument == 'alpha':
        raise ValueError("Use ion_mass for mass of an alpha particle")
    elif argument in ['p', 'proton']:
        raise ValueError("Use plasmapy.constants.m_p or ion_mass('p') to get "
                         "proton mass")

    symbol = isotope_symbol(argument, mass_numb)
    atomic_mass = Isotopes[symbol]['atomic_mass']

    return atomic_mass


def ion_mass(argument, Z=None, mass_numb=None):
    """Returns the mass of an ion by finding the standard atomic
    weight of an element or the atomic mass of an isotope, and then
    accounting for the change in mass due to loss of electrons from
    ionization.

    Parameters
    ----------
    argument: string or integer
        Returns the mass of an ion by finding the standard atomic
        weight of an element or the atomic mass of an isotope, and then
        accounting for the change in mass due to loss of electrons from
        ionization.

    Z: integer (optional)
        The ionization state of the ion (defaulting to a charge of Z=1)

    mass_numb: integer (optional)
        The mass number of an isotope.

    Returns
    -------
    m_i: Quantity
        The mass of a single ion of the isotope or element with charge state Z.

    Raises
    ------
    ValueError
        If the ionization state exceeds the atomic number.

    UserWarning
        If a mass was inputted and it is outside of the range of known isotopes
        or electrons/positrons.

    See also
    --------
    standard_atomic_mass : returns the conventional atomic mass of an element
        based on terrestrial values and assuming the atom is neutral.

    isotope_mass : returns the mass of an isotope (if available) assuming the
        atom is neutral.

    Notes
    -----
    This function in general finds the mass of an isotope (or the standard
    atomic weight based on terrestrial values if a unique isotope cannot be
    identified), and then substracts the mass of Z electrons.  If Z is not
    provided as an input, then this function assumes that the ion is singly
    ionized.

    Specific values are returns for protons, deuterons, tritons, alpha
    particles, and positrons.

    Calling ion_mass('H') does not return the mass of a proton but instea
    uses hydrogen's standard atomic weight based on terrestrial values.  To
    get the mass of a proton, use ion_mass('p').

    This function can accept a Quantity in units of mass.  If the Quantity
    is close to the mass of an electron or positron, it will return the
    mass of an electron or positron to full known precision.  If the
    Quantity is within the mass range of known isotopes, it will return the
    mass that is inputted to it.

    Examples
    --------
    >>> ion_mass('p')  # proton
    <Constant name='Proton mass' value=1.672621777e-27 uncertainty=7.4e-35
    unit='kg' reference='CODATA 2010'>
    >>> ion_mass('H')  # assumes terrestrial abundance of D
    <Quantity 1.672912294077e-27 kg>
    >>> ion_mass('H') == ion_mass('p')
    False
    >>> ion_mass('P')  # phosphorus
    <Quantity 5.143222638917872e-26 kg>
    >>> ion_mass('He-4', 2)
    <Quantity 6.64465723e-27 kg>
    >>> ion_mass('T')
    <Quantity 3.343583719e-27 kg>
    >>> ion_mass(26, Z=1, mass_numb=56)
    <Quantity 9.288122788133088e-26 kg>
    >>> ion_mass('Fe-56')
    <Quantity 9.288122788133088e-26 kg>
    >>> ion_mass(9.11e-31*u.kg)
    <Constant name='Electron mass' value=9.10938291e-31 uncertainty=4e-38 unit='kg' reference='CODATA 2010'>
    >>> ion_mass(1.67e-27*u.kg)
    <Quantity 1.67e-27 kg>

    """

    if isinstance(argument, u.Quantity) and Z is None and mass_numb is None:

        try:
            m_i = argument.to(u.kg)
        except:
            raise UnitConversionError("If the ion in given as a Quantity, "
                                      "then it must have units of mass.")

        if np.isclose(m_i.value, const.m_e.value, atol=1e-33):  # positrons
            return const.m_e
        elif 1.66e-27 <= m_i.value < 7e-25:  # mass range of known isotopes
            return m_i
        else:
            raise UserWarning("The mass that was inputted to ion_mass and is "
                              "being returned from ion_mass is outside of the "
                              "range of known isotopes or electrons/ions.")

    if isinstance(argument, str) and \
            str(argument).lower() in ['e+', 'positron']:
        return const.m_e

    if isinstance(argument, str) and \
            str(argument).lower() in ['e+', 'positron']:
        return const.m_e

    if atomic_number(argument) == 0:
        raise ValueError("Use isotope_mass or m_n to get mass of neutron")

    if isinstance(argument, str):
        argument, Z_from_arg = __extract_charge_state(argument)
    else:
        Z_from_arg = None

    if Z is None and Z_from_arg is None:
        Z = 1
    elif Z is not None and Z_from_arg is not None and Z != Z_from_arg:
        raise ValueError("Inconsistent charge state information in ion_mass")
    elif Z is None and Z_from_arg is not None:
        Z = Z_from_arg
    elif Z == Z_from_arg:
        raise UserWarning("Redundant charge state information in ion_mass")

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
        raise ValueError("The ionization state cannot exceed the "
                         "atomic number in ion_mass")

    if argument == 'alpha' or element_symbol(argument) == 'He' and Z == 2:
        return 6.644657230e-27*u.kg
    elif argument in ['p', 'p+'] or str(argument).lower() in \
            ['proton', 'protium']:  # Not H because conv atom weight has some D
        return const.m_p
    elif atomic_number(argument) == 1:
        if isinstance(argument, str) and '-1' in str(argument):
            return const.m_p
        elif argument == 1 and mass_numb == 1:
            return const.m_p

    try:
        isotope_symb = isotope_symbol(argument, mass_numb)
        if isotope_symb == 'D' and Z == 1:
            return 3.343583719e-27 * u.kg
        elif isotope_symb == 'T' and Z == 1:
            return 5.007356665e-27 * u.kg
        atomic_mass = isotope_mass(isotope_symb)
    except Exception:
        try:
            atomic_mass = standard_atomic_weight(argument)
        except Exception:
            errormessage = "No isotope mass or standard atomic weight is " +\
                "available to get ion mass for " + str(argument)
            if isinstance(mass_numb, int):
                errormessage += " with mass number " + str(mass_numb)

            raise ValueError(errormessage)

    m_i = (atomic_mass - Z*const.m_e).to(u.kg)

    return m_i


def nuclear_binding_energy(argument, mass_numb=None):
    """Returns the nuclear binding energy associated with an isotope.

    Parameters
    ----------
    argument: string or integer
        A string representing an element or isotope, or an integer representing
        the atomic number of an element.

    mass_numb: integer
        The mass number of an isotope, which is required if and only if the
        first argument can only be used

    Returns
    -------
    binding_energy: Quantity
        The binding energy of the nucleus.

    Raises
    ------
    ValueError
        If the isotope cannot be uniquely determined from the inputs.

    Examples
    --------
    >>> nuclear_binding_energy('Fe-56')
    <Quantity 7.674002331034521e-11 J>
    >>> nuclear_binding_energy(26, 56)
    <Quantity 7.674002331034521e-11 J>
    >>> nuclear_binding_energy('p')  # proton
    <Quantity 0.0 J>

    >>> from astropy import units as u
    >>> before = nuclear_binding_energy("D") + nuclear_binding_energy("T")
    >>> after = nuclear_binding_energy("alpha")
    >>> (after - before).to(u.MeV)  # released energy from D + T --> alpha + n
    <Quantity 17.589296207151556 MeV>

    """

    if argument == 'n' and mass_numb is None or mass_numb == 1:
        return 0.0 * u.J

    isotopic_symbol = isotope_symbol(argument, mass_numb)
    isotopic_mass = isotope_mass(isotopic_symbol)
    number_of_protons = atomic_number(argument)

    if mass_numb is None:
        mass_numb = mass_number(argument)
    number_of_neutrons = mass_numb - number_of_protons

    if number_of_protons == 1 and number_of_neutrons == 0:
        binding_energy = 0.0 * u.J
    else:
        mass_of_nucleons = (number_of_protons * const.m_p +
                            number_of_neutrons * const.m_n)
        mass_defect = mass_of_nucleons - isotopic_mass
        binding_energy = mass_defect * const.c**2

    return binding_energy.to(u.J)


def energy_from_nuclear_reaction(reaction):
    """Returns the released energy from a nuclear fusion reaction.

    Parameters
    ----------
    reaction: string
        A string representing the reaction, like "D + T --> alpha + n"

    Returns
    -------
    energy: the change in nuclear binding energy, which will be positive
    if the reaction releases and negative if the reaction is energetically
    unfavorable.

    Raises
    ------
    ValueError:
        If the input is not a valid reaction, there is insufficient
        information to determine an isotope, or if the number of nucleons
        is not conserved during the reaction.

    TypeError:
        If the input is not a string.

    See also
    --------
    nuclear_binding_energy : finds the binding energy of an isotope

    Examples
    --------
    >>> energy_from_nuclear_reaction("D + T --> alpha + n")
    <Quantity 17.589296207151556 MeV>
    >>> triple_alpha1 = 'alpha + He-4 --> Be-8'
    >>> triple_alpha2 = 'Be-8 + alpha --> carbon-12'
    >>> energy_triplealpha1 = energy_from_nuclear_reaction(triple_alpha1)
    >>> energy_triplealpha2 = energy_from_nuclear_reaction(triple_alpha2)
    >>> print(energy_triplealpha1, '\n', energy_triplealpha2)
    -0.09183948324626812 MeV
    7.366586766240317 MeV
    >>> energy_triplealpha1.to(u.J)
    <Quantity -1.471430677988809e-14 J>
    >>> energy_triplealpha2.cgs
    <Quantity -1.4714306779888094e-07 erg>
    >>> energy_from_nuclear_reaction('alpha + alpha --> 2alpha')
    <Quantity 0.0 MeV>

    """

    import re

    def _get_isotopes_list(side):
        """Parse a side of a reaction to get a list of the isotopes."""
        pre_list = re.split(' \+ ', side)
        isotopes_list = []
        for item in pre_list:
            item = item.strip()
            try:
                if item == 'n':
                    symbol = 'n'
                else:
                    symbol = isotope_symbol(item)
                isotopes_list.append(symbol)
            except Exception:
                try:
                    multiplier_string = ''
                    while item[0].isdigit():
                        multiplier_string += item[0]
                        item = item[1:]
                    symbol = isotope_symbol(item)
                    for i in range(0, int(multiplier_string)):
                        isotopes_list.append(symbol)
                except Exception:
                    raise
        return isotopes_list

    def _mass_number_of_list(isotopes_list):
        """Find the total number of nucleons in a list of isotopes."""
        mass_numb = 0
        for symbol in isotopes_list:
            mass_numb += mass_number(symbol)
        return mass_numb

    def _add_binding_energies(isotopes_list):
        """Finds the total binding energy from a list of isotopes."""
        total_binding_energy = 0.0*u.MeV
        for isotope in isotopes_list:
            total_binding_energy += nuclear_binding_energy(isotope)
        return total_binding_energy

    if not isinstance(reaction, str):
        raise TypeError("The input of energy_from_nuclear_reaction "
                        "must be a string representing the reaction (e.g., "
                        "'D + T --> He + n')")

    try:
        LHS, RHS = re.split('-+>', reaction)
    except Exception:
        raise ValueError("The left and right hand sides of the reaction "
                         "should be separated by '-->'")

    reactants = _get_isotopes_list(LHS)
    products = _get_isotopes_list(RHS)

    mass_num_reactants = _mass_number_of_list(reactants)
    mass_num_products = _mass_number_of_list(products)

    if mass_num_reactants != mass_num_products:
        raise ValueError("Mass numbers on LHS and RHS do not match for "
                         "reaction "+reaction)

    binding_energy_before = _add_binding_energies(reactants)
    binding_energy_after = _add_binding_energies(products)
    energy = binding_energy_after - binding_energy_before

    return energy.to(u.MeV)


def known_isotopes(argument):
    """Return a list of all known isotopes of an element in order of mass
    number.

    Parameters
    ----------
    argument: integer or string
        A string or integer representing an atomic number or element, or a
        string represnting an isotope.

    Returns
    -------
    sorted_isotopes: list of strings or empty list
        List of all known isotopes of an element, sorted from lowest mass
        number to highest mass number.  This list includes all isotopes
        that have been discovered.

    See also
    --------
    common_isotopes : returns a list of all isotopes of an element with
        isotopic compositions greater than zero, sorted from highest
        abundance to lowest abundance.

    stable_isotopes : returns a list of all stable isotopes of an element,
        sorted from lowest to highest mass number, or an empty list if an
        element has no stable isotopes.

    Examples
    --------
    >>> known_isotopes('H')
    ['H-1', 'D', 'T', 'H-4', 'H-5', 'H-6', 'H-7']
    >>> known_isotopes('helium')
    ['He-3', 'He-4', 'He-5', 'He-6', 'He-7', 'He-8', 'He-9', 'He-10']

    """

    element = element_symbol(argument)

    isotopes_of_element = []
    for isotope in Isotopes.keys():
        if element + '-' in isotope and isotope[0:len(element)] == element:
            isotopes_of_element.append(isotope)

    if element == 'H':
        isotopes_of_element.insert(1, 'D')
        isotopes_of_element.insert(2, 'T')

    mass_numbers = [mass_number(isotope) for isotope in isotopes_of_element]

    sorted_isotopes = [mass_number for (isotope, mass_number) in
                       sorted(zip(mass_numbers, isotopes_of_element))]

    return sorted_isotopes


def common_isotopes(argument):
    """Returns a list of isotopes of an element with isotopic compositions
    greater than zero, or an empty list if no isotopes have nonzero
    isotopic compositions.

    Parameters
    ----------
    argument: integer or string
        A string or integer representing an atomic number or element, or a
        string represnting an isotope.

    Returns
    -------
    sorted_isotopes: list of strings or empty list
        List of all isotopes of an element with isotopic compositions greater
        than zero, sorted from most abundant to least abundant.  If no isotopes
        have isotopic compositions greater than zero, this function returns
        an empty list.

    See also
    --------
    known_isotopes : returns a list of all isotopes of an element that have
        been discovered, sorted from lowest to highest mass number.

    stable_isotopes : returns a list of all stable isotopes of an element,
        sorted from lowest to highest mass number, or an empty list if an
        element has no stable isotopes.

    Examples
    --------
    >>> common_isotopes('H')
    ['H-1', 'D']
    >>> common_isotopes(44)
    ['Ru-102', 'Ru-104', 'Ru-101', 'Ru-99', 'Ru-100', 'Ru-96', 'Ru-98']
    >>> common_isotopes('beryllium')
    ['Be-9']
    >>> common_isotopes('Pb-209')
    ['Pb-208', 'Pb-206', 'Pb-207', 'Pb-204']
    >>> common_isotopes(118)
    []

    """

    isotopes = known_isotopes(argument)

    CommonIsotopes = [isotope for isotope in isotopes if
                      'isotopic_composition' in Isotopes[isotope].keys()]

    isotopic_compositions = [Isotopes[isotope]['isotopic_composition']
                             for isotope in CommonIsotopes]

    sorted_isotopes = [iso_comp for (isotope, iso_comp) in
                       sorted(zip(isotopic_compositions, CommonIsotopes))]

    sorted_isotopes.reverse()

    return sorted_isotopes


def stable_isotopes(argument):
    """Returns a list of all stable isotopes of an element, or an empty list
    if an element has no stable isotopes.

    Parameters
    ----------
    argument: integer or string
        A string or integer representing an atomic number or element, or a
        string represnting an isotope.

    Returns
    -------
    StableIsotopes: list of strings or empty list
        List of all stable isotopes of an element, sorted from lowest mass
        number.  If an element has no stable isotopes, this function returns
        an empty list.

    See also
    --------
    known_isotopes : returns a list of all isotopes of an element that have
        been discovered, sorted from lowest to highest mass number.

    common_isotopes : returns a list of all isotopes of an element with
        isotopic compositions greater than zero, sorted from highest
        abundance to lowest abundance.

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

    """

    KnownIsotopes = known_isotopes(argument)

    StableIsotopes = [isotope for isotope in KnownIsotopes if
                      Isotopes[isotope]['is_stable']]

    return StableIsotopes


def isotopic_composition(argument, mass_numb=None):
    """Returns the isotopic composition if known, and otherwise zero.

    Parameters
    ----------
    argument: string or integer
        A string representing an element or isotope, or an integer representing
        the atomic number of an element.

    mass_numb: integer
        The mass number of an isotope, which is required if and only if the
        first argument can only be used

    Returns
    -------
    iso_comp: float
        The isotopic composition

    Raises
    ------
    ValueError
        Invalid isotope input, or the input corresponded to neutrons.

    Examples
    --------
    >>> isotopic_composition('Pb-208')
    0.524
    >>> isotopic_composition('hydrogen', 1)
    0.999885
    >>> isotopic_composition(118, 294)  # Og-294
    0.0

    """

    try:
        isotope = isotope_symbol(argument, mass_numb)
    except Exception:
        raise ValueError("Invalid isotope in isotopic_composition.")

    if isotope == 'n':
        raise ValueError("Neutrons do not have an isotopic composition.")

    try:
        iso_comp = Isotopes[isotope]['isotopic_composition']
    except Exception:
        iso_comp = 0.0

    return iso_comp


def charge_state(argument):
    """Returns the charge state of an ion.

    Parameters
    ----------
    argument : string
        String representing an element or isotope followed by charge state
        information.

    Returns
    -------
    Z : integer
        The charge state, or None if it is not available.

    Raises
    ------
    ValueError:
        If the charge state or isotope information is invalid, or the charge
        state exceeds the atomic number.

    Notes
    -----
    This function supports two formats for the charge state information.

    The first format is a string that has information for the element
    or isotope at the beginning, a space in between, and the charge
    state information in the form of an integer followed by a plus or
    minus sign, or a plus or minus sign followed by an integer.

    The second format is a string containing element information at
    the beginning, following by one or more plus or minus signs.

    Examples
    --------
    >>> charge_state('Fe-56 2+')
    2
    >>> charge_state('He -2)
    -2
    >>> charge_state('H+')
    1
    >>> charge_state('N-14++')
    2

    """

    argument, Z = __extract_charge_state(argument)

    try:
        atomic_numb = atomic_number(argument)
    except Exception:
        raise ValueError("Invalid element or isotope information in "
                         "charge_state")

    if Z is not None and Z > atomic_numb:
        raise ValueError("The charge state cannot be greater than the atomic "
                         "number.")

    if Z is not None and Z < -atomic_numb - 1:
        raise UserWarning("Element " + element_symbol(argument) + " has a "
                          "charge of " + str(Z) + " which is "
                          "unlikely to occur in nature.")

    return Z


def __extract_charge_state(argument):
    """Splits strings containing element or isotope and charge state
    information into a string without the charge state information and
    the charge state as an integer (or None if no charge state
    information is given).

    Parameters
    ----------
    argument : string
        String containing information for an element or isotope in any of
        the allowed formats, followed by charge state information.

    Returns
    -------
    argument : string
        The original string with charge state information removed.

    Z : integer
        The charge state of an ion (e.g., this will return 1 if one electron
        has been removed and -1 if one electron has been gained)

    Notes
    -----
    If the argument is not a string, this function will return the
    original argument and None.

    This function does not check that the charge state is valid.

    Examples
    --------
    >>> isotope, Z = __extract_charge_state('Fe-56+++')
    >>> print(isotope)
    'Fe-56'
    >>> print(Z)
    3
    >>> __extract_charge_state('D +1')
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
            raise ValueError("The following input does not have valid charge "
                             "state information: " + argument + ion_info)

    elif argument.endswith(('-', '+')):  # For cases like 'Fe++' or 'Si-'

        char = argument[-1]
        match = re.match(r"["+char+"]*", argument[::-1])
        charge_state = match.span()[1]

        if char == '-':
            charge_state = -charge_state

        if argument.count(char) != np.abs(charge_state):
            raise ValueError("The following input does not have valid charge "
                             "state information: " + argument)

        argument = argument[0:len(argument)-match.span()[1]]

    else:
        charge_state = None

    if charge_state is not None and charge_state < -6:
        raise UserWarning("Element " + element_symbol(argument) + " has a "
                          "charge of " + str(charge_state) + " which is "
                          "unlikely to occur in nature.")

    return argument, charge_state
