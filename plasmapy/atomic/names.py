r"""Functions that deal with string representations of atomic symbols
and numbers."""

import re
import warnings
from typing import (Union, Optional, Any, Tuple)

from .elements import (_atomic_symbols, _atomic_symbols_dict, _Elements)

from .isotopes import _Isotopes

from .parsing import (_is_special_particle, _get_standard_symbol)

from ..utils import (AtomicWarning,
                     InvalidElementError,
                     InvalidIsotopeError,
                     InvalidIonError,
                     AtomicError,
                     InvalidParticleError,
                     ChargeError)

# TODO: Create an ion_symbol function
# TODO: Create a particle_symbol function


def atomic_symbol(argument: Union[str, int]) -> str:
    r"""Returns the atomic symbol.

    Parameters
    ----------

    argument: string or integer
        A string representing an element, isotope, or ion; or an
        integer representing an atomic number.

    Returns
    -------

    symbol: string
        The atomic symbol of the element, isotope, or nucleon.

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

    if _is_special_particle(argument):
        raise InvalidElementError(f"{argument} is not a valid element.")

    try:
        argument, Z = _extract_charge_state(argument)
    except InvalidParticleError:
        raise InvalidParticleError("Invalid charge in atomic_symbol")

    if not isinstance(argument, (str, int)):
        raise TypeError("The first argument in atomic_symbol must be either "
                        "a string representing an element or isotope, or an "
                        "integer representing the atomic number (or 0 for "
                        "neutrons).")

    if isinstance(argument, str) and argument.isdigit():
        argument = int(argument)

    if isinstance(argument, int):

        try:
            element = _atomic_symbols[argument]
        except KeyError:
            raise InvalidParticleError(f"{argument} is an invalid atomic "
                                       "number in atomic_symbol.")

    elif _is_hydrogen(argument):
        element = 'H'
    elif _is_alpha(argument):
        element = 'He'
    elif isinstance(argument, str):

        if argument.count('-') == 1:
            dash_position = argument.find('-')
            mass_numb = argument[dash_position+1:]
            if not mass_numb.isdigit():
                raise InvalidParticleError("Invalid isotope format in "
                                           "atomic_symbol")
            argument = argument[:dash_position]
        else:
            mass_numb = ''

        if argument.lower() in _atomic_symbols_dict.keys():
            element = _atomic_symbols_dict[argument.lower()]
        elif argument in _atomic_symbols.values():
            element = argument.capitalize()
        else:
            raise InvalidParticleError(f"{argument} is an invalid argument "
                                       "for atomic_symbol")

        if mass_numb.isdigit():

            isotope = element.capitalize() + '-' + mass_numb

            if isotope not in _Isotopes.keys():
                raise InvalidParticleError(
                    "The input in atomic_symbol corresponding "
                    f"to {isotope} is not a valid isotope.")

    if Z is not None and \
            Z > _Elements[element]['atomic_number']:
        raise InvalidParticleError("Cannot have an ionization state greater "
                                   "than the atomic number.")

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

    InvalidIsotopeError
        If the argument is a valid particle but not a valid isotope.

    InvalidParticleError
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

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

    # TODO: Remove this functionality when particle_symbol comes online
    if _is_neutron(argument, mass_numb):
        return 'n'

    if _is_special_particle(argument):
        raise InvalidIsotopeError("The argument {argument} does not "
                                  "correspond to a valid isotope in "
                                  "isotope_symbol.")

    try:
        element = atomic_symbol(argument)
    except InvalidParticleError:
        raise InvalidParticleError(
            f"The argument {argument} to isotope_symbol is not a valid "
            f"particle.")
    except InvalidElementError:
        raise InvalidIsotopeError(f"The argument {argument} to isotope_symbol "
                                  f"does not correspond to a valid element or "
                                  f"isotope.")

    # If the argument is already in our standard form for an isotope,
    # return the argument.

    if mass_numb is None and argument in _Isotopes.keys():
        return argument

    if isinstance(argument, str):
        argument, Z = _extract_charge_state(argument)

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
            raise InvalidParticleError(
                "Invalid atomic number in isotope_symbol")
        if mass_numb is None:
            raise InvalidIsotopeError("Insufficient information to determine "
                                      "isotope in isotope_symbol.")

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
            raise InvalidIsotopeError(
                "Insufficient information to determine the mass"
                " number from the inputs to isotope_symbol.")

        if mass_numb is not None and mass_numb_from_arg is not None:
            if mass_numb == mass_numb_from_arg:
                warnings.warn(
                    "Redundant mass number information in isotope_symbol "
                    f"from inputs: {argument}, {mass_numb}", AtomicWarning)
            else:  # coveralls: ignore
                raise InvalidParticleError(
                    "Contradictory mass number information in isotope_symbol.")

        if mass_numb_from_arg is not None:
            mass_numb = mass_numb_from_arg

    isotope = f"{element}-{mass_numb}"

    if isotope == 'H-2':
        isotope = 'D'
    elif isotope == 'H-3':
        isotope = 'T'

    if _Elements[element]['atomic_number'] > mass_numb:
        raise InvalidParticleError("The atomic number cannot exceed the mass "
                                   "number in isotope_symbol.")

    return isotope


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
        name = _Elements[element]["name"]
    except InvalidElementError:
        raise InvalidElementError("Invalid element in element_name.")
    except InvalidParticleError:
        raise InvalidParticleError("Invalid particle in element_name.")
    except TypeError:
        raise TypeError("Invalid input to element_name.")

    return name


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

    InvalidParticleError
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

    argument = _get_standard_symbol(argument)

    if argument in ['n', 'antineutron'] or 'nu_' in argument:
        return argument, 0
    if argument in ['e-', 'mu-', 'tau-', 'p-']:
        return argument, -1
    elif argument in ['e+', 'mu+', 'tau+', 'p']:
        return argument, 1

    if argument.lower() == 'alpha':
        return 'He-4', 1
    elif 'alpha' in argument.lower():
        raise InvalidParticleError("Invalid representation of an alpha"
                                   " particle")

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
            charge_state = sign * int(charge)
            check3 = True
        except Exception:
            check3 = False

        if not (check1 and check2 and check3):
            raise InvalidParticleError(
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
            raise InvalidParticleError("Invalid charge state information")

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
            raise InvalidParticleError("Invalid charge state of hydrogen")

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

    argument, Z_from_arg = _extract_charge_state(arg)

    if (Z is None) == (Z_from_arg is None):
        return False
    else:
        if Z is None:
            Z = Z_from_arg

    try:
        isotope = isotope_symbol(arg, mass_numb)
    except Exception:
        return False

    return isotope == 'H-1' and Z == 1


def _is_alpha(arg: Any) -> bool:
    r"""Returns True if the argument corresponds to an alpha particle,
    and False otherwise."""

    if not isinstance(arg, str):
        return False

    if arg.lower() == 'alpha':
        return True
    elif 'alpha' in arg.lower():
        raise InvalidParticleError(
            f"{arg} is an invalid representation of an alpha particle")
    else:
        arg, Z = _extract_charge_state(arg)

        if Z != 2 or arg[-2:] != '-4':
            return False
        else:

            dash_position = arg.find('-')
            arg = arg[:dash_position]

            return arg.lower() == 'helium' or arg == 'He'
