r"""Functions that deal with string representations of atomic symbols
and numbers."""

import re
import warnings
from typing import (Optional, Any, Tuple)
from .parsing import _dealias_particle_aliases
from .particle_class import Particle
from .particle_input import particle_input
from ..utils import (AtomicWarning, InvalidParticleError)

# The @particle_input decorator takes the inputs for a function or
# method and passes through the corresponding instance of the Particle
# class via the argument annotated with Particle.  If the argument is
# named element, isotope, or ion; then this decorator will raise an
# InvalidElementError, InvalidIsotopeError, or InvalidIonError,
# respectively.  The Particle class constructor will raise an
# InvalidParticleError if the input does not correspond to a valid
# particle.


@particle_input
def atomic_symbol(element: Particle) -> str:
    r"""Returns the atomic symbol.

    Parameters
    ----------

    element: string, integer, or Particle
        A string representing an element, isotope, or ion; or an
        integer or string representing an atomic number.

    Returns
    -------

    symbol: string
        The atomic symbol of the element, isotope, or ion.

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

    element_name : returns the name of an element.

    isotope_symbol : returns the isotope symbol instead of the atomic symbol.

    ion_symbol : returns the ion symbol instead of the atomic symbol.

    particle_symbol : returns the symbol of any valid particle.

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
    return element.element


@particle_input
def isotope_symbol(isotope: Particle, mass_numb: int = None) -> str:
    r"""Returns the symbol representing an isotope.

    Parameters
    ----------

    isotope: string, integer, or Particle
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    mass_numb: integer or string, optional
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

    atomic_symbol : returns the atomic symbol instead of the isotope symbol

    ion_symbol : returns the ion symbol instead of the isotope symbol.

    particle_symbol : returns the symbol of any valid particle.


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
    return isotope.isotope


@particle_input
def ion_symbol(ion: Particle, mass_numb: int = None, Z: int = None) -> str:
    r"""Returns the symbol representing an ion.

    Parameters
    ----------

    ion: int, str, or Particle
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    mass_numb: integer or string
        An integer or string representing the mass number of the ion.

    Z: integer or string
        An integer or string representing the integer charge of the ion.

    Returns
    -------

    symbol: string
        The ion symbol. The result will generally be returned as
        something like 'He-4 2+', 'D 1+', or 'p+'.

    Raises
    ------

    InvalidIonError
        If the arguments correspond to a valid particle but not a valid ion.

    InvalidParticleError
        If arguments do not correspond to a valid particle or contradictory
        information is provided.

    TypeError
        If ion is not a string, integer, or Particle; or if either of
        mass_numb or Z is not an integer or a string representing an integer.

    Warns
    -----

    AtomicWarning
        If redundant mass number or charge information is provided.

    See also
    --------

    atomic_symbol : returns the atomic symbol instead of the ion symbol

    isotope_symbol : returns the isotope symbol instead of the ion symbol

    particle_symbol : returns the symbol of any valid particle

    Notes
    -----

    This function returns the symbol of the element rather than the
    symbol of an isotope.  For example, 'deuterium', 'T', or
    'hydrogen-2' will yield 'H'; 'alpha' will yield 'He'; and
    'iron-56' or 'Fe-56' will yield 'Fe'.

    Examples
    --------
    >>> ion_symbol('alpha')
    'He-4 2+'
    >>> ion_symbol(79, mass_numb=197, Z=12)
    'Au-197 12+'
    >>> ion_symbol('proton')
    'p+'
    >>> ion_symbol('D', Z=1)
    'D 1+'
    """
    return ion.ion


@particle_input
def particle_symbol(particle: Particle, mass_numb: int = None,
                    Z: int = None) -> str:
    r"""Returns the symbol of a particle.

    Parameters
    ----------

    particle: int, str, or Particle
        A string representing a particle, element, isotope, or ion or an
        integer representing an atomic number

    mass_numb: integer or string
        An integer or string representing the mass number of an isotope.

    Z: integer or string
        An integer or string representing the integer charge of an ion.

    Returns
    -------

    symbol: string
        The particle symbol, containing charge and mass number information
        when available. The result will generally be returned as
        something like 'e-', 'Fe', 'He-4 2+', 'D', 'n', 'mu-', or 'p+'.

    Raises
    ------

    InvalidParticleError
        If arguments do not correspond to a valid particle or contradictory
        information is provided.

    TypeError
        If ion is not a string, integer, or Particle; or if either of
        mass_numb or Z is not an integer or a string representing an integer.

    Warns
    -----

    AtomicWarning
        If redundant mass number or charge information is provided.

    See also
    --------

    atomic_symbol : returns the atomic symbol instead

    isotope_symbol : returns the isotope symbol instead

    ion_symbol : returns the ion symbol instead

    Examples
    --------
    >>> particle_symbol('electron')
    'e-'
    >>> particle_symbol('proton')
    'p+'
    >>> particle_symbol('alpha')
    'He-4 2+'
    >>> particle_symbol('H-1', Z=-1)
    'H-1 1-'
    """
    return particle.particle


@particle_input
def element_name(element: Particle) -> str:
    r"""Returns the name of an element.

    Parameters
    ----------

    argument : string, integer, or Particle
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

    isotope_symbol : returns the isotope symbol

    particle_symbol : returns the symbol of any valid particle

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
    return element.element_name


def _is_neutron(argument: Any, mass_numb: int = None) -> bool:
    r"""Returns True if the argument corresponds to a neutron, and
    False otherwise."""

    warnings.warn("_is_neutron is deprecated.", DeprecationWarning)

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

    warnings.warn("_is_hydrogen is deprecated.", DeprecationWarning)

    case_sensitive_aliases = ['p', 'p+', 'H', 'D', 'T']

    case_insensitive_aliases = ['proton', 'protium', 'deuterium',
                                'deuteron', 'triton', 'tritium', 'hydrogen']

    for mass_numb in range(1, 8):
        case_sensitive_aliases.append(f'H-{mass_numb}')
        case_insensitive_aliases.append(f'hydrogen-{mass_numb}')

    if isinstance(argument, str):

        argument, Z = _extract_integer_charge(argument)

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

    warnings.warn("_is_electron is deprecated.", DeprecationWarning)

    if not isinstance(arg, str):
        return False

    return arg in ['e', 'e-'] or arg.lower() == 'electron'


def _is_positron(arg: Any) -> bool:
    r"""Returns True if the argument corresponds to a positron, and False
    otherwise."""

    warnings.warn("_is_positron is deprecated.", DeprecationWarning)

    if not isinstance(arg, str):
        return False

    return arg == 'e+' or arg.lower() == 'positron'


def _is_antiproton(arg: Any) -> bool:
    r"""Returns True if the argument corresponds to an antiproton, and
    False otherwise."""

    warnings.warn("_is_antiproton is deprecated.", DeprecationWarning)

    if not isinstance(arg, str):
        return False

    return arg == 'p-' or arg.lower() == 'antiproton'


def _is_antineutron(arg: Any) -> bool:
    r"""Returns True if the argument corresponds to an antineutron, and
    False otherwise."""

    warnings.warn("_is_antineutron is deprecated.", DeprecationWarning)

    if not isinstance(arg, str):
        return False

    return arg.lower() == 'antineutron'


def _is_proton(arg: Any, Z: int = None, mass_numb: int = None) -> bool:
    r"""Returns True if the argument corresponds to a proton, and
    False otherwise.  This function returns False for 'H-1' if no
    charge state is given."""

    warnings.warn("_is_proton is deprecated.", DeprecationWarning)

    argument, Z_from_arg = _extract_integer_charge(arg)

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

    warnings.warn("_is_alpha is deprecated.", DeprecationWarning)

    if not isinstance(arg, str):
        return False

    if arg.lower() == 'alpha':
        return True
    elif 'alpha' in arg.lower():
        raise InvalidParticleError(
            f"{arg} is an invalid representation of an alpha particle")
    else:
        arg, Z = _extract_integer_charge(arg)

        if Z != 2 or arg[-2:] != '-4':
            return False
        else:

            dash_position = arg.find('-')
            arg = arg[:dash_position]

            return arg.lower() == 'helium' or arg == 'He'
