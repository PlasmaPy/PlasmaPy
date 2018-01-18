import numpy as np
import re
import warnings
from typing import (Union, Optional, Any, Tuple)

from .elements import (_atomic_symbols, _atomic_symbols_dict, _Elements)

from .isotopes import _Isotopes

from .particles import (_is_special_particle,
                        _get_standard_symbol,
                        _case_sensitive_aliases,
                        _case_insensitive_aliases)

from ..utils import (AtomicWarning,
                     InvalidElementError,
                     InvalidIsotopeError,
                     InvalidIonError,
                     AtomicError,
                     InvalidParticleError,
                     ChargeError)


def _parse_and_check_atomic_input(
        argument: Union[str, int],
        mass_numb: int = None,
        Z: int = None):
    r"""Parses information about a particle into a dictionary
    containing standard symbols, while checking to make sure
    that the particle is valid.
    Parameters
    ----------
    argument : string or integer
        String containing information for an element, isotope, or ion
        in any of the allowed formats; or an integer representing an
        atomic number.
    mass_numb : integer, optional
        The mass number of an isotope.
    Z : integer, optional
        The integer charge of an ion.
    Returns
    -------
    nomenclature_dict : dict
        A dictionary containing information about the element, isotope,
        or ion.  The key 'symbol' corresponds to the particle symbol
        containing the most information, 'element' corresponds to the
        atomic symbol, 'isotope' corresponds to the isotope symbol,
        'ion' corresponds to the ion symbol, 'mass_numb' corresponds
        to the mass number, and 'Z' corresponds to the integer charge.
        The corresponding items will be given by None if the necessary
        information is not provided.
    Raises
    ------
    InvalidParticleError
        If the arguments do not correspond to a valid particle or
        antiparticle.
    InvalidElementError
        If the particle is valid but does not correspond to an element,
        ion, or isotope.
    """

    def _atomic_number_to_symbol(atomic_numb: int):
        r"""Returns the atomic symbol associated with an integer
        representing an atomic number, or raises an InvalidParticleError
        if the atomic number does not represent a known element."""

        if atomic_numb in _atomic_symbols.keys():
            element = _atomic_symbols[atomic_numb]
            return element
        else:
            raise InvalidParticleError

    def _extract_charge(arg: str):
        r"""Receives a string representing an element, isotope, or ion.
        Returns a tuple containing a string that should represent an
        element or isotope, and either an integer representing the
        charge or None if no charge information is provided.  Raises
        an InvalidParticleError if charge information is inputted
        incorrectly."""

        if arg.count(' ') == 1:  # Cases like 'H 1-' and 'Fe-56 1+'
            isotope_info, charge_info = arg.split(' ')

            sign_indicator_only_on_one_end = (
                    charge_info.endswith(('-', '+')) ^
                    charge_info.startswith(('-', '+')))

            just_one_sign_indicator = (
                    (charge_info.count('-') == 1 and
                     charge_info.count('+') == 0) or
                    (charge_info.count('-') == 0 and
                     charge_info.count('+') == 1))

            if not sign_indicator_only_on_one_end and just_one_sign_indicator:
                raise InvalidParticleError

            if '-' in charge_info:
                sign = -1
            elif '+' in charge_info:
                sign = 1

            charge_str = charge_info.strip('+-')

            try:
                Z_from_arg = sign * int(charge_str)
            except ValueError:
                raise InvalidParticleError

        elif arg.endswith(('-', '+')):  # Cases like 'H-' and 'Pb-209+++'
            char = arg[-1]
            match = re.match(f"[{char}]*", arg[::-1])
            Z_from_arg = match.span()[1]
            isotope_info = arg[0:len(arg) - match.span()[1]]

            if char == '-':
                Z_from_arg = -Z_from_arg
            if isotope_info.endswith(('-', '+')):
                raise InvalidParticleError
        else:
            isotope_info = arg
            Z_from_arg = None

        return isotope_info, Z_from_arg

    def _extract_mass_number(isotope_info: str):
        r"""Receives a string representing an element or isotope.
        Returns a tuple containing a string that should represent
        an element, and either an integer representing the mass
        number or None if no mass number is available.  Raises an
        InvalidParticleError if the mass number information is
        inputted incorrectly."""

        if isotope_info == 'D':
            element_info, mass_numb = 'H', 2
        elif isotope_info == 'T':
            element_info = 'H'
            mass_numb = 3
        elif '-' not in isotope_info:
            element_info = isotope_info
            mass_numb = None
        elif isotope_info.count('-') == 1:
            element_info, mass_numb_str = isotope_info.split('-')
            try:
                mass_numb = int(mass_numb_str)
            except ValueError:
                raise InvalidParticleError
        else:
            element_info = isotope_info
            mass_numb = None

        return element_info, mass_numb

    def _get_element(element_info: str) -> str:
        r"""Receives a string representing an element's symbol or
        name, and returns a string representing the atomic symbol."""

        if element_info.lower() in _atomic_symbols_dict.keys():
            element = _atomic_symbols_dict[element_info.lower()]
        elif element_info in _atomic_symbols.values():
            element = element_info
        else:
            raise InvalidParticleError

        return element

    def _reconstruct_isotope_symbol(element: str, mass_numb: int) -> str:
        r"""Receives a string representing an atomic symbol and an
        integer representing a mass number.  Returns the isotope symbol
        or None if no mass number information is available."""

        if mass_numb is not None:
            isotope = f"{element}-{mass_numb}"

            if isotope == 'H-2':
                isotope = 'D'
            elif isotope == 'H-3':
                isotope = 'T'

            if isotope not in _Isotopes.keys():
                print(isotope)
                raise InvalidParticleError
        else:
            isotope = None

        return isotope

    def _reconstruct_ion_symbol(
            element: str, isotope: int = None, Z: int = None):
        r"""Receives a string representing an atomic symbol and/or a
        string representing an isotope, and an integer representing the
        integer charge.  Returns a string representing the ion symbol,
        or None if no charge information is available."""

        if Z is not None:
            if Z < 0:
                sign = '-'
            else:
                sign = '+'

            if isotope is None:
                base = element
            else:
                base = isotope

            ion = f"{base} {np.abs(Z)}{sign}"
        else:
            ion = None

        if ion == 'H-1 1+':
            ion = 'p'

        return ion

    argument = _get_standard_symbol(argument)  # Deals with aliases

    if _is_special_particle(argument):
        if (mass_numb is not None) or (Z is not None):
            raise InvalidParticleError
        else:
            raise InvalidElementError

    if isinstance(argument, int):
        element = _atomic_number_to_symbol(argument)
        Z_from_arg = None
        mass_numb_from_arg = None
    elif isinstance(argument, str):
        isotope_info, Z_from_arg = _extract_charge(argument)
        element_info, mass_numb_from_arg = _extract_mass_number(isotope_info)
        element = _get_element(element_info)
    else:
        raise InvalidParticleError

    # Check the validity and consistency of mass numbers

    if mass_numb is not None and mass_numb_from_arg is not None:
        if mass_numb != mass_numb_from_arg:
            raise InvalidParticleError
        else:
            warnings.warn("Redundant...", AtomicWarning)

    if mass_numb_from_arg is not None:
        mass_numb = mass_numb_from_arg

    # Check the validity and consistency of integer charge

    if Z is not None and Z_from_arg is not None:
        if Z != Z_from_arg:
            raise InvalidParticleError
        else:
            warnings.warn("Redundant...", AtomicWarning)

    if Z_from_arg is not None:
        Z = Z_from_arg

    if isinstance(Z, int):
        if Z > _Elements[element]['atomic_number']:
            raise InvalidParticleError
        elif Z < -3:
            warnings.warn("Unrealistically negative charge state...")

    isotope = _reconstruct_isotope_symbol(element, mass_numb)
    ion = _reconstruct_ion_symbol(element, isotope, Z)

    if ion:
        symbol = ion
    elif isotope:
        symbol = isotope
    else:
        symbol = element

    nomenclature_dict = {
        'symbol': symbol,
        'element': element,
        'isotope': isotope,
        'ion': ion,
        'mass_numb': mass_numb,
        'Z': Z,
    }

    return nomenclature_dict
