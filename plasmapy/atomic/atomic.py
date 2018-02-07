r"""Functions that retrieve or are related to elemental or isotopic data."""

import warnings
from typing import (Union, Optional, List, Any)

from astropy import units as u, constants as const
from astropy.units import Quantity

from .elements import _Elements
from .isotopes import _Isotopes
from .particle_class import Particle
from .particle_input import particle_input

from ..utils import (
    MissingAtomicDataError,
    InvalidParticleError,
    InvalidElementError,
    InvalidIsotopeError,
    InvalidIonError,
)

from .symbols import atomic_symbol


@particle_input
def atomic_number(element: Particle) -> str:
    r"""Returns the number of protons in an atom, isotope, or ion.

    Parameters
    ----------

    element: string or Particle
        A string representing an element, isotope, or ion; or an
        instance of the Particle class.

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
    return element.atomic_number


@particle_input
def mass_number(isotope: Particle) -> int:
    r"""Get the mass number (the number of protons and neutrons) of an
    isotope.

    Parameters
    ----------

    isotope : string or Particle
        A string representing an isotope or a neutron; or an instance of
        the Particle class.

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
    >>> mass_number("alpha")
    4

    """
    return isotope.mass_number


@particle_input(must_be='element', exclude={'isotope', 'ion'})
def standard_atomic_weight(element: Particle) -> Quantity:
    r"""Returns the standard (conventional) atomic weight of an element
    based on the relative abundances of isotopes in terrestrial
    environments.

    Parameters
    ----------

    element: string, integer, or Particle
        A string representing an element or an integer representing an
        atomic number, or an instance of the Particle class

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
    return element.standard_atomic_weight


@particle_input(must_be='isotope', exclude='ion')
def isotope_mass(isotope: Particle,
                 mass_numb: int = None) -> Quantity:
    r"""Return the mass of an isotope.

    Parameters
    ----------

    isotope : string, integer, or Particle
        A string representing an element, isotope, or ion; an integer
        representing an atomic number; or an instance of the Particle
        class.

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
    >>> isotope_mass("He", mass_numb=4)
    <Quantity 4.00260325 u>

    """
    return isotope.mass


@particle_input
def ion_mass(particle: Particle,
             Z: int = None,
             mass_numb: int = None) -> Quantity:
    r"""Returns the mass of an ion by finding the standard atomic
    weight of an element or the atomic mass of an isotope, and then
    accounting for the change in mass due to loss of electrons from
    ionization.

    Parameters
    ----------

    particle: string, integer, or Particle
        A string representing an element, isotope, or ion; an integer
        representing an atomic number; or an instance of the Particle
        class.

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

    InvalidParticleError
        If the arguments do not correspond to a valid particle.

    InvalidIonError
        If the argument represents a particle other than an ion, the
        ionization state exceeds the atomic number, or no isotope mass
        or standard atomic weight is available.

    MissingAtomicDataError
        If the standard atomic weight or isotopic mass is not known.

    See also
    --------

    standard_atomic_weight : returns the conventional atomic mass of an
        element based on terrestrial values and assuming the atom is
        neutral.

    isotope_mass : returns the mass of an isotope (if available)
        assuming the atom is neutral.

    Notes
    -----

    This function in general finds the mass of an isotope (or the
    standard atomic weight based on terrestrial values if a unique
    isotope cannot be identified), and then subtracts the mass of Z
    electrons.

    If Z is not provided as an input, then ion_mass assumes that the ion
    is singly ionized while issuing a DeprecationWarning.

    Specific values are returned for some special particles such as
    protons, deuterons, tritons, and positrons.

    Calling ion_mass('H') does not return the mass of a proton but
    instead uses hydrogen's standard atomic weight based on
    terrestrial values.  To get the mass of a proton, use
    ion_mass('p').

    Examples
    --------

    >>> print(ion_mass('p').si.value)
    1.672621898e-27
    >>> ion_mass('H+')  # assumes terrestrial abundance of D
    <Quantity 1.67291241e-27 kg>
    >>> ion_mass('H+') == ion_mass('p')
    False
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
    >>> ion_mass('e+')
    <<class 'astropy.constants.codata2014.CODATA2014'> name='Electron mass' value=9.10938356e-31 uncertainty=1.1e-38 unit='kg' reference='CODATA 2014'>
    """

    # TODO: Remove deprecated functionality elsewhere in the code

    if particle.ion or particle.particle in {'e+'}:
        return particle.mass
    elif particle.particle == 'n':
        raise InvalidIonError
    elif particle.isotope in ['H-1']:
        warnings.warn("Use 'p+' instead of 'H-1' to refer to protons",
                      DeprecationWarning)
        return const.m_p
    elif particle.element:
        warnings.warn("The assumption that particles are singly ionized in"
                      "ion_mass is deprecated.", DeprecationWarning)
        return particle.mass - const.m_e
    elif particle.mass is not None:
        warnings.warn('The use of ion_mass for particles besides ions and '
                      'positrons is deprecated.', DeprecationWarning)
        return particle.mass
    else:
        raise InvalidIonError


@particle_input(exclude={'neutrino', 'antineutrino'})
def particle_mass(particle: Particle, *, Z: int = None,
                  mass_numb: int = None) -> Quantity:
    r"""Returns the mass of a particle.

    Parameters
    ----------

    particle: string, integer, or Particle
        A string representing an element, isotope, ion, or special
        particle; an integer representing an atomic number; or an
        instance of the Particle class.

    Z: integer (optional)
        The ionization state of the ion (defaulting to a charge of
        Z=1)

    mass_numb: integer (optional)
        The mass number of an isotope.

    Returns
    -------

    mass: Quantity
        The mass of the particle.

    Raises
    ------
    TypeError
        The argument is not a string, integer, or Quantity.

    InvalidParticleError
        If the argument does not correspond to a valid particle.

    MissingAtomicDataError
        If the standard atomic weight, the isotope mass, or the particle
        mass is not available.

    See also
    --------

    standard_atomic_weight : returns the conventional atomic mass of an
        element based on terrestrial values and assuming the atom is
        neutral.

    isotope_mass : returns the mass of an isotope (if available)
        assuming the atom is neutral.

    ion_mass : returns the mass of an ion of an element or isotope,
        accounting for reduction of mass from the neutral state due to
        different numbers of electrons.

    Notes
    -----
    This function will return the ion mass for ions, the isotope mass
    for isotopes (when available), the standard atomic weight for
    elements (when available), or the mass of special particles, as
    appropriate.

    The masses of neutrinos are not available because primarily upper
    limits are presently known.

    """
    return particle.mass


@particle_input
def isotopic_abundance(isotope: Particle,
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
    return isotope.isotopic_abundance


@particle_input
def integer_charge(particle: Particle) -> int:
    r"""Returns the integer charge of a particle.

    Parameters
    ----------

    particle : string
        String representing a particle.

    Returns
    -------

    Z : integer
        The charge as a multiple of the elementary charge.

    Raises
    ------

    InvalidParticleError
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    ChargeError
        If charge information for the particle is not available.

    AtomicWarning
        If the input represents an ion with an integer charge that is
        less than or equal to -3, which is unlikely to occur in nature.

    Notes
    -----

    This function supports two formats for integer charge information.

    The first format is a string that has information for the element
    or isotope at the beginning, a space in between, and the integer
    charge information in the form of an integer followed by a plus or
    minus sign, or a plus or minus sign followed by an integer.

    The second format is a string containing element information at
    the beginning, following by one or more plus or minus signs.

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
    return particle.integer_charge


@particle_input(must_be={'charged', 'uncharged'}, any=True)
def electric_charge(particle: Particle) -> Quantity:
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

    >>> electric_charge('p+')
    <<class 'astropy.constants.codata2014.EMCODATA2014'> name='Electron charge' value=1.6021766208e-19 uncertainty=9.8e-28 unit='C' reference='CODATA 2014'>
    >>> electric_charge('H-')
    <Quantity -1.60217662e-19 C>

    """
    return particle.charge


@particle_input
def is_stable(particle: Particle, mass_numb: int = None) -> bool:
    r"""Returns true for stable isotopes and particles and false otherwise.

    Parameters
    ----------

    particle: integer or string
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

    >>> is_stable("H-1")
    True
    >>> is_stable("tritium")
    False
    >>> is_stable("e-")
    True
    >>> is_stable("tau+")
    False

    """
    if particle.element and not particle.isotope:
        raise InvalidIsotopeError(
            "The input to is_stable must be either an isotope or a "
            "special particle."
        )
    return particle.is_category('stable')


@particle_input
def half_life(particle: Particle, mass_numb: int = None) -> Quantity:
    r"""Returns the half-life in seconds for unstable isotopes and
    particles, and numpy.inf in seconds for stable isotopes and particles.

    Parameters
    ----------

    particle: integer, string, or Particle
        A string representing an isotope or particle, an integer
        representing an atomic number, or an instance of the Particle
        class.

    mass_numb: integer
        The mass number of an isotope.

    Returns
    -------

    half_life_sec: astropy Quantity
        The half-life of the isotope or particle in units of seconds.

    Raises:
    -------

    InvalidParticleError
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    MissingAtomicDataError
        If no half-life data is available for the isotope.

    TypeError
        The argument is not an integer or string or the mass number is
        not an integer.

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
    return particle.half_life


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

    This list returns both natural and artificially produced isotopes.

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

    def common_isotopes_for_element(
            argument: Union[str, int],
            most_common_only: Optional[bool]) -> List[str]:

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


def _is_electron(arg: Any) -> bool:
    r"""Returns True if the argument corresponds to an electron, and False
    otherwise."""

    if not isinstance(arg, str):
        return False

    return arg in ['e', 'e-'] or arg.lower() == 'electron'


def periodic_table_period(argument: Union[str, int]) -> int:
    r"""Returns the periodic table period.

    Parameters
    ----------

    argument: string or integer
        Atomic number (either integer or string), atomic symbol (e.g. "H",
        string), or element name (e.g. "Francium", string).

    Returns
    -------

    period: integer
        The the periodic table period of the element.

    Raises
    ------

    TypeError:
        If the argument is not a string or integer.

    See also
    --------

    periodic_table_group : returns periodic table group of element.

    periodic_table_block : returns periodic table block of element.

    Examples
    --------

    >>> periodic_table_period(5)
    2
    >>> periodic_table_period("5")
    2
    >>> periodic_table_period("Au")
    6
    >>> periodic_table_period("nitrogen")
    2

    """
    if not isinstance(argument, (str, int)):
        raise TypeError("The argument to periodic_table_period must be " +
                        "either a string representing the element or its " +
                        "symbol, or an integer representing its atomic " +
                        "number.")
    symbol = atomic_symbol(argument)
    period = _Elements[symbol]["period"]
    return period


def periodic_table_group(argument: Union[str, int]) -> int:
    r"""Returns the periodic table group.

    Parameters
    ----------

    argument: string or integer
        Atomic number (either integer or string), atomic symbol (e.g. "H",
        string), or element name (e.g. "Francium", string).

    Returns
    -------

    group: integer
        The periodic table group of the element.

    Raises
    ------

    TypeError:
        If the argument is not a string or integer.

    See also
    --------

    periodic_table_period : returns periodic table period of element.

    periodic_table_block : returns periodic table block of element.

    Examples
    --------

    >>> periodic_table_group(18)
    18
    >>> periodic_table_group(24)
    6
    >>> periodic_table_group("Al")
    13
    >>> periodic_table_group("neon")
    18
    >>> periodic_table_group("BARIUM")
    2

    """
    if not isinstance(argument, (str, int)):
        raise TypeError("The argument to periodic_table_group must be " +
                        "either a string representing the element or its " +
                        "symbol, or an integer representing its atomic " +
                        "number.")
    symbol = atomic_symbol(argument)
    group = _Elements[symbol]["group"]
    return group


def periodic_table_block(argument: Union[str, int]) -> str:
    r"""Returns the periodic table block.

    Parameters
    ----------

    argument: string or integer
        Atomic number (either integer or string), atomic symbol (e.g. "H",
        string), or element name (e.g. "Francium", string).

    Returns
    -------

    block: string
        The periodic table block of the element.

    Raises
    ------

    TypeError:
        If the argument is not a string or integer.

    See also
    --------

    periodic_table_period : returns periodic table period of element.

    periodic_table_group : returns periodic table group of element.

    Examples
    --------

    >>> periodic_table_block(66)
    'f'
    >>> periodic_table_block(72)
    'd'
    >>> periodic_table_block("Tl")
    'p'
    >>> periodic_table_block("thallium")
    'p'
    >>> periodic_table_block("FRANCIUM")
    's'

    """
    if not isinstance(argument, (str, int)):
        raise TypeError("The argument to periodic_table_block must be " +
                        "either a string representing the element or its " +
                        "symbol, or an integer representing its atomic " +
                        "number.")
    symbol = atomic_symbol(argument)
    block = _Elements[symbol]["block"]
    return block


def periodic_table_category(argument: Union[str, int]) -> str:
    r"""Returns the periodic table category.

    Parameters
    ----------

    argument: string or integer
        Atomic number (either integer or string), atomic symbol (e.g. "H",
        string), or element name (e.g. "Francium", string).

    Returns
    -------

    category: string
        The periodic table category of the element.

    Raises
    ------

    TypeError:
        If the argument is not a string or integer.

    See also
    --------

    periodic_table_period : returns periodic table period of element.

    periodic_table_group : returns periodic table group of element.

    Examples
    --------

    >>> periodic_table_category(82)
    'Post-transition metals'
    >>> periodic_table_category("85")
    'Halogens'
    >>> periodic_table_category("Ra")
    'Alkaline earth metals'
    >>> periodic_table_category("rhodium")
    'Transition metals'

    """
    if not isinstance(argument, (str, int)):
        raise TypeError("The argument to periodic_table_category must be " +
                        "either a string representing the element or its " +
                        "symbol, or an integer representing its atomic " +
                        "number.")
    symbol = atomic_symbol(argument)
    category = _Elements[symbol]["category"]
    return category
