"""Functions that retrieve or are related to elemental or isotopic data."""

from numbers import Integral, Real
from typing import Any, List, Optional, Union

import astropy.constants as const
import astropy.units as u
import numpy as np

from plasmapy.particles.elements import _Elements
from plasmapy.particles.exceptions import (
    InvalidElementError,
    InvalidIsotopeError,
    InvalidParticleError,
    MissingAtomicDataError,
)
from plasmapy.particles.isotopes import _Isotopes
from plasmapy.particles.particle_class import Particle
from plasmapy.particles.particle_input import particle_input
from plasmapy.particles.symbols import atomic_symbol

__all__ = [
    "atomic_number",
    "mass_number",
    "standard_atomic_weight",
    "particle_mass",
    "isotopic_abundance",
    "integer_charge",
    "electric_charge",
    "is_stable",
    "half_life",
    "known_isotopes",
    "common_isotopes",
    "stable_isotopes",
    "reduced_mass",
    "periodic_table_period",
    "periodic_table_group",
    "periodic_table_block",
    "periodic_table_category",
]


@particle_input
def atomic_number(element: Particle) -> Integral:
    """
    Return the number of protons in an atom, isotope, or ion.

    Parameters
    ----------
    element: `str` or `~plasmapy.particles.Particle`
        A string representing an element, isotope, or ion; or an
        instance of the `~plasmapy.particles.Particle` class.

    Returns
    -------
    atomic_number: `int`
        The atomic number of an element.

    Raises
    ------
    `~plasmapy.utils.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    `TypeError`
        If the argument is not a `str`.

    See Also
    --------
    `~plasmapy.particles.mass_number` : returns the mass number (the total
        number of protons and neutrons) of an isotope.

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
def mass_number(isotope: Particle) -> Integral:
    """Get the mass number (the number of protons and neutrons) of an
    isotope.

    Parameters
    ----------
    isotope : `str` or `~plasmapy.particles.Particle`
        A string representing an isotope or a neutron; or an instance of
        the `plasmapy.particles.Particle` class.

    Returns
    -------
    mass_number : `int`
       The total number of protons plus neutrons in a nuclide.

    Raises
    ------
    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    `~plasmapy.utils.InvalidIsotopeError`
        If the argument does not correspond to a valid isotope.

    `TypeError`
        The argument is not a `str`.


    See Also
    --------
    `~plasmapy.particles.atomic_number` : returns the number of protons in
        an isotope or element

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


@particle_input(exclude={"isotope", "ion"})
def standard_atomic_weight(element: Particle) -> u.Quantity:
    """Return the standard (conventional) atomic weight of an element
    based on the relative abundances of isotopes in terrestrial
    environments.

    Parameters
    ----------
    element: `str`, `int`, or `~plasmapy.particles.Particle`
        A string representing an element or an integer representing an
        atomic number, or an instance of the Particle class.

    Returns
    -------
    atomic_weight: `~astropy.units.Quantity`
        The standard atomic weight of an element based on values from
        NIST.

    Raises
    ------
    `~plasmapy.utils.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    `TypeError`
        If the argument is not a `str` or `int`.

    See Also
    --------
    particle_mass

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
    >>> standard_atomic_weight("H")
    <Quantity 1.6738233e-27 kg>
    >>> standard_atomic_weight("lead")
    <Quantity 3.440636e-25 kg>

    """
    # TODO: Put in ReST links into above docstring
    return element.standard_atomic_weight


@particle_input(exclude={"neutrino", "antineutrino"})
def particle_mass(
    particle: Particle, *, Z: Integral = None, mass_numb: Integral = None
) -> u.Quantity:
    """
    Return the mass of a particle.

    Parameters
    ----------
    particle: `str`, `int`, or `~plasmapy.particles.Particle`
        A string representing an element, isotope, ion, or special
        particle; an integer representing an atomic number; or an
        instance of the Particle class.

    Z: `int`, optional, keyword-only
        The ionization state of the ion.

    mass_numb: `int`, optional, keyword-only
        The mass number of an isotope.

    Returns
    -------
    mass: `~astropy.units.Quantity`
        The mass of the particle.

    Raises
    ------
    `TypeError`
        The argument is not a string, integer, or Quantity.

    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    `~plasmapy.utils.MissingAtomicDataError`
        If the standard atomic weight, the isotope mass, or the particle
        mass is not available.

    See Also
    --------
    ~plasmapy.particles.standard_atomic_weight

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
def isotopic_abundance(isotope: Particle, mass_numb: Optional[Integral] = None) -> Real:
    """
    Return the isotopic abundances if known, and otherwise zero.

    Parameters
    ----------
    argument: `str` or `int`
        A string representing an element or isotope, or an integer
        representing the atomic number of an element.

    mass_numb: `int`, optional
        The mass number of an isotope, which is required if and only
        if the first argument can only be used.

    Returns
    -------
    iso_comp: `float`
        The relative isotopic abundance in the terrestrial environment.

    Raises
    ------
    `~plasmapy.utils.InvalidIsotopeError`
        If the argument is a valid particle but not a valid isotope.

    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    `TypeError`
        If the argument is not a `str` or `int`.

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

    """
    return isotope.isotopic_abundance


@particle_input(any_of={"charged", "uncharged"})
def integer_charge(particle: Particle) -> Integral:
    """Return the integer charge of a particle.

    Parameters
    ----------
    particle : `str`
        String representing a particle.

    Returns
    -------
    Z : `int`
        The charge as a multiple of the elementary charge.

    Raises
    ------
    `~plasmapy.particles.InvalidParticleError`
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    `~plasmapy.particles.ChargeError`
        If charge information for the particle is not available.

    `~plasmapy.particles.AtomicWarning`
        If the input represents an ion with an integer charge that is
        less than or equal to ``-3``, which is unlikely to occur in
        nature.

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


@particle_input(any_of={"charged", "uncharged"})
def electric_charge(particle: Particle) -> u.Quantity:
    """
    Return the electric charge (in coulombs) of a particle.

    Parameters
    ----------
    particle : `str`
        String representing an element or isotope followed by integer
        charge information.

    Returns
    -------
    charge: `~astropy.units.Quantity`
        The electric charge in coulombs.

    Raises
    ------
    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    `~plasmapy.utils.ChargeError`
        If charge information for the particle is not available.

    `~plasmapy.utils.AtomicWarning`
        If the input represents an ion with an integer charge that is
        below ``-3``.

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
    <<class 'astropy.constants.codata...'> name='Electron charge' ...>
    >>> electric_charge('H-')
    <Quantity -1.60217662e-19 C>

    """
    return particle.charge


@particle_input
def is_stable(particle: Particle, mass_numb: Optional[Integral] = None) -> bool:
    """
    Return `True` for stable isotopes and particles and `False` for
    unstable isotopes.

    Parameters
    ----------
    particle: `int`, `str`, or `~plasmapy.particles.Particle`
        A string representing an isotope or particle, or an integer
        representing an atomic number.

    mass_numb: `int`, optional
        The mass number of the isotope.

    Returns
    -------
    is_stable: `bool`
        `True` if the isotope is stable, `False` if it is unstable.

    Raises
    ------
    `~plasmapy.utils.InvalidIsotopeError`
        If the arguments correspond to a valid element but not a
        valid isotope.

    `~plasmapy.utils.InvalidParticleError`
        If the arguments do not correspond to a valid particle.

    `TypeError`
        If the argument is not a `str` or `int`.

    `~plasmapy.utils.MissingAtomicDataError`
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
            "The input to is_stable must be either an isotope or a special particle."
        )
    return particle.is_category("stable")


@particle_input(any_of={"stable", "unstable", "isotope"})
def half_life(particle: Particle, mass_numb: Optional[Integral] = None) -> u.Quantity:
    """
    Return the half-life in seconds for unstable isotopes and particles,
    and numpy.inf in seconds for stable isotopes and particles.

    Parameters
    ----------
    particle: `int`, `str`, or `~plasmapy.particles.Particle`
        A string representing an isotope or particle, an integer
        representing an atomic number, or an instance of the Particle
        class.

    mass_numb: `int`, optional
        The mass number of an isotope.

    Returns
    -------
    half_life_sec: `~astropy.units.Quantity`
        The half-life of the isotope or particle in units of seconds.

    Raises
    ------
    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    `~plasmapy.utils.MissingAtomicDataError`
        If no half-life data is available for the isotope.

    `TypeError`
        The argument is not an `int` or `str` or the mass number is
        not an `int`.

    Notes
    -----
    Accurate half-life data is not known for all isotopes. Some isotopes
    may have upper or lower limits on the half-life, in which case this
    function will return a string with that information and issue a
    `~plasmapy.utils.MissingAtomicDataWarning`.  When no isotope
    information is available, then this function raises a
    `~plasmapy.utils.MissingAtomicDataError`.

    Examples
    --------
    >>> half_life('T')
    <Quantity 3.888e+08 s>
    >>> half_life('n')
    <Quantity 881.5 s>
    >>> half_life('H-1')
    <Quantity inf s>

    """
    return particle.half_life


def known_isotopes(argument: Union[str, Integral] = None) -> List[str]:
    """Return a list of all known isotopes of an element, or a list
    of all known isotopes of every element if no input is provided.

    Parameters
    ----------
    argument: `int` or `str`, optional
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    Returns
    -------
    isotopes_list: `list` containing `str` items or an empty `list`
        List of all of the isotopes of an element that have been
        discovered, sorted from lowest mass number to highest mass
        number.  If no argument is provided, then a list of all known
        isotopes of every element will be returned that is sorted by
        atomic number, with entries for each element sorted by mass
        number.

    Raises
    ------
    `~plasmapy.utils.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    `TypeError`
        If the argument is not a `str` or `int`.

    Notes
    -----
    This list returns both natural and artificially produced isotopes.

    See Also
    --------
    `~plasmapy.particles.common_isotopes` : returns isotopes with non-zero
        isotopic abundances.

    `~plasmapy.particles.stable_isotopes` : returns isotopes that are
        stable against radioactive decay.

    Examples
    --------
    >>> known_isotopes('H')
    ['H-1', 'D', 'T', 'H-4', 'H-5', 'H-6', 'H-7']
    >>> known_isotopes('helium 1+')
    ['He-3', 'He-4', 'He-5', 'He-6', 'He-7', 'He-8', 'He-9', 'He-10']
    >>> known_isotopes()[0:10]
    ['H-1', 'D', 'T', 'H-4', 'H-5', 'H-6', 'H-7', 'He-3', 'He-4', 'He-5']
    >>> len(known_isotopes())  # the number of known isotopes
    3352

    """

    # TODO: Allow Particle objects representing elements to be inputs

    def known_isotopes_for_element(argument):
        element = atomic_symbol(argument)
        isotopes = []
        for isotope in _Isotopes.keys():
            if element + "-" in isotope and isotope[0 : len(element)] == element:
                isotopes.append(isotope)
        if element == "H":
            isotopes.insert(1, "D")
            isotopes.insert(2, "T")
        mass_numbers = [mass_number(isotope) for isotope in isotopes]
        sorted_isotopes = [
            mass_number
            for (isotope, mass_number) in sorted(zip(mass_numbers, isotopes))
        ]
        return sorted_isotopes

    if argument is not None:
        try:
            element = atomic_symbol(argument)
            isotopes_list = known_isotopes_for_element(element)
        except InvalidElementError:
            raise InvalidElementError(
                "known_isotopes is unable to get "
                f"isotopes from an input of: {argument}"
            )
        except InvalidParticleError:
            raise InvalidParticleError("Invalid particle in known_isotopes.")
    elif argument is None:
        isotopes_list = []
        for atomic_numb in range(1, len(_Elements.keys()) + 1):
            isotopes_list += known_isotopes_for_element(atomic_numb)

    return isotopes_list


def common_isotopes(
    argument: Union[str, Integral] = None, most_common_only: bool = False
) -> List[str]:
    """
    Return a list of isotopes of an element with an isotopic abundances
    greater than zero, or if no input is provided, a list of all such
    isotopes for every element.

    Parameters
    ----------
    argument: `int` or `str`, optional
        A string or integer representing an atomic number or element,
        or a string representing an isotope.

    most_common_only: `bool`
        If set to `True`, return only the most common isotope.

    Returns
    -------
    isotopes_list: `list` of `str` or empty `list`
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
    `~plasmapy.utils.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    `TypeError`
        If the argument is not a string or integer.

    Notes
    -----
    The isotopic abundances are based on the terrestrial environment
    and may not be appropriate for space and astrophysical applications.

    See Also
    --------
    `~plasmapy.utils.known_isotopes` : returns a list of isotopes that
        have been discovered.

    `~plasmapy.utils.stable_isotopes` : returns isotopes that are stable
        against radioactive decay.

    `~plasmapy.utils.isotopic_abundance` : returns the relative isotopic
         abundance.

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

    # TODO: Allow Particle objects representing elements to be inputs

    def common_isotopes_for_element(
        argument: Union[str, int], most_common_only: Optional[bool]
    ) -> List[str]:

        isotopes = known_isotopes(argument)

        CommonIsotopes = [
            isotope for isotope in isotopes if "abundance" in _Isotopes[isotope].keys()
        ]

        isotopic_abundances = [
            _Isotopes[isotope]["abundance"] for isotope in CommonIsotopes
        ]

        sorted_isotopes = [
            iso_comp
            for (isotope, iso_comp) in sorted(zip(isotopic_abundances, CommonIsotopes))
        ]

        sorted_isotopes.reverse()

        if most_common_only and len(sorted_isotopes) > 1:
            sorted_isotopes = sorted_isotopes[0:1]

        return sorted_isotopes

    if argument is not None:

        try:
            element = atomic_symbol(argument)
            isotopes_list = common_isotopes_for_element(element, most_common_only)
        except InvalidParticleError:
            raise InvalidParticleError("Invalid particle")
        except InvalidElementError:
            raise InvalidElementError(
                "common_isotopes is unable to get isotopes "
                f"from an input of: {argument}"
            )

    elif argument is None:
        isotopes_list = []
        for atomic_numb in range(1, 119):
            isotopes_list += common_isotopes_for_element(atomic_numb, most_common_only)

    return isotopes_list


def stable_isotopes(
    argument: Union[str, Integral] = None, unstable: bool = False
) -> List[str]:
    """
    Return a list of all stable isotopes of an element, or if no input is
    provided, a list of all such isotopes for every element.

    Parameters
    ----------
    argument: `int` or `str`
        A string or integer representing an atomic number or element,
        or a string representing an isotope.

    unstable: `bool`
        If set to `True`, this function will return a list of the
        unstable isotopes instead of the stable isotopes.

    Returns
    -------
    StableIsotopes: `list` of strings or empty list
        List of all stable isotopes of an element, sorted from lowest
        mass number.  If an element has no stable isotopes, this
        function returns an empty list.

    Raises
    ------
    `~plasmapy.utils.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    `TypeError`
        If the argument is not a string or integer.

    Notes
    -----
    There are 254 isotopes for which no radioactive decay has been
    observed.  It is possible that some isotopes will be discovered to
    be unstable but with extremely long half-lives.  For example,
    bismuth-209 was recently discovered to have a half-life of about
    1.9e19 years.  However, such isotopes can be regarded as virtually
    stable for most applications.

    See Also
    --------
    `~plasmapy.particles.known_isotopes` : returns a list of isotopes that
        have been discovered

    `~plasmapy.particles.common_isotopes` : returns isotopes with non-zero
        isotopic abundances

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

    Find unstable isotopes using the ``unstable`` keyword.

    >>> stable_isotopes('U', unstable=True)[:5]  # only first five
    ['U-217', 'U-218', 'U-219', 'U-220', 'U-221']

    """

    # TODO: Allow Particle objects representing elements to be inputs

    def stable_isotopes_for_element(
        argument: Union[str, int], stable_only: Optional[bool]
    ) -> List[str]:
        KnownIsotopes = known_isotopes(argument)
        StableIsotopes = [
            isotope
            for isotope in KnownIsotopes
            if _Isotopes[isotope]["stable"] == stable_only
        ]
        return StableIsotopes

    if argument is not None:
        try:
            element = atomic_symbol(argument)
            isotopes_list = stable_isotopes_for_element(element, not unstable)
        except InvalidParticleError:
            raise InvalidParticleError("Invalid particle in stable_isotopes")
        except InvalidElementError:
            raise InvalidElementError(
                "stable_isotopes is unable to get isotopes "
                f"from an input of: {argument}"
            )
    elif argument is None:
        isotopes_list = []
        for atomic_numb in range(1, 119):
            isotopes_list += stable_isotopes_for_element(atomic_numb, not unstable)

    return isotopes_list


def reduced_mass(test_particle, target_particle) -> u.Quantity:
    """
    Find the reduced mass between two particles.

    Parameters
    ----------
    test_particle, target_particle : `str`, `int`, `~plasmapy.particles.Particle`,
    `~astropy.units.Quantity`, or `~astropy.constants.Constant`

        The test particle as represented by a string, an integer
        representing atomic number, a `~plasmapy.particles.Particle`
        object, or a `~astropy.units.Quantity` or
        `~astropy.constants.Constant` with units of mass.

    Return
    -------
    reduced_mass : `~astropy.units.Quantity`
        The reduced mass between the test particle and target particle.

    Raises
    ------
    `~plasmapy.utils.InvalidParticleError`
        If either particle is invalid.

    `~astropy.units.UnitConversionError`
        If an argument is a `~astropy.units.Quantity` or
        `~astropy.units.Constant` but does not have units of mass.

    `~plasmapy.utils.MissingAtomicDataError`
        If the mass of either particle is not known.

    `TypeError`
        If either argument is not a `str`, `int`,
        `~plasmapy.particles.Particle`, `~astropy.units.Quantity`, or
        `~astropy.constants.Constant`.

    Example
    -------
    >>> from astropy import units as u
    >>> reduced_mass('p+', 'e-')
    <Quantity 9.104425e-31 kg>
    >>> reduced_mass(5.4e-27 * u.kg, 8.6e-27 * u.kg)
    <Quantity 3.31714286e-27 kg>

    """

    # TODO: Add discussion on reduced mass and its importance to docstring
    # TODO: Add equation for reduced mass to docstring

    def get_particle_mass(particle) -> u.Quantity:
        """Return the mass of a particle.

        Take a representation of a particle and returns the mass in
        kg.  If the input is a `~astropy.units.Quantity` or
        `~astropy.constants.Constant` with units of mass already, then
        this returns that mass converted to kg.
        """
        try:
            if isinstance(particle, (u.Quantity, const.Constant)):
                return particle.to(u.kg)
            if not isinstance(particle, Particle):
                particle = Particle(particle)
            return particle.mass.to(u.kg)
        except u.UnitConversionError as exc1:
            raise u.UnitConversionError(f"Incorrect units in reduced_mass.") from exc1
        except MissingAtomicDataError:
            raise MissingAtomicDataError(
                f"Unable to find the reduced mass because the mass of "
                f"{particle} is not available."
            ) from None

    test_mass = get_particle_mass(test_particle)
    target_mass = get_particle_mass(target_particle)

    return (test_mass * target_mass) / (test_mass + target_mass)


def periodic_table_period(argument: Union[str, Integral]) -> Integral:
    """
    Return the periodic table period.

    Parameters
    ----------
    argument: `str` or `int`
        Atomic number (either integer or string), atomic symbol (e.g.,
        ``"H"``, string), or element name (e.g. ``"Francium"``, string).

    Returns
    -------
    period: `int`
        The periodic table period of the element.

    Raises
    ------
    `TypeError`
        If the argument is not a `str` or `int`.

    See Also
    --------
    `~plasmapy.particles.periodic_table_group` : returns periodic table
        group of element.

    `~plasmapy.particles.periodic_table_block` : returns periodic table
        block of element.

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
    # TODO: Implement @particle_input
    if not isinstance(argument, (str, Integral)):
        raise TypeError(
            "The argument to periodic_table_period must be either a "
            "string representing the element or its symbol, or an "
            "integer representing its atomic number."
        )
    symbol = atomic_symbol(argument)
    period = _Elements[symbol]["period"]
    return period


def periodic_table_group(argument: Union[str, Integral]) -> Integral:
    """
    Return the periodic table group.

    Parameters
    ----------
    argument: `str` or `int`
        Atomic number (either integer or string), atomic symbol (e.g.,
        ``"H"``, string), or element name (e.g., ``"francium"``,
        string).

    Returns
    -------
    group: `int`
        The periodic table group of the element.

    Raises
    ------
    `TypeError`
        If the argument is not a `str` or `int`.

    See Also
    --------
    `~plasmapy.particles.periodic_table_period` : returns periodic table
        period of element.

    `~plasmapy.particles.periodic_table_block` : returns periodic table
        block of element.

    `~plasmapy.particles.periodic_table_category` : returns periodic table
        category of element.

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
    >>> periodic_table_group("barium")
    2

    """
    # TODO: Implement @particle_input
    if not isinstance(argument, (str, Integral)):
        raise TypeError(
            "The argument to periodic_table_group must be "
            "either a string representing the element or its "
            "symbol, or an integer representing its atomic number."
        )
    symbol = atomic_symbol(argument)
    group = _Elements[symbol]["group"]
    return group


def periodic_table_block(argument: Union[str, Integral]) -> str:
    """
    Return the periodic table block.

    Parameters
    ----------
    argument: `str` or `int`
        Atomic number (either integer or string), atomic symbol (e.g.,
        ``"H"``, string), or element name (e.g., ``"francium"``,
        string).

    Returns
    -------
    block: `str`
        The periodic table block of the element.

    Raises
    ------
    `TypeError`
        If the argument is not a `str` or `int`.

    See Also
    --------
    `~plasmapy.particles.periodic_table_period` : returns periodic table
        period of element.

    `~plasmapy.particles.periodic_table_group` : returns periodic table
        group of element.

    `~plasmapy.particles.periodic_table_category` : returns periodic table
        category of element.

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
    >>> periodic_table_block("francium")
    's'

    """
    # TODO: Implement @particle_input
    if not isinstance(argument, (str, Integral)):
        raise TypeError(
            "The argument to periodic_table_block must be "
            "either a string representing the element or its "
            "symbol, or an integer representing its atomic number."
        )
    symbol = atomic_symbol(argument)
    block = _Elements[symbol]["block"]
    return block


def periodic_table_category(argument: Union[str, Integral]) -> str:
    """
    Return the periodic table category.

    Parameters
    ----------
    argument: `str` or `int`
        Atomic number (either integer or string), atomic symbol (e.g.,
        ``"H"``, string), or element name (e.g., ``"francium"``,
        string).

    Returns
    -------
    category: `str`
        The periodic table category of the element.

    Raises
    ------
    `TypeError`
        If the argument is not a `str` or `int`.

    See Also
    --------
    `~plasmapy.particles.periodic_table_period` : returns periodic table
        period of element.

    `~plasmapy.particles.periodic_table_group` : returns periodic table
        group of element.

    `~plasmapy.particles.periodic_table_block` : returns periodic table
        block of element.

    Examples
    --------
    >>> periodic_table_category(82)
    'post-transition metal'
    >>> periodic_table_category("85")
    'halogen'
    >>> periodic_table_category("Ra")
    'alkaline earth metal'
    >>> periodic_table_category("rhodium")
    'transition metal'

    """
    # TODO: Implement @particle_input
    if not isinstance(argument, (str, Integral)):
        raise TypeError(
            "The argument to periodic_table_category must be "
            "either a string representing the element or its "
            "symbol, or an integer representing its atomic number."
        )
    symbol = atomic_symbol(argument)
    category = _Elements[symbol]["category"]
    return category


def _is_electron(arg: Any) -> bool:
    """
    Return `True` if the argument corresponds to an electron, and
    `False` otherwise.
    """
    # TODO: Remove _is_electron from all parts of code.

    if not isinstance(arg, str):
        return False

    return arg in ["e", "e-"] or arg.lower() == "electron"
