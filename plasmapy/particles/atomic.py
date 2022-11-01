"""Functions that retrieve or are related to elemental or isotopic data."""

__all__ = [
    "atomic_number",
    "charge_number",
    "common_isotopes",
    "electric_charge",
    "half_life",
    "ionic_levels",
    "isotopic_abundance",
    "is_stable",
    "known_isotopes",
    "mass_number",
    "particle_mass",
    "periodic_table_block",
    "periodic_table_category",
    "periodic_table_group",
    "periodic_table_period",
    "reduced_mass",
    "stable_isotopes",
    "standard_atomic_weight",
]

import astropy.constants as const
import astropy.units as u

from numbers import Integral, Real
from typing import Any, List, Optional, Union

from plasmapy.particles import _elements, _isotopes
from plasmapy.particles.decorators import particle_input
from plasmapy.particles.exceptions import (
    ChargeError,
    InvalidElementError,
    InvalidIsotopeError,
    InvalidParticleError,
    MissingParticleDataError,
)
from plasmapy.particles.particle_class import Particle, ParticleLike
from plasmapy.particles.particle_collections import ParticleList
from plasmapy.particles.symbols import atomic_symbol

__all__.sort()


@particle_input
def atomic_number(element: Particle) -> Integral:
    """
    Return the number of protons in an atom, isotope, or ion.

    Parameters
    ----------
    element : |atom-like|
        A string representing an element, isotope, or ion; or an
        instance of the `~plasmapy.particles.particle_class.Particle` class.

    Returns
    -------
    `int`
        The atomic number of an element.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    See Also
    --------
    ~plasmapy.particles.atomic.mass_number

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
    isotope : |atom-like|
        A string representing an isotope or a neutron; or an instance of
        the `plasmapy.particles.particle_class.Particle` class.

    Returns
    -------
    `int`
       The total number of protons plus neutrons in a nuclide.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    `~plasmapy.particles.exceptions.InvalidIsotopeError`
        If the argument does not correspond to a valid isotope.

    See Also
    --------
    ~plasmapy.particles.atomic.atomic_number

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
    element : |atom-like|
        A string representing an element or an integer representing an
        atomic number, or an instance of the Particle class.

    Returns
    -------
    `~astropy.units.Quantity`
        The standard atomic weight of an element based on values from
        NIST.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    See Also
    --------
    ~plasmapy.particles.atomic.particle_mass

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
    particle: Particle,
    *,
    mass_numb: Optional[Integral] = None,
    Z: Optional[Integral] = None,
) -> u.kg:
    """
    Return the mass of a particle.

    Parameters
    ----------
    particle : |particle-like|
        A string representing an element, isotope, ion, or special
        particle; an integer representing an atomic number; or a
        |Particle|.

    mass_numb : integer, optional, |keyword-only|
        The mass number of an isotope.

    Z : integer, optional, |keyword-only|
        The |charge number| of an ion or neutral atom.

    Returns
    -------
    `~astropy.units.Quantity`
        The mass of the particle.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    `~plasmapy.particles.exceptions.MissingParticleDataError`
        If the standard atomic weight, the isotope mass, or the particle
        mass is not available.

    See Also
    --------
    ~plasmapy.particles.atomic.standard_atomic_weight

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
    isotope : |atom-like|
        A string representing an element or isotope, or an integer
        representing the atomic number of an element.

    mass_numb : integer, optional
        The mass number of an isotope.

    Returns
    -------
    `float`
        The relative isotopic abundance in the terrestrial environment.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidIsotopeError`
        If the argument is a valid particle but not a valid isotope.

    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

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
def charge_number(particle: Particle) -> Integral:
    """Return the charge number of a particle.

    Parameters
    ----------
    particle : |particle-like|
        String representing a particle.

    Returns
    -------
    `int`
        The charge as a multiple of the elementary charge.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    `~plasmapy.particles.exceptions.ChargeError`
        If charge information for the particle is not available.

    `~plasmapy.particles.exceptions.ParticleWarning`
        If the input represents an ion with a charge number that is
        less than or equal to ``-3``, which is unlikely to occur in
        nature.

    Notes
    -----
    This function supports two formats for charge number information.

    The first format is a string that has information for the element
    or isotope at the beginning, a space in between, and the charge
    number information in the form of an integer followed by a plus or
    minus sign, or a plus or minus sign followed by an integer.

    The second format is a string containing element information at
    the beginning, following by one or more plus or minus signs.

    Examples
    --------
    >>> charge_number('Fe-56 2+')
    2
    >>> charge_number('He -2')
    -2
    >>> charge_number('H+')
    1
    >>> charge_number('N-14++')
    2
    """
    return particle.charge_number


@particle_input(any_of={"charged", "uncharged"})
def electric_charge(particle: Particle) -> u.C:
    """
    Return the electric charge (in coulombs) of a particle.

    Parameters
    ----------
    particle : |particle-like|
        String representing an element or isotope followed by integer
        charge information.

    Returns
    -------
    `~astropy.units.Quantity`
        The electric charge in coulombs.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    `~plasmapy.particles.exceptions.ChargeError`
        If charge information for the particle is not available.

    `~plasmapy.particles.exceptions.ParticleWarning`
        If the input represents an ion with a charge number that is
        below ``-3``.

    Notes
    -----
    This function supports two formats for charge number information.

    The first format is a string that has information for the element
    or isotope at the beginning, a space in between, and the charge
    number information in the form of an integer followed by a plus or
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
    particle : |particle-like|
        A string representing an isotope or particle, or an integer
        representing an atomic number.

    mass_numb : integer, optional
        The mass number of the isotope.

    Returns
    -------
    `bool`
        `True` if the isotope is stable, `False` if it is unstable.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidIsotopeError`
        If the arguments correspond to a valid element but not a
        valid isotope.

    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the arguments do not correspond to a valid particle.

    `~plasmapy.particles.exceptions.MissingParticleDataError`
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
    and |inf| seconds for stable isotopes and particles.

    Parameters
    ----------
    particle : |particle-like|
        A string representing an isotope or particle, an integer
        representing an atomic number, or an instance of the |Particle|
        class.

    mass_numb : integer, optional
        The mass number of an isotope.

    Returns
    -------
    `~astropy.units.Quantity`
        The half-life of the isotope or particle in units of seconds.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    `~plasmapy.particles.exceptions.MissingParticleDataError`
        If no half-life data is available for the isotope.

    Notes
    -----
    Accurate half-life data is not known for all isotopes. Some isotopes
    may have upper or lower limits on the half-life, in which case this
    function will return a string with that information and issue a
    `~plasmapy.particles.exceptions.MissingParticleDataError`.  When no
    isotope information is available, then this function raises a
    `~plasmapy.particles.exceptions.MissingParticleDataError`.

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
    """
    Return a list of all known isotopes of an element, or a list of all
    known isotopes of every element if no input is provided.

    Parameters
    ----------
    argument : |atom-like|
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    Returns
    -------
    `list` of `str`
        List of all the isotopes of an element that have been
        discovered, sorted from low to high mass number. If no argument
        is provided, then a list of all known isotopes of every element
        will be returned that is sorted first by atomic number and
        second by mass number.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    Notes
    -----
    This list returns both natural and artificially produced isotopes.

    See Also
    --------
    ~plasmapy.particles.atomic.common_isotopes :
        Returns isotopes with non-zero isotopic abundances.

    ~plasmapy.particles.atomic.stable_isotopes :
        Returns isotopes that are stable against radioactive decay.

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
        isotopes = [
            isotope
            for isotope in _isotopes.data_about_isotopes
            if f"{element}-" in isotope and isotope[: len(element)] == element
        ]

        if element == "H":
            isotopes.insert(1, "D")
            isotopes.insert(2, "T")
        mass_numbers = [mass_number(isotope) for isotope in isotopes]
        return [
            mass_number
            for (isotope, mass_number) in sorted(zip(mass_numbers, isotopes))
        ]

    if argument is not None:
        try:
            element = atomic_symbol(argument)
            isotopes_list = known_isotopes_for_element(element)
        except InvalidElementError as ex:
            raise InvalidElementError(
                "known_isotopes is unable to get "
                f"isotopes from an input of: {argument}"
            ) from ex
        except InvalidParticleError as ex:
            raise InvalidParticleError("Invalid particle in known_isotopes.") from ex
    elif argument is None:
        isotopes_list = []
        for atomic_numb in range(1, len(_elements.data_about_elements) + 1):
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
    argument : |atom-like|, optional
        A string or integer representing an atomic number or element,
        or a string representing an isotope.

    most_common_only : `bool`
        If set to `True`, return only the most common isotope.

    Returns
    -------
    `list` of `str`
        List of all isotopes of an element with isotopic abundances
        greater than zero, sorted from most abundant to least
        abundant.  If no isotopes have isotopic abundances greater
        than zero, this function will return an empty list.  If no
        arguments are provided, then a list of all common isotopes of
        all elements will be provided that is sorted first by low to
        high atomic number and second by most abundant to least abundant
        isotope.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    Notes
    -----
    The isotopic abundances are based on the terrestrial environment
    and may not be appropriate for space and astrophysical applications.

    See Also
    --------
    ~plasmapy.particles.atomic.known_isotopes :
        Returns a list of isotopes that have been discovered.

    ~plasmapy.particles.atomic.stable_isotopes :
        Returns isotopes that are stable against radioactive decay.

    ~plasmapy.particles.atomic.isotopic_abundance :
        Returns the relative isotopic abundance.

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
            isotope
            for isotope in isotopes
            if "abundance" in _isotopes.data_about_isotopes[isotope]
        ]

        isotopic_abundances = [
            _isotopes.data_about_isotopes[isotope]["abundance"]
            for isotope in CommonIsotopes
        ]

        sorted_isotopes = [
            iso_comp
            for (isotope, iso_comp) in sorted(zip(isotopic_abundances, CommonIsotopes))
        ]

        sorted_isotopes.reverse()

        if most_common_only and len(sorted_isotopes) > 1:
            sorted_isotopes = sorted_isotopes[:1]

        return sorted_isotopes

    if argument is not None:

        try:
            element = atomic_symbol(argument)
            isotopes_list = common_isotopes_for_element(element, most_common_only)
        except InvalidParticleError as ex:
            raise InvalidParticleError("Invalid particle") from ex
        except InvalidElementError as ex:
            raise InvalidElementError(
                "common_isotopes is unable to get isotopes "
                f"from an input of: {argument}"
            ) from ex

    elif argument is None:
        isotopes_list = []
        for atomic_numb in range(1, 119):
            isotopes_list += common_isotopes_for_element(atomic_numb, most_common_only)

    return isotopes_list


def stable_isotopes(
    argument: Optional[ParticleLike] = None, unstable: bool = False
) -> List[str]:
    """
    Return a list of all stable isotopes of an element, or if no input is
    provided, a list of all such isotopes for every element.

    Parameters
    ----------
    argument : |atom-like|
        A string or integer representing an atomic number or element,
        or a string representing an isotope.

    unstable : `bool`
        If set to `True`, this function will return a list of the
        unstable isotopes instead of the stable isotopes.

    Returns
    -------
    `list` of `str`
        List of all stable isotopes of an element, sorted from low to
        high mass number.  If an element has no stable isotopes, this
        function returns an empty list.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    Notes
    -----
    There are 254 isotopes for which no radioactive decay has been
    observed.  It is possible that some isotopes will be discovered to
    be unstable but with extremely long half-lives.  For example,
    bismuth-209 was recently discovered to have a half-life of about
    :math:`1.9 × 10^{19}` years.  However, such isotopes can be regarded
    as virtually stable for most applications.

    See Also
    --------
    ~plasmapy.particles.atomic.known_isotopes :
        Returns a list of isotopes that have been discovered.

    ~plasmapy.particles.atomic.common_isotopes :
        Returns isotopes with non-zero isotopic abundances.

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
        return [
            isotope
            for isotope in KnownIsotopes
            if _isotopes.data_about_isotopes[isotope]["stable"] == stable_only
        ]

    if argument is not None:
        try:
            element = atomic_symbol(argument)
            isotopes_list = stable_isotopes_for_element(element, not unstable)
        except InvalidParticleError as ex:
            raise InvalidParticleError("Invalid particle in stable_isotopes") from ex
        except InvalidElementError as ex:
            raise InvalidElementError(
                "stable_isotopes is unable to get isotopes "
                f"from an input of: {argument}"
            ) from ex
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
    test_particle, target_particle : |particle-like| or |Quantity|
        The test particle as represented by a string, an integer
        representing atomic number, a `~plasmapy.particles.particle_class.Particle`
        object, or a `~astropy.units.Quantity` or
        `~astropy.constants.Constant` with units of mass.

    Returns
    -------
    `~astropy.units.Quantity`
        The reduced mass between the test particle and target particle.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidParticleError`
        If either particle is invalid.

    `~astropy.units.UnitConversionError`
        If an argument is a `~astropy.units.Quantity` or
        `~astropy.constants.Constant` but does not have units of mass.

    `~plasmapy.particles.exceptions.MissingParticleDataError`
        If the mass of either particle is not known.

    Examples
    --------
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

        except u.UnitConversionError as exc1:
            raise u.UnitConversionError("Incorrect units in reduced_mass.") from exc1
        except MissingParticleDataError:
            raise MissingParticleDataError(
                f"Unable to find the reduced mass because the mass of "
                f"{particle} is not available."
            ) from None
        else:
            return particle.mass.to(u.kg)

    test_mass = get_particle_mass(test_particle)
    target_mass = get_particle_mass(target_particle)

    return (test_mass * target_mass) / (test_mass + target_mass)


def periodic_table_period(argument: Union[str, Integral]) -> Integral:
    """
    Return the periodic table period.

    Parameters
    ----------
    argument : |atom-like|
        Atomic number (either integer or string), atomic symbol (e.g.,
        ``"H"``, string), or element name (e.g. ``"Francium"``, string).

    Returns
    -------
    `int`
        The periodic table period of the element.

    See Also
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
    # TODO: Implement @particle_input
    if not isinstance(argument, (str, Integral)):
        raise TypeError(
            "The argument to periodic_table_period must be either a "
            "string representing the element or its symbol, or an "
            "integer representing its atomic number."
        )
    symbol = atomic_symbol(argument)
    return _elements.data_about_elements[symbol]["period"]


def periodic_table_group(argument: Union[str, Integral]) -> Integral:
    """
    Return the periodic table group.

    Parameters
    ----------
    argument : |atom-like|
        Atomic number (either integer or string), atomic symbol (e.g.,
        ``"H"``, string), or element name (e.g., ``"francium"``,
        string).

    Returns
    -------
    `int`
        The periodic table group of the element.

    See Also
    --------
    periodic_table_period : returns periodic table period of element.

    periodic_table_block : returns periodic table block of element.

    periodic_table_category : returns periodic table category of element.

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
    return _elements.data_about_elements[symbol]["group"]


def periodic_table_block(argument: Union[str, Integral]) -> str:
    """
    Return the periodic table block.

    Parameters
    ----------
    argument : |atom-like|
        Atomic number (either integer or string), atomic symbol (e.g.,
        ``"H"``, string), or element name (e.g., ``"francium"``,
        string).

    Returns
    -------
    `str`
        The periodic table block of the element.

    See Also
    --------
    ~plasmapy.particles.atomic.periodic_table_period :
        Returns periodic table period of element.

    ~plasmapy.particles.atomic.periodic_table_group :
        Returns periodic table group of element.

    ~plasmapy.particles.atomic.periodic_table_category :
        Returns periodic table category of element.

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
    return _elements.data_about_elements[symbol]["block"]


def periodic_table_category(argument: Union[str, Integral]) -> str:
    """
    Return the periodic table category.

    Parameters
    ----------
    argument : |atom-like|
        Atomic number (either integer or string), atomic symbol (e.g.,
        ``"H"``, string), or element name (e.g., ``"francium"``,
        string).

    Returns
    -------
    `str`
        The periodic table category of the element.

    See Also
    --------
    periodic_table_period : returns periodic table period of element.

    periodic_table_group : returns periodic table group of element.

    periodic_table_block : returns periodic table block of element.

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
    return _elements.data_about_elements[symbol]["category"]


@particle_input(any_of={"element", "isotope", "ion"})
def ionic_levels(
    particle: Particle,
    min_charge: Integral = 0,
    max_charge: Optional[Integral] = None,
) -> ParticleList:
    """
    Return a |ParticleList| that includes different ionic levels of a
    base atom.

    Parameters
    ----------
    particle : |atom-like|
        Representation of an element, ion, or isotope.

    min_charge : integer, optional, default: ``0``
        The starting charge number.

    max_charge : integer, optional
        The ending charge number, which will be included in the
        |ParticleList|.  Defaults to the atomic number of ``particle``.

    Returns
    -------
    `~plasmapy.particles.particle_collections.ParticleList`
        The ionic levels of the atom provided from ``min_charge`` to
        ``max_charge``.

    Examples
    --------
    >>> from plasmapy.particles import ionic_levels
    >>> ionic_levels("He")
    ParticleList(['He 0+', 'He 1+', 'He 2+'])
    >>> ionic_levels("Fe-56", min_charge=13, max_charge=15)
    ParticleList(['Fe-56 13+', 'Fe-56 14+', 'Fe-56 15+'])
    """
    base_particle = Particle(particle.isotope or particle.element)

    if max_charge is None:
        max_charge = particle.atomic_number

    if not min_charge <= max_charge <= particle.atomic_number:
        raise ChargeError(
            f"Need min_charge ({min_charge}) "
            f"≤ max_charge ({max_charge}) "
            f"≤ atomic number ({base_particle.atomic_number})."
        )

    return ParticleList(
        [Particle(base_particle, Z=Z) for Z in range(min_charge, max_charge + 1)]
    )


def _is_electron(arg: Any) -> bool:
    """
    Return `True` if the argument corresponds to an electron, and
    `False` otherwise.
    """
    # TODO: Remove _is_electron from all parts of code.

    return (
        arg in ["e", "e-"] or arg.lower() == "electron"
        if isinstance(arg, str)
        else False
    )
