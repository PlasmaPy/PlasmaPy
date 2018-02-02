r"""The Particle class and particle_input decorator."""

import numpy as np
import warnings
from typing import (Union, Set, Tuple, List, Optional)
from astropy import units as u, constants as const
import collections

from ..utils import (
    AtomicError,
    AtomicWarning,
    InvalidParticleError,
    InvalidElementError,
    InvalidIsotopeError,
    ChargeError,
    MissingAtomicDataError,
)

from .parsing import (
    _dealias_particle_aliases,
    _parse_and_check_atomic_input,
    _invalid_particle_errmsg,
)

from .elements import _Elements
from .isotopes import _Isotopes

from .special_particles import _Particles, ParticleZoo


class Particle:
    r"""A class for individual particles or antiparticles.

    Parameters
    ----------
    argument : str or int
        A string representing a particle, element, isotope, or ion; or
        an integer representing the atomic number of an element.

    mass_numb : int, optional
        The mass number of an isotope or nuclide.

    Z : int, optional
        The integer charge of the particle.

    Raises
    ------
    InvalidParticleError
        Raised when the particle input does not correspond to a valid particle
        or is contradictory.

    InvalidElementError
        For when an attribute is being accessed that requires information
        about an element, but the particle is not an element, isotope, or ion.

    InvalidIsotopeError
        For when an attribute is being accessed that requires information
        about an isotope or nuclide, but the particle is not an isotope (or
        an ion of an isotope).

    ChargeError
        For when either the charge or integer_charge attributes is being
        accessed but the charge information for the particle is not
        available.

    TypeError
        For when any of the arguments or keywords is not of the required
        type.

    Examples
    --------

    >>> proton = Particle('p+')
    >>> electron = Particle('e-')
    >>> neutron = Particle('neutron')
    >>> deuteron = Particle('D', Z=1)
    >>> alpha = Particle('He', mass_numb=2, Z=2)
    >>> positron = Particle('positron')
    >>> proton.element
    'H'
    >>> alpha.isotope
    'He-4'
    >>> deuterium.ion
    'D 1+'
    >>> positron.particle
    'e+'
    >>> electron.integer_charge
    -1
    >>> proton.spin
    0.5
    >>> alpha.atomic_number
    2
    >>> deuteron.mass_number
    2
    >>> deuteron.binding_energy.to('MeV')
    <Quantity 2.73557907 MeV>
    >>> alpha.charge
    <Quantity 3.20435324e-19 C>
    >>> neutron.half_life
    <Quantity 881.5 s>
    >>> Particle('C-14').half_life
    <Quantity 5730. yr>

    The `is_category` attribute may be used to determine whether or not
    a certain particle is or is not a member of a certain category.

    >>> electron.is_category('lepton')
    True
    >>> proton.is_category('baryon', exclude='charged')
    False
    >>> neutron.is_category({'matter', 'baryon'}, exclude={'charged'})
    True

    """

    def __init__(self,
                 argument: Union[str, int],
                 mass_numb: int = None,
                 Z: int = None):
        r"""Initializes a Particle object by setting all necessary private
        attributes."""

        if not isinstance(argument, (int, str)):
            raise TypeError(
                "The first positional argument when creating a Particle "
                "object must be either an integer or string.")

        if mass_numb is not None and not isinstance(mass_numb, int):
            raise TypeError("mass_numb is not an integer")

        if Z is not None and not isinstance(Z, int):
            raise TypeError("Z is not an integer.")

        self._attributes = collections.defaultdict(lambda: None)

        # Use this set to keep track of particle categories such as 'lepton'
        # for use with the is_category method later on.

        self._categories = set()

        # If the argument corresponds to one of the numerous case-sensitive or
        # case-insensitive aliases for particles, return the standard symbol.
        # Otherwise, just return the original argument.

        particle = _dealias_particle_aliases(argument)

        if particle in _Particles.keys():  # special particles

            special_particle_keys = [
                'particle',
                'name',
                'spin',
                'class',
                'lepton number',
                'baryon number',
                'integer charge',
                'half-life',
                'mass',
            ]

            for key in special_particle_keys:
                self._attributes[key] = _Particles[particle].get(key, None)

            particle_taxonomy_dict = ParticleZoo._taxonomy_dict
            categories = particle_taxonomy_dict.keys()

            for category in categories:
                if particle in particle_taxonomy_dict[category]:
                    self._categories.add(category)

            if particle == 'p+':

                # Protons are a special case amongst special cases, since they
                # are both a special particle and correspond to an element and
                # an isotope.  Protons are not as weird as electrons, though.
                # Electrons are weird.

                proton_keys_and_vals = [
                    ('atomic symbol', 'H'),
                    ('atomic number', 1),
                    ('mass number', 1),
                    ('element name', 'hydrogen'),
                    ('isotope', 'H-1'),
                    ('ion', 'H-1'),
                ]

                for key, val in proton_keys_and_vals:
                    self._attributes[key] = val

                self._categories.update({'element', 'isotope', 'ion'})

                if mass_numb is not None or Z is not None:
                    warnings.warn(
                        "Redundant mass number or charge information.",
                        AtomicWarning)

            elif mass_numb is not None or Z is not None:
                raise InvalidParticleError(
                    "The keywords 'mass_numb' and 'Z' cannot be used when "
                    "creating Particle objects for special particles. To "
                    f"create a Particle object for {self._name}s, "
                    f"use:  Particle({repr(self._attributes['particle'])})")

        else:  # elements, isotopes, and ions (besides protons)
            try:

                atomic_nomenclature_dict = _parse_and_check_atomic_input(
                    argument, mass_numb=mass_numb, Z=Z)

                element_keys = [
                    'particle',
                    'element',
                    'isotope',
                    'ion',
                    'mass number',
                    'Z',
                ]

                for key in element_keys:
                    self._attributes[key] = \
                        atomic_nomenclature_dict.get(key, None)

            except Exception as exc:
                errmsg = _invalid_particle_errmsg(
                    argument, mass_numb=mass_numb, Z=Z)
                raise InvalidParticleError(errmsg) from exc

            element = self._attributes['element']
            isotope = self._attributes['isotope']
            ion = self._attributes['ion']

            if element:
                self._categories.add('element')
            if isotope:
                self._categories.add('isotope')
            if ion:
                self._categories.add('ion')

            # Element properties

            self._attributes['atomic number'] = \
                _Elements[element]['atomic_number']
            self._attributes['element name'] = _Elements[element]['name']

            self._attributes['element name'] = _Elements[element]['name']
            self._attributes['baryon number'] = self._attributes['mass number']

            # For the moment, set the lepton number to zero for elements,
            # isotopes, and ions.  The lepton number will probably come up
            # primarily during

            self._attributes['lepton number'] = 0

            if isotope:
                if _Isotopes[isotope]['is_stable']:
                    self._attributes['half-life'] = np.inf * u.s
                else:
                    self._attributes['half-life'] = \
                        _Isotopes[isotope].get('half_life', None)
            elif element and not isotope:
                self._attributes['half-life'] = None

            if ion == 'He-4 2+':
                self._attributes['spin'] = 0
                self._categories.add('boson')

            # Set the masses

            if element and not isotope:
                try:
                    self._attributes['standard atomic weight'] = \
                        _Elements[element]['atomic_mass']
                except (KeyError, u.UnitConversionError):
                    self._attributes['standard atomic weight'] = None
            elif isotope:
                self._attributes['isotope mass'] = \
                    _Isotopes[isotope].get('atomic_mass', None)
            if element and not isotope and not ion:
                self._attributes['mass'] = \
                    self._attributes['standard atomic weight']
            elif isotope and not ion:
                self._attributes['mass'] = self._attributes['isotope mass']
            elif ion and isotope and self._attributes['isotope mass']:
                self._attributes['mass'] = \
                    self._attributes['isotope mass'] \
                    - self._attributes['integer charge'] * const.m_e
            elif ion and not isotope and \
                    self._attributes['standard atomic weight']:
                self._attributes['mass'] = \
                    self._attributes['standard atomic weight'] \
                    - self._attributes['integer charge'] * const.m_e
            else:
                self._attributes['mass'] = None

        # Set the charge

        if self._attributes['integer charge'] is not None:
            self._attributes['electric charge'] = \
                self._attributes['integer charge'] * const.e.si

        if self._attributes['integer charge']:
            self._categories.add('charged')

        self._element_errmsg = (
            f"The particle '{self.particle}' is not an element, so "
            f"this attribute is not available.")

        self._isotope_errmsg = (
            f"The particle '{self.particle}' does not have an "
            f"isotope specified, so this attribute is not available.")

        self._ion_errmsg = (
            f"The particle '{self.particle}' is not an ion, so this"
            f"attribute is not available.")

    def __repr__(self) -> str:
        r"""Returns a string of the call that would recreate this object."""
        return f'Particle("{self.particle}")'

    def __str__(self) -> str:
        r"""Returns a string of the particle symbol."""
        return f"{self.particle}"

    def __eq__(self, other) -> bool:
        r"""Returns True when comparing two Particle objects that correspond
        to the same particle, and False when the two objects differ."""
        try:
            if self.__dict__ == other.__dict__:
                return True
            else:
                return False
        except Exception:
            return False

    def __ne__(self, other) -> bool:
        r"""Returns True when the two objects differ, and False when
        comparing two Particle objects that correspond to the same particle."""
        return not self.__eq__(other)

    @property
    def particle(self) -> str:
        r"""Returns the particle symbol."""
        return self._attributes['particle']

    @property
    def element(self) -> Optional[str]:
        r"""Returns the atomic symbol if the particle corresponds to an
        element, and None otherwise."""
        return self._attributes['atomic_symbol']

    @property
    def isotope(self) -> Optional[str]:
        r"""Returns the isotope symbol if the particle corresponds to an
        isotope, and None otherwise."""
        return self._attributes['isotope']

    @property
    def ion(self) -> Optional[str]:
        r"""Returns the ion symbol if the particle corresponds to an ion,
        and None otherwise."""
        return self._attributes['ion']

    @property
    def element_name(self) -> str:
        r"""Returns the name of the element corresponding to this particle,
        or raises an InvalidElementError if the particle does not correspond
        to an element."""
        if not self.element:
            raise InvalidElementError(self._element_errmsg)
        return self._attributes['element name']

    @property
    def integer_charge(self) -> int:
        r"""Returns the integer charge of the particle, or raises a ChargeError
        if the charge has not been specified."""
        if self._attributes['integer charge'] is None:
            raise ChargeError(
                f"The charge of particle {self.particle} has not been "
                f"specified.")
        return self._attributes['integer charge']

    @property
    def charge(self) -> u.C:
        r"""Returns the electric charge as a Quantity in units of coulombs,
        or raises a ChargeError if the charge has not been specified."""
        if self._electric_charge is None:
            raise ChargeError(
                f"The charge of particle {self.particle} has not been "
                f"specified.")
        return self._electric_charge

    @property
    def mass(self) -> u.kg:
        r"""Returns the mass of the element, isotope, ion, particle, or
        antiparticle; or raises a MissingAtomicDataError if the mass
        is unavailable (e.g., if the particle is a neutrino).

        For special particles, this attribute will return the standard value of
        the mass of the particle.

        If the particle is an element and not an isotope or ion, then this
        attribute will return the standard atomic weight if available. If the
        particle is an isotope but not an ion, then this attribute will return
        the isotopic mass. If this particle is an ion, then this attribute will
        return the mass of the element or isotope (as just described) minus the
        integer charge times the electron mass."""
        if self._attributes['mass'] is None:
            raise MissingAtomicDataError(
                f"The mass of particle '{self.particle}' is unavailable.")
        return self._attributes['mass']

    @property
    def standard_atomic_weight(self) -> u.u:
        r"""Returns the standard atomic weight of an element if available.
        Raises a MissingAtomicDataError if the particle is an element for
        which the standard_atomic_weight is unavailable.  Raises an
        InvalidElementError if the particle is not an element."""
        if self.isotope or self.ion or not self.element:
            raise InvalidElementError(self._element_errmsg)
        if self._attributes['standard atomic weight'] is None:
            raise MissingAtomicDataError(
                f"The standard atomic weight of {self.element} is "
                f"unavailable.")
        return self._attributes['standard atomic weight']

    @property
    def nuclide_mass(self) -> u.kg:
        r"""Returns the mass of the nucleus of an isotope.  This attribute
        raises an InvalidIsotopeError if the particle is not an isotope or
        neutron, or a MissingAtomicDataError if the isotope mass is
        not available."""
        if self.particle in ['H-1', 'p+']:
            _nuclide_mass = const.m_p
        elif self.particle == 'n':
            _nuclide_mass = const.m_n
        elif not self.isotope:
            raise InvalidIsotopeError(self._isotope_errmsg)
        else:
            try:
                _nuclide_mass = self._isotope_mass \
                                - self._atomic_number * const.m_e
            except KeyError:  # coveralls: ignore
                raise MissingAtomicDataError(
                    f"The mass of a {self.isotope} nuclide is not available.")
        return _nuclide_mass

    @property
    def atomic_number(self) -> int:
        r"""Returns the atomic number of the element corresponding to this
        particle, or raises an InvalidElementError if the particle does not
        correspond to an element."""
        if not self.element:
            raise InvalidElementError(self._element_errmsg)
        return self._attributes['atomic number']

    @property
    def mass_number(self) -> int:
        r"""Returns the mass number of the isotope corresponding to this
        particle, or raises an InvalidIsotopeError if the particle does not
        correspond to an isotope."""
        if not self.isotope:
            raise InvalidIsotopeError(self._isotope_errmsg)
        return self._attributes['mass number']

    @property
    def baryon_number(self) -> int:
        r"""Returns the number of protons plus neutrons minus the number of
        antiprotons and antineutrons in the particle, or raises an
        AtomicError if the baryon number is unavailable.  The baryon number
        is equivalent to the mass number for isotopes."""
        if self._attributes['baryon number'] is None:  # coveralls: ignore
            raise AtomicError(
                f"The baryon number for '{self.particle}' is not "
                f"available.")
        return self._attributes['baryon number']

    @property
    def lepton_number(self) -> int:
        r"""Returns 1 for leptons, -1 for antileptons, and 0 for
        nuclides/isotopes; or raises an AtomicError if the lepton number is
        not available.  This attribute does not count the electrons in
        an atom or ion."""
        if self._attributes['lepton number'] is None:  # coveralls: ignore
            raise AtomicError(
                f"The lepton number for {self.particle} is not available.")
        return self._attributes['lepton number']

    @property
    def binding_energy(self) -> u.J:
        r"""Returns the nuclear binding energy, or raises an
        InvalidIsotopeError if the particle is not a nucleon or isotope."""

        if self._attributes['baryon number'] == 1:
            return 0 * u.J

        if not self.isotope:
            raise InvalidIsotopeError(
                f"The nuclear binding energy may only be calculated for "
                f"nucleons and isotopes.")

        number_of_protons = self.atomic_number
        number_of_neutrons = self.mass_number - self.atomic_number

        mass_of_protons = number_of_protons * const.m_p
        mass_of_neutrons = number_of_neutrons * const.m_n

        mass_of_nucleons = mass_of_protons + mass_of_neutrons
        mass_of_nuclide = self.mass - const.m_e * self.atomic_number

        mass_defect = mass_of_nucleons - mass_of_nuclide
        nuclear_binding_energy = mass_defect * const.c ** 2

        return nuclear_binding_energy.to(u.J)

    @property
    def half_life(self) -> u.s:
        r"""Returns the half-life of the particle, or raises a
        MissingAtomicDataError if the half-life is unavailable."""
        if self.element and not self.isotope:
            raise InvalidIsotopeError(self._isotope_errmsg)
        if self._attributes['half-life'] is None:
            raise MissingAtomicDataError(
                f"The half-life of '{self.particle}' is not available.")
        return self._attributes['half-life']

    @property
    def spin(self) -> Union[int, float]:
        r"""Returns the spin of the particle, or raises a
        MissingAtomicDataError if the spin is not available."""
        if self._attributes['spin'] is None:
            raise MissingAtomicDataError(
                f"The spin of particle '{self.particle}' is unavailable.")
        return self._attributes['spin']

    def reduced_mass(self, other, Z=None, mass_numb=None) -> u.kg:
        r"""Finds the reduced mass between two particles, or will raise a
        MissingAtomicDataError if either particle's mass is unavailable or
        an AtomicError for any other errors.  The other particle may be
        represented by another Particle object, a Quantity with units of,
        mass, or a string of the other particle's symbol (in conjunction
        with keywords Z and mass_numb)."""

        try:
            mass_this = self.mass.to(u.kg)
        except MissingAtomicDataError:
            raise MissingAtomicDataError(
                f"Unable to find the reduced mass because the mass of "
                f"{self.particle} is not available.") from None

        if isinstance(other, (str, int)):
                other = Particle(other, Z=Z, mass_numb=mass_numb)

        if isinstance(other, Particle):
            try:
                mass_that = other.mass.to(u.kg)
            except MissingAtomicDataError:
                raise MissingAtomicDataError(
                    f"Unable to find the reduced mass because the mass of "
                    f"{other.particle} is not available.") from None
        else:
            try:
                mass_that = other.to(u.kg)
            except Exception as exc:  # coveralls: ignore
                raise AtomicError(
                    f"{other} must be either a Particle or a Quantity or "
                    f"Constant with units of mass in order to calculate "
                    f"reduced mass.") from exc

        return (mass_this * mass_that) / (mass_this + mass_that)

    def is_category(self,
                    *categories,
                    any: bool = False,
                    exclude: Union[str, Set, Tuple, List] = set()
                    ) -> bool:
        r"""Returns True if the particle is in all of the inputted categories.

        If any is True, then the particle will return True if any of the listed
        categories are met.

        If the exclude keyword is set, then is_category will return False if
        the particle is in any of the excluded categories, whether or not the
        the particle matches the other criteria.

        The valid categories are: 'lepton', 'antilepton', 'baryon',
        'antibaryon', 'fermion', 'boson', 'neutrino', 'antineutrino', 'matter',
        'antimatter', 'element', 'isotope', 'ion', and 'charged'.

        Examples
        --------
        >>> neutron = Particle('n')
        >>> neutron.is_category({'fermion'})
        True
        >>> neutron.is_category('fermion', 'lepton', 'boson', any=True)
        True
        >>> neutron.is_category(['baryon', 'matter'], exclude=['fermion'])
        False
        >>> neutron.is_category([], exclude={'boson'})
        True

        """

        def _make_into_set(arg: Union[str, Set, Tuple, List]) -> Set[str]:
                r"""Turns the input (a string, set, tuple, or list) into
                a set containing the items in input."""
                if len(arg) == 0:
                    return set()

                if isinstance(arg, set):
                    return arg

                if isinstance(arg, str):
                    return {arg}

                if isinstance(arg[0], (tuple, list, set)):
                    return set(arg[0])
                else:
                    return set(arg)

        if not isinstance(any, bool):
            raise TypeError(
                f"The keyword any in {self.__repr__()}.is_category must be "
                f"set to either True or False.")
        elif any and len(categories) == 0:
            raise AtomicError(
                f"The keyword 'any' to {self.__repr__()}.is_category "
                f"cannot be set to True if no categories to be matched "
                f"are inputted.")

        categories = _make_into_set(categories)
        exclude = _make_into_set(exclude)

        # If valid_categories is changed, remember to change the docstring
        # for the Particle class.

        valid_categories = {
            'lepton', 'antilepton', 'fermion', 'boson', 'baryon', 'neutrino',
            'antineutrino', 'element', 'isotope', 'ion', 'matter',
            'antimatter', 'stable', 'unstable', 'charged',
        }

        if categories - valid_categories:
            raise AtomicError(
                f"The following categories in {self.__repr__()}.is_category "
                f"are not valid categories: {categories - valid_categories}")

        if exclude - valid_categories:
            raise AtomicError(
                f"The following categories to be excluded in "
                f"{self.__repr__()}.is_category are not valid categories: "
                f"{exclude - valid_categories}")

        if exclude & categories:
            raise AtomicError(
                f"The following are duplicate categories in "
                f"{self.__repr__()}.is_category: {categories & exclude}")

        if exclude & self._categories:
            return False

        if any and categories & self._categories:
            return True

        return categories <= self._categories
