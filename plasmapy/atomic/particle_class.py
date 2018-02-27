"""The Particle class."""

import numpy as np
import warnings
from typing import (Union, Set, Tuple, List, Optional)
import collections

import astropy.units as u
import astropy.constants as const

from ..utils import (
    AtomicError,
    AtomicWarning,
    InvalidParticleError,
    InvalidElementError,
    InvalidIsotopeError,
    ChargeError,
    MissingAtomicDataError,
    MissingAtomicDataWarning,
)

from .parsing import (
    _dealias_particle_aliases,
    _parse_and_check_atomic_input,
    _invalid_particle_errmsg,
)

from .elements import _Elements
from .isotopes import _Isotopes

from .special_particles import (_Particles, ParticleZoo, _special_ion_masses)

_PeriodicTable = collections.namedtuple(
    "periodic_table", ['group', 'category', 'block', 'period'])

_classification_categories = {
    'lepton',
    'antilepton',
    'fermion',
    'boson',
    'baryon',
    'neutrino',
    'antineutrino',
    'matter',
    'antimatter',
    'stable',
    'unstable',
    'charged',
    'uncharged',
}

_periodic_table_categories = {
    'nonmetal',
    'metal',
    'alkali metal',
    'alkaline earth metal',
    'metalloid',
    'transition metal',
    'post-transition metal',
    'halogen',
    'noble gas',
    'actinide',
    'lanthanide',
}

_atomic_property_categories = {
    'element',
    'isotope',
    'ion',
}

_specific_particle_categories = {
    'electron',
    'positron',
    'proton',
    'neutron',
}

_valid_categories = _periodic_table_categories | _classification_categories \
    | _atomic_property_categories | _specific_particle_categories


class Particle:
    """
    A class for an individual particle or antiparticle.

    Parameters
    ----------
    argument : `str` or `int`
        A string representing a particle, element, isotope, or ion; or
        an integer representing the atomic number of an element.

    mass_numb : `int`, optional
        The mass number of an isotope or nuclide.

    Z : `int`, optional
        The integer charge of the particle.

    Raises
    ------
    `~plasmapy.utils.InvalidParticleError`
        Raised when the particle input does not correspond to a valid
        particle or is contradictory.

    `~plasmapy.utils.InvalidElementError`
        For when an attribute is being accessed that requires
        information about an element, but the particle is not an
        element, isotope, or ion.

    `~plasmapy.utils.InvalidIsotopeError`
        For when an attribute is being accessed that requires
        information about an isotope or nuclide, but the particle is not
        an isotope (or an ion of an isotope).

    `~plasmapy.utils.ChargeError`
        For when either the charge or integer_charge attributes is being
        accessed but the charge information for the particle is not
        available.

    `TypeError`
        For when any of the arguments or keywords is not of the required
        type.

    Examples
    --------
    >>> proton = Particle('p+')
    >>> electron = Particle('e-')
    >>> neutron = Particle('neutron')
    >>> deuteron = Particle('D', Z=1)
    >>> alpha = Particle('He', mass_numb=4, Z=2)
    >>> positron = Particle('positron')
    >>> proton.element
    'H'
    >>> alpha.isotope
    'He-4'
    >>> deuteron.ion
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
    <Quantity 2.22456652 MeV>
    >>> alpha.charge
    <Quantity 3.20435324e-19 C>
    >>> neutron.half_life
    <Quantity 881.5 s>
    >>> Particle('C-14').half_life.to(u.year)
    <Quantity 5730. yr>

    The `periodic_table` attribute provides category, period, block, and
    group information if the particle corresponds to an element.

    >>> sulfur = Particle('S')
    >>> sulfur.periodic_table.category
    'nonmetal'
    >>> sulfur.periodic_table.block
    'p'
    >>> sulfur.periodic_table.period
    3
    >>> sulfur.periodic_table.group
    16

    The `is_category` attribute may be used to determine whether or not
    a certain particle is or is not a member of a certain category.

    >>> electron.is_category('lepton')
    True
    >>> proton.is_category('baryon', exclude='charged')
    False
    >>> neutron.is_category({'matter', 'baryon'}, exclude={'charged'})
    True
    >>> positron.is_category('antilepton', exclude='baryon')
    True

    Valid categories include `'lepton'`, `'antilepton'`, `'fermion'`,
    `'boson'`, `'baryon'`, `'neutrino'`, `'antineutrino'`, `'element'`,
    `'isotope'`, `'ion'`, `'matter'`, `'antimatter'`, `'nonmetal'`,
    `'metal'`, `'alkali metal'`, `'alkaline earth metal'`,
    `'metalloid'`, `'transition metal'`, `'post-transition metal',
    `'halogen'`, `'noble gas'`, `'actinide'`, `'lanthanide'`, 'stable'`,
    `'unstable'`, `'charged'`, `'uncharged'`, `'electron'`,
    `'positron'`, `'proton'`, and `'neutron'`.

    """

    def __init__(self,
                 argument: Union[str, int],
                 mass_numb: int = None,
                 Z: int = None):
        """Initialize a Particle object and set private attributes."""

        if not isinstance(argument, (int, str)):
            raise TypeError(
                "The first positional argument when creating a Particle "
                "object must be either an integer or string.")

        if mass_numb is not None and not isinstance(mass_numb, int):
            raise TypeError("mass_numb is not an integer")

        if Z is not None and not isinstance(Z, int):
            raise TypeError("Z is not an integer.")

        self._attributes = collections.defaultdict(lambda: None)
        attributes = self._attributes

        # Use this set to keep track of particle categories such as
        # 'lepton' for use with the is_category method later on.

        self._categories = set()

        # If the argument corresponds to one of the case-sensitive or
        # case-insensitive aliases for particles, return the standard
        # symbol. Otherwise, return the original argument.

        particle = _dealias_particle_aliases(argument)

        attributes['particle'] = particle

        if particle in _Particles.keys():  # special particles
            special_particle_keys = [
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
                attributes[key] = _Particles[particle].get(key, None)

            particle_taxonomy_dict = ParticleZoo._taxonomy_dict
            categories = particle_taxonomy_dict.keys()

            for category in categories:
                if particle in particle_taxonomy_dict[category]:
                    self._categories.add(category)

            if attributes['name'] in _specific_particle_categories:
                self._categories.add(attributes['name'])

            if particle == 'p+':

                # Protons are a special case amongst special cases, since they
                # are both a special particle and correspond to an element and
                # an isotope.  Protons are not as weird as electrons, though.
                # Electrons are weird.

                proton_keys_and_vals = [
                    ('element', 'H'),
                    ('atomic number', 1),
                    ('mass number', 1),
                    ('element name', 'hydrogen'),
                    ('isotope', 'H-1'),
                    ('ion', 'p+'),
                    ('mass', const.m_p),
                    ('integer charge', 1),
                ]

                self._periodic_table = _PeriodicTable(
                    group=_Elements['H']['group'],
                    period=_Elements['H']['period'],
                    block=_Elements['H']['block'],
                    category=_Elements['H']['category'],
                )

                for key, val in proton_keys_and_vals:
                    attributes[key] = val

                self._categories.update({'element', 'isotope', 'ion'})

            if mass_numb is not None or Z is not None:
                if particle == 'p+' and (mass_numb == 1 or Z == 1):
                    warnings.warn(
                        "Redundant mass number or charge information.",
                        AtomicWarning)
                else:
                    raise InvalidParticleError(
                        "The keywords 'mass_numb' and 'Z' cannot be used when "
                        "creating Particle objects for special particles. To "
                        f"create a Particle object for "
                        f"{attributes['name']}s, "
                        f"use:  Particle({repr(attributes['particle'])})")

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
                    'integer charge',
                ]

                for key in element_keys:
                    attributes[key] = \
                        atomic_nomenclature_dict.get(key, None)

            except Exception as exc:
                errmsg = _invalid_particle_errmsg(
                    argument, mass_numb=mass_numb, Z=Z)
                raise InvalidParticleError(errmsg) from exc

            element = attributes['element']
            isotope = attributes['isotope']
            ion = attributes['ion']

            if element:
                self._categories.add('element')
            if isotope:
                self._categories.add('isotope')
            if ion:
                self._categories.add('ion')

            # Element properties

            attributes['atomic number'] = \
                _Elements[element]['atomic_number']
            attributes['element name'] = _Elements[element]['name']

            attributes['element name'] = _Elements[element]['name']
            attributes['baryon number'] = attributes['mass number']

            # For the moment, set the lepton number to zero for elements,
            # isotopes, and ions.  The lepton number will probably come up
            # primarily during nuclear reactions.

            attributes['lepton number'] = 0

            if isotope:
                if _Isotopes[isotope]['stable']:
                    attributes['half-life'] = np.inf * u.s
                else:
                    attributes['half-life'] = \
                        _Isotopes[isotope].get('half-life', None)

            if ion == 'He-4 2+':
                attributes['spin'] = 0
                self._categories.add('boson')

            if element and not isotope:
                attributes['standard atomic weight'] = \
                    _Elements[element].get('atomic_mass', None)

            if isotope:
                attributes['isotope mass'] = \
                    _Isotopes[isotope].get('mass', None)
                attributes['isotopic abundance'] = \
                    _Isotopes[isotope].get('abundance', 0.0)

            if ion in _special_ion_masses.keys():
                attributes['mass'] = _special_ion_masses[ion]

            self._periodic_table = _PeriodicTable(
                group=_Elements[element]['group'],
                period=_Elements[element]['period'],
                block=_Elements[element]['block'],
                category=_Elements[element]['category'],
            )

            self._categories.add(self._periodic_table.category)

        if attributes['integer charge'] == 1:
            attributes['charge'] = const.e.si
        elif attributes['integer charge'] is not None:
            attributes['charge'] = \
                attributes['integer charge'] * const.e.si

        if attributes['integer charge']:
            self._categories.add('charged')
        elif attributes['integer charge'] == 0:
            self._categories.add('uncharged')

        if attributes['half-life'] is not None:
            if attributes['half-life'] == np.inf * u.s:
                self._categories.add('stable')
            else:
                self._categories.add('unstable')

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
        """Return a string of the call that would recreate this object."""
        return f'Particle("{self.particle}")'

    def __str__(self) -> str:
        """Return a string of the particle symbol."""
        return f"{self.particle}"

    def __eq__(self, other) -> bool:
        """
        Return True when comparing two Particle objects that correspond
        to the same particle, and False when the two objects differ.
        """

        if not isinstance(other, Particle):
            return False

        try:
            if self.particle == other.particle:
                return True
            else:
                return False
        except Exception:  # coveralls: ignore
            return False

    def __ne__(self, other) -> bool:
        """
        Return `True` when the two objects differ, and False when
        comparing two Particle objects that correspond to the same
        particle.
        """
        return not self.__eq__(other)

    @property
    def particle(self) -> str:
        """Return the particle symbol."""
        return self._attributes['particle']

    @property
    def element(self) -> Optional[str]:
        """
        Return the atomic symbol if the particle corresponds to an
        element, and `None` otherwise.
        """
        return self._attributes['element']

    @property
    def isotope(self) -> Optional[str]:
        """
        Return the isotope symbol if the particle corresponds to an
        isotope, and `None` otherwise.
        """
        return self._attributes['isotope']

    @property
    def ion(self) -> Optional[str]:
        """
        Return the ion symbol if the particle corresponds to an ion,
        and `None` otherwise.
        """
        return self._attributes['ion']

    @property
    def element_name(self) -> str:
        """
        Return the name of the element corresponding to this
        particle, or raises an `~plasmapy.utils.InvalidElementError` if
        the particle does not correspond to an element.
        """
        if not self.element:
            raise InvalidElementError(self._element_errmsg)
        return self._attributes['element name']

    @property
    def integer_charge(self) -> int:
        """
        Return the integer charge of the particle, or raises a
        `~plasmapy.utils.ChargeError` if the charge has not been
        specified.
        """
        if self._attributes['integer charge'] is None:
            raise ChargeError(
                f"The charge of particle {self.particle} has not been "
                f"specified.")
        return self._attributes['integer charge']

    @property
    def charge(self) -> u.Quantity:
        """
        Return the electric charge as a `~astropy.units.Quantity`
        in units of coulombs, or raises a `~plasmapy.utils.ChargeError`
        if the charge has not been specified.
        """
        if self._attributes['charge'] is None:
            raise ChargeError(
                f"The charge of particle {self.particle} has not been "
                f"specified.")
        if self._attributes['integer charge'] == 1:
            return const.e

        return self._attributes['charge']

    @property
    def standard_atomic_weight(self) -> u.Quantity:
        """
        Return the standard atomic weight of an element if
        available.  Raises a `~plasmapy.utils.MissingAtomicDataError` if
        the particle is an element for which the standard atomic weight
        is unavailable.  Raises an `~plasmapy.utils.InvalidElementError`
        if the particle is not an element.
        """
        if self.isotope or self.ion or not self.element:
            raise InvalidElementError(self._element_errmsg)
        if self._attributes['standard atomic weight'] is None:
            raise MissingAtomicDataError(
                f"The standard atomic weight of {self.element} is "
                f"unavailable.")
        return self._attributes['standard atomic weight']

    @property
    def nuclide_mass(self) -> u.Quantity:
        """
        Return the mass of the nucleus of an isotope.  This
        attribute raises an `~plasmapy.utils.InvalidIsotopeError` if the
        particle is not an isotope or neutron, or a
        `~plasmapy.utils.MissingAtomicDataError` if the isotope mass is
        not available.
        """

        if self.particle in ['H-1', 'p+']:
            return const.m_p
        elif self.particle == 'n':
            return const.m_n
        elif self.particle in ['D', 'D 1+']:
            return _special_ion_masses['D 1+']
        elif self.particle in ['T', 'T 1+']:
            return _special_ion_masses['T 1+']

        if not self.isotope:
            raise InvalidIsotopeError(self._isotope_errmsg)

        base_mass = self._attributes['isotope mass']

        if base_mass is None:
            raise MissingAtomicDataError(
                f"The mass of a {self.isotope} nuclide is not available.")

        _nuclide_mass = \
            self._attributes['isotope mass'] - self.atomic_number * const.m_e

        return _nuclide_mass.to(u.kg)

    @property
    def mass(self) -> u.Quantity:
        """
        Return the mass of the element, isotope, ion, particle, or
        antiparticle; or raises a
        `~plasmapy.utils.MissingAtomicDataError` if the mass is
        unavailable (e.g., if the particle is a neutrino).

        For special particles, this attribute will return the standard
        value of the mass of the particle.

        If the particle is an element and not an isotope or ion, then
        this attribute will return the standard atomic weight if
        available. If the particle is an isotope but not an ion, then
        this attribute will return the isotopic mass. If this particle
        is an ion, then this attribute will return the mass of the
        element or isotope (as just described) minus the integer charge
        times the electron mass.

        """

        if self._attributes['mass'] is not None:
            return self._attributes['mass']

        if self.ion:

            if self.isotope:
                base_mass = self._attributes['isotope mass']
            else:
                base_mass = self._attributes['standard atomic weight']

            if base_mass is None:
                raise MissingAtomicDataError(
                    f"The mass of ion '{self.ion}' is not available.")

            mass = base_mass - self.integer_charge * const.m_e

            return mass.to(u.kg)

        if self.element:

            if self.isotope:
                mass = self._attributes['isotope mass']
            else:
                mass = self._attributes['standard atomic weight']

            if mass is not None:
                return mass

        raise MissingAtomicDataError(
            f"The mass of {self.particle} is not available.")

    @property
    def atomic_number(self) -> int:
        """
        Return the atomic number of the element corresponding to
        this particle, or raises an `~plasmapy.utils.InvalidElementError`
        if the particle does not correspond to an element.
        """
        if not self.element:
            raise InvalidElementError(self._element_errmsg)
        return self._attributes['atomic number']

    @property
    def mass_number(self) -> int:
        """
        Return the mass number of the isotope corresponding to this
        particle, or raises an `~plasmapy.utils.InvalidIsotopeError` if
        the particle does not correspond to an isotope.
        """
        if not self.isotope:
            raise InvalidIsotopeError(self._isotope_errmsg)
        return self._attributes['mass number']

    @property
    def neutron_number(self) -> int:
        """
        Return the number of neutrons of the isotope corresponding
        to this particle, or raises an `~plasmapy.utils.InvalidIsotopeError` if
        the particle does not correspond to an isotope.
        """
        return self.mass_number - self.atomic_number

    @property
    def isotopic_abundance(self) -> u.Quantity:
        """
        Return the isotopic abundance of an isotope if available,
        or raises a `~plasmapy.utils.MissingAtomicDataError`.
        """
        if not self.isotope:
            raise InvalidIsotopeError(self._isotope_errmsg)
        # TODO: Raise exception when element doesn't occur naturally
        return self._attributes.get('isotopic abundance', 0.0)

    @property
    def baryon_number(self) -> int:
        """
        Return the number of protons plus neutrons minus the number
        of antiprotons and antineutrons in the particle, or raises an
        `~plasmapy.utils.AtomicError` if the baryon number is
        unavailable. The baryon number is equivalent to the mass number
        for isotopes.
        """
        if self._attributes['baryon number'] is None:  # coveralls: ignore
            raise AtomicError(
                f"The baryon number for '{self.particle}' is not "
                f"available.")
        return self._attributes['baryon number']

    @property
    def lepton_number(self) -> int:
        """
        Return 1 for leptons, -1 for antileptons, and 0 for
        nuclides/isotopes; or raises an `~plasmapy.utils.AtomicError` if
        the lepton number is not available.  This attribute does not
        count the electrons in an atom or ion.
        """
        if self._attributes['lepton number'] is None:  # coveralls: ignore
            raise AtomicError(
                f"The lepton number for {self.particle} is not available.")
        return self._attributes['lepton number']

    @property
    def binding_energy(self) -> u.Quantity:
        """
        Return the nuclear binding energy, or raises an
        `~plasmapy.utils.InvalidIsotopeError` if the particle is not a
        nucleon or isotope.
        """

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

        mass_defect = mass_of_nucleons - self.nuclide_mass
        nuclear_binding_energy = mass_defect * const.c ** 2

        return nuclear_binding_energy.to(u.J)

    @property
    def half_life(self) -> u.Quantity:
        """
        Return the half-life of the particle, or raises a
        `~plasmapy.utils.MissingAtomicDataError` if the half-life is
        unavailable.
        """
        if self.element and not self.isotope:
            raise InvalidIsotopeError(self._isotope_errmsg)

        if isinstance(self._attributes['half-life'], str):
            warnings.warn(
                f"The half-life for {self.particle} is not known precisely; "
                "returning string with estimated value.",
                MissingAtomicDataWarning)

        if self._attributes['half-life'] is None:
            raise MissingAtomicDataError(
                f"The half-life of '{self.particle}' is not available.")
        return self._attributes['half-life']

    @property
    def spin(self) -> Union[int, float]:
        """
        Return the spin of the particle, or raises a
        `~plasmapy.utils.MissingAtomicDataError` if the spin is not
        available.
        """
        if self._attributes['spin'] is None:
            raise MissingAtomicDataError(
                f"The spin of particle '{self.particle}' is unavailable.")
        return self._attributes['spin']

    @property
    def periodic_table(self):
        """
        Return a `~collections.namedtuple` to access category, period,
        group, and block information about an element.

        If the particle is not an element, isotope, or ion, then this
        method will raise an `~plasmapy.utils.InvalidElementError`.
        """
        if self.element:
            return self._periodic_table
        else:
            raise InvalidElementError(self._element_errmsg)

    def reduced_mass(self, other, Z=None, mass_numb=None) -> u.Quantity:
        """
        Find the reduced mass between two particles, or will raise
        a `~plasmapy.utils.MissingAtomicDataError` if either particle's
        mass is unavailable or an `~plasmapy.utils.AtomicError` for any
        other errors.  The other particle may be represented by another
        Particle object, a `~astropy.units.Quantity` with units of mass,
        or a string of the other particle's symbol (in conjunction with
        keywords `Z` and `mass_numb`).

        """

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

    @property
    def categories(self) -> Set[str]:
        """Return the particle's categories."""
        return self._categories

    def is_category(self,
                    *category_tuple,
                    require: Union[str, Set, Tuple, List] = set(),
                    any_of: Union[str, Set, Tuple, List] = set(),
                    exclude: Union[str, Set, Tuple, List] = set(),
                    ) -> bool:
        """
        Determine if the particle meets categorization criteria.

        Return `True` if the particle is in all of the inputted
        categories, and `False` the particle is not.

        Required categories may be entered as positional arguments,
        including as a `list`, `set`, or `tuple` of required categories.
        These may also be included using the `require` keyword argument.
        This method will return `False` if the particle is not in all of
        the required categories.

        If categories are inputted using the `any_of` keyword argument,
        then this method will return `False` if the particle is not of
        any of the categories in `any_of`.

        If the `exclude` keyword is set, then this method will return
        `False` if the particle is in any of the excluded categories,
        whether or not the particle matches the other criteria.

        """

        def become_set(arg: Union[str, Set, Tuple, List]) -> Set[str]:
                """Change the argument into a `set`."""
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

        if category_tuple != () and require != set():
            raise AtomicError(
                "No positional arguments are allowed if the `require` keyword "
                "is set in is_category.")

        require = become_set(category_tuple) if category_tuple else \
            become_set(require)

        exclude = become_set(exclude)
        any_of = become_set(any_of)

        if not require and not exclude and not any_of:
            return _valid_categories

        invalid_categories = (require | exclude | any_of) - _valid_categories

        duplicate_categories = \
            require & exclude | exclude & any_of | require & any_of

        categories_and_adjectives = [
            (invalid_categories, 'invalid'),
            (duplicate_categories, 'duplicated')]

        for problem_categories, adjective in categories_and_adjectives:
            if problem_categories:
                raise AtomicError(
                    f"The following categories in {self.__repr__()}"
                    f".is_category are {adjective}: {problem_categories}")

        if exclude and exclude & self._categories:
            return False

        if any_of and not any_of & self._categories:
            return False

        return require <= self._categories

    def is_electron(self) -> bool:
        """
        Returns True if the particle is an electron.
        """
        return self.particle == "e-"
