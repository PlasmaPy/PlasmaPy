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
    InvalidIonError,
    ChargeError,
    MissingAtomicDataError,
    MissingAtomicDataWarning,
)

from .parsing import (
    _dealias_particle_aliases,
    _parse_and_check_atomic_input,
    _invalid_particle_errmsg,
)

from .elements import _Elements, _PeriodicTable
from .isotopes import _Isotopes

from .special_particles import (
    _Particles,
    ParticleZoo,
    _special_ion_masses,
    _antiparticles,
)

__all__ = [
    "Particle",
]

_classification_categories = {
    'lepton',
    'antilepton',
    'fermion',
    'boson',
    'antibaryon',
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

_valid_categories = (
    _periodic_table_categories
    | _classification_categories
    | _atomic_property_categories
    | _specific_particle_categories
)


def _category_errmsg(particle, category: str) -> str:
    """
    Return an error message when an attribute raises an
    `~plasmapy.utils.InvalidElementError`,
    `~plasmapy.utils.InvalidIonError`, or
    `~plasmapy.utils.InvalidIsotopeError`.
    """
    article = 'an' if category[0] in 'aeiouAEIOU' else 'a'
    errmsg = (
        f"The particle {particle} is not {article} {category}, "
        f"so this attribute is not available.")
    return errmsg


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
    `TypeError`
        For when any of the arguments or keywords is not of the required
        type.

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
        For when either the `~plasmapy.atomic.Particle.charge` or
        `~plasmapy.atomic.Particle.integer_charge` attributes is
        being accessed but the charge information for the particle is
        not available.

    `~plasmapy.utils.AtomicError`
        Raised for attempts at converting a
        `~plasmapy.atomic.Particle` object to a `bool`.

    Examples
    --------
    Particles may be defined using a wide variety of aliases:

    >>> proton = Particle('p+')
    >>> electron = Particle('e-')
    >>> neutron = Particle('neutron')
    >>> deuteron = Particle('D', Z=1)
    >>> alpha = Particle('He', mass_numb=4, Z=2)
    >>> positron = Particle('positron')
    >>> hydrogen = Particle(1)  # atomic number

    The `particle` attribute returns the particle's symbol in the
    standard form.

    >>> positron.particle
    'e+'

    The `atomic_symbol`, `isotope_symbol`, and `ionic_symbol` attributes
    return the symbols for each of these different types of particles.

    >>> proton.element
    'H'
    >>> alpha.isotope
    'He-4'
    >>> deuteron.ionic_symbol
    'D 1+'

    The `ionic_symbol` attribute works for neutral atoms if charge
    information is available.

    >>> deuterium = Particle("D", Z=0)
    >>> deuterium.ionic_symbol
    'D 0+'

    If the particle doesn't belong to one of those categories, then
    these attributes return `None`.

    >>> positron.element is None
    True

    The attributes of a `~plasmapy.atomic.Particle` instance may be used
    to test whether or not a particle is an element, isotope, or ion.

    >>> True if positron.element else False
    False
    >>> True if deuterium.isotope else False
    True
    >>> True if Particle('alpha').is_ion else False
    True

    Many of the attributes return physical properties of a particle.

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
    >>> deuteron.electron_number
    0
    >>> alpha.neutron_number
    2

    If a `~plasmapy.atomic.Particle` instance represents an elementary
    particle, then the unary `~` (invert) operator may be used to return
    the particle's antiparticle.

    >>> ~electron
    Particle("e+")
    >>> ~proton
    Particle("p-")
    >>> ~positron
    Particle("e-")

    The `~plasmapy.atomic.particle_class.Particle.categories` attribute
    and `~plasmapy.atomic.particle_class.Particle.is_category` method
    may be used to find and test particle membership in categories.

    Valid particle categories include: `'actinide'`, `'alkali
    metal'`, `'alkaline earth metal'`, `'antibaryon'`,
    `'antilepton'`, `'antimatter'`, `'antineutrino'`, `'baryon'`,
    `'boson'`, `'charged'`, `'electron'`, `'element'`,
    `'fermion'`, `'halogen'`, `'ion'`, `'isotope'`,
    `'lanthanide'`, `'lepton'`, `'matter'`, `'metal'`,
    `'metalloid'`, `'neutrino'`, `'neutron'`, `'noble gas'`,
    `'nonmetal'`, `'positron'`, `'post-transition metal'`,
    `'proton'`, `'stable'`, `'transition metal'`, `'uncharged'`,
    and `'unstable'`.

"""

    def __init__(self, argument: Union[str, int], mass_numb: int = None, Z: int = None):
        """
        Initialize a `~plasmapy.atomic.Particle` object and set private
        attributes.
        """

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
        categories = self._categories

        # If the argument corresponds to one of the case-sensitive or
        # case-insensitive aliases for particles, return the standard
        # symbol. Otherwise, return the original argument.

        particle = _dealias_particle_aliases(argument)

        if particle in _Particles.keys():  # special particles

            attributes['particle'] = particle

            for attribute in _Particles[particle].keys():
                attributes[attribute] = _Particles[particle][attribute]

            particle_taxonomy = ParticleZoo._taxonomy_dict
            all_categories = particle_taxonomy.keys()

            for category in all_categories:
                if particle in particle_taxonomy[category]:
                    categories.add(category)

            if attributes['name'] in _specific_particle_categories:
                categories.add(attributes['name'])

            if particle == 'p+':
                categories.update({'element', 'isotope', 'ion'})

            if mass_numb is not None or Z is not None:
                if particle == 'p+' and (mass_numb == 1 or Z == 1):
                    warnings.warn("Redundant mass number or charge information.", AtomicWarning)
                else:
                    raise InvalidParticleError(
                        "The keywords 'mass_numb' and 'Z' cannot be used when "
                        "creating Particle objects for special particles. To "
                        f"create a Particle object for {attributes['name']}s, "
                        f"use:  Particle({repr(attributes['particle'])})")

        else:  # elements, isotopes, and ions (besides protons)
            try:
                nomenclature = _parse_and_check_atomic_input(argument, mass_numb=mass_numb, Z=Z)
            except Exception as exc:
                errmsg = _invalid_particle_errmsg(argument, mass_numb=mass_numb, Z=Z)
                raise InvalidParticleError(errmsg) from exc

            for key in nomenclature.keys():
                attributes[key] = nomenclature[key]

            element = attributes['element']
            isotope = attributes['isotope']
            ion = attributes['ion']

            if element:
                categories.add('element')
            if isotope:
                categories.add('isotope')
            if self.element and self._attributes['integer charge']:
                categories.add('ion')

            # Element properties

            Element = _Elements[element]

            attributes['atomic number'] = Element['atomic number']
            attributes['element name'] = Element['element name']

            # Set the lepton number to zero for elements, isotopes, and
            # ions.  The lepton number will probably come up primarily
            # during nuclear reactions.

            attributes['lepton number'] = 0

            if isotope:

                Isotope = _Isotopes[isotope]

                attributes['baryon number'] = Isotope['mass number']
                attributes['isotope mass'] = Isotope.get('mass', None)
                attributes['isotopic abundance'] = Isotope.get('abundance', 0.0)

                if Isotope['stable']:
                    attributes['half-life'] = np.inf * u.s
                else:
                    attributes['half-life'] = Isotope.get('half-life', None)

            if element and not isotope:
                attributes['standard atomic weight'] = Element.get('atomic mass', None)

            if ion in _special_ion_masses.keys():
                attributes['mass'] = _special_ion_masses[ion]

            attributes['periodic table'] = _PeriodicTable(
                group=Element['group'],
                period=Element['period'],
                block=Element['block'],
                category=Element['category'],
            )

            categories.add(Element['category'])

        if attributes['integer charge'] == 1:
            attributes['charge'] = const.e.si
        elif attributes['integer charge'] is not None:
            attributes['charge'] = attributes['integer charge'] * const.e.si

        if attributes['integer charge']:
            categories.add('charged')
        elif attributes['integer charge'] == 0:
            categories.add('uncharged')

        if attributes['half-life'] is not None:
            if attributes['half-life'] == np.inf * u.s:
                categories.add('stable')
            else:
                categories.add('unstable')

        self.__name__ = self.__repr__()

    def __repr__(self) -> str:
        """Return a call string that would recreate this object.

        Examples
        --------
        >>> lead = Particle('lead')
        >>> repr(lead)
        'Particle("Pb")'

        """
        return f'Particle("{self.particle}")'

    def __str__(self) -> str:
        """Return the particle's symbol."""
        return self.particle

    def __eq__(self, other) -> bool:
        """
        Determine if two objects correspond to the same particle.

        This method will return `True` if `other` is an identical
        `~plasmapy.atomic.Particle` instance or a `str` representing the
        same particle, and return `False` if `other` is a different
        `~plasmapy.atomic.Particle` or a `str` representing a different
        particle.

        If `other` is not a `str` or `~plasmapy.atomic.Particle`
        instance, then this method will raise a `TypeError`.  If
        `other.particle` equals `self.particle` but the attributes
        differ, then this method will raise a
        `~plasmapy.utils.AtomicError`.

        Examples
        --------
        >>> electron = Particle('e-')
        >>> positron = Particle('e+')
        >>> electron == positron
        False
        >>> electron == 'e-'
        True

        """
        if isinstance(other, str):
            try:
                other_particle = Particle(other)
                return self.particle == other_particle.particle
            except InvalidParticleError as exc:
                raise InvalidParticleError(
                    f"{other} is not a particle and cannot be "
                    f"compared to {self}.") from exc

        if not isinstance(other, self.__class__):
            raise TypeError(
                f"The equality of a Particle object with a {type(other)} is undefined.")

        no_particle_attr = 'particle' not in dir(self) or 'particle' not in dir(other)
        no_attributes_attr = '_attributes' not in dir(self) or '_attributes' not in dir(other)

        if no_particle_attr or no_attributes_attr:  # coveralls: ignore
            raise TypeError(f"The equality of {self} with {other} is undefined.")

        same_particle = self.particle == other.particle

        # The following two loops are a hack to enable comparisons
        # between defaultdicts.  By accessing all of the defined keys in
        # each of the defaultdicts, this makes sure that
        # self._attributes and other._attributes have the same keys.

        # TODO: create function in utils to account for equality between
        # defaultdicts, and implement it here

        for attribute in self._attributes.keys():
            other._attributes[attribute]

        for attribute in other._attributes.keys():
            self._attributes[attribute]

        same_attributes = self._attributes == other._attributes

        if same_particle and not same_attributes:  # coveralls: ignore
            raise AtomicError(
                f"{self} and {other} should be the same Particle, but "
                f"have differing attributes.\n\n"
                f"The attributes of {self} are:\n\n{self._attributes}\n\n"
                f"The attributes of {other} are:\n\n{other._attributes}\n")

        return same_particle

    def __ne__(self, other) -> bool:
        """
        Test whether or not two objects are different particles.

        This method will return `False` if `other` is an identical
        `~plasmapy.atomic.Particle` instance or a `str` representing the
        same particle, and return `True` if `other` is a different
        `~plasmapy.atomic.Particle` or a `str` representing a different
        particle.

        If `other` is not a `str` or `~plasmapy.atomic.Particle`
        instance, then this method will raise a `TypeError`.  If
        `other.particle` equals `self.particle` but the attributes
        differ, then this method will raise a
        `~plasmapy.utils.AtomicError`.

        """
        return not self.__eq__(other)

    def __bool__(self):
        """
        Raise an `~plasmapy.utils.AtomicError` because Particle objects
        do not have a truth value.
        """
        raise AtomicError("The truthiness of a Particle object is not defined.")

    def __invert__(self):
        """
        Return the corresponding antiparticle, or raise an
        `~plasmapy.utils.AtomicError` if the particle is not an
        elementary particle.
        """
        return self.antiparticle

    @property
    def particle(self) -> Optional[str]:
        """
        Return the particle's symbol.

        Examples
        --------
        >>> electron = Particle('electron')
        >>> electron.particle
        'e-'

        """
        return self._attributes['particle']

    @property
    def antiparticle(self):
        """
        Return the corresponding antiparticle, or raise an
        `~plasmapy.utils.AtomicError` if the particle is not an
        elementary particle.

        This attribute may be accessed by using the unary operator `~`
        acting on a `~plasma.atomic.Particle` instance.

        Examples
        --------
        >>> electron = Particle('e-')
        >>> electron.antiparticle
        Particle("e+")

        >>> antineutron = Particle('antineutron')
        >>> ~antineutron
        Particle("n")

        """
        if self.particle in _antiparticles.keys():
            return Particle(_antiparticles[self.particle])
        else:
            raise AtomicError(
                "The unary operator can only be used for elementary "
                "particles and antiparticles.")

    @property
    def element(self) -> Optional[str]:
        """
        Return the atomic symbol if the particle corresponds to an
        element, and `None` otherwise.

        Examples
        --------
        >>> alpha = Particle('alpha')
        >>> alpha.element
        'He'

        """
        return self._attributes['element']

    @property
    def isotope(self) -> Optional[str]:
        """
        Return the isotope symbol if the particle corresponds to an
        isotope, and `None` otherwise.

        Examples
        --------
        >>> alpha = Particle('alpha')
        >>> alpha.isotope
        'He-4'

        """
        return self._attributes['isotope']

    @property
    def ionic_symbol(self) -> Optional[str]:
        """
        Return the ionic symbol if the particle corresponds to an ion or
        neutral atom, and `None` otherwise.

        Examples
        --------
        >>> deuteron = Particle('deuteron')
        >>> deuteron.ionic_symbol
        'D 1+'
        >>> hydrogen_atom = Particle('H', Z=0)
        >>> hydrogen_atom.ionic_symbol
        'H 0+'

        """
        return self._attributes['ion']

    @property
    def element_name(self) -> str:
        """
        Return the name of the element corresponding to this
        particle, or raise an `~plasmapy.utils.InvalidElementError` if
        the particle does not correspond to an element.

        Examples
        --------
        >>> tritium = Particle('T')
        >>> tritium.element_name
        'hydrogen'

        """
        if not self.element:
            raise InvalidElementError(_category_errmsg(self, 'element'))
        return self._attributes['element name']

    @property
    def integer_charge(self) -> int:
        """
        Return the particle's integer charge.

        This attribute will raise a `~plasmapy.utils.ChargeError` if the
        charge has not been specified.

        Examples
        --------
        >>> muon = Particle('mu-')
        >>> muon.integer_charge
        -1

        """
        if self._attributes['integer charge'] is None:
            raise ChargeError(f"The charge of particle {self} has not been specified.")
        return self._attributes['integer charge']

    @property
    def charge(self) -> u.Quantity:
        """
        Return the particle's electron charge in coulombs.

        This attribute will raise a `~plasmapy.utils.ChargeError` if the
        charge has not been specified.

        Examples
        --------
        >>> electron = Particle('e-')
        >>> electron.charge
        <Quantity -1.60217662e-19 C>

        """
        if self._attributes['charge'] is None:
            raise ChargeError(f"The charge of particle {self} has not been specified.")
        if self._attributes['integer charge'] == 1:
            return const.e.si

        return self._attributes['charge']

    @property
    def standard_atomic_weight(self) -> u.Quantity:
        """
        Return an element's standard atomic weight in kg.

        If the particle is isotope or ion or not an element, this
        attribute will raise an `~plasmapy.utils.InvalidElementError`.

        If the element does not have a defined stsandard atomic weight,
        this attribute will raise a
        `~plasmapy.utils.MissingAtomicDataError`.

        Examples
        --------
        >>> oxygen = Particle('O')
        >>> oxygen.standard_atomic_weight
        <Quantity 2.65669641e-26 kg>

        """
        if self.isotope or self.is_ion or not self.element:
            raise InvalidElementError(_category_errmsg(self, 'element'))
        if self._attributes['standard atomic weight'] is None:  # coveralls: ignore
            raise MissingAtomicDataError(
                f"The standard atomic weight of {self} is unavailable.")
        return self._attributes['standard atomic weight'].to(u.kg)

    @property
    def nuclide_mass(self) -> u.Quantity:
        """
        Return the mass of the bare nucleus of an isotope or a neutron.

        This attribute will raise a
        `~plasmapy.utils.InvalidIsotopeError` if the particle is not an
        isotope or neutron, or a
        `~plasmapy.utils.MissingAtomicDataError` if the isotope mass is
        not available.

        Examples
        --------
        >>> deuterium = Particle('D')
        >>> deuterium.nuclide_mass
        <Quantity 3.34358372e-27 kg>

        """

        if self.isotope == 'H-1':
            return const.m_p
        elif self.isotope == 'D':
            return _special_ion_masses['D 1+']
        elif self.isotope == 'T':
            return _special_ion_masses['T 1+']
        elif self.particle == 'n':
            return const.m_n

        if not self.isotope:
            raise InvalidIsotopeError(_category_errmsg(self, 'isotope'))

        base_mass = self._attributes['isotope mass']

        if base_mass is None:  # coveralls: ignore
            raise MissingAtomicDataError(f"The mass of a {self.isotope} nuclide is not available.")

        _nuclide_mass = self._attributes['isotope mass'] - self.atomic_number * const.m_e

        return _nuclide_mass.to(u.kg)

    @property
    def mass(self) -> u.Quantity:
        """
        Return the mass of the particle in kilograms.

        If the particle is an element and not an isotope or ion, then
        this attribute will return the standard atomic weight, if
        available.

        If the particle is an isotope but not an ion, then this
        attribute will return the isotopic mass, including bound
        electrons.

        If the particle is an ion, then this attribute will return the
        mass of the element or isotope (as just described) minus the
        product of the integer charge and the electron mass.

        For special particles, this attribute will return the standard
        value for the particle's mass.

        If the mass is unavailable (e.g., for neutrinos or elements with
        no standard atomic weight), then this attribute will raise a
        `~plasmapy.utils.MissingAtomicDataError`.

        Examples
        --------
        >>> Particle('He').mass
        <Quantity 6.64647688e-27 kg>
        >>> Particle('He+').mass
        <Quantity 6.64556594e-27 kg>
        >>> Particle('He-4 +1').mass
        <Quantity 6.64556803e-27 kg>
        >>> Particle('alpha').mass
        <Quantity 6.64465709e-27 kg>

        """

        if self._attributes['mass'] is not None:
            return self._attributes['mass'].to(u.kg)

        if self.is_ion:

            if self.isotope:
                base_mass = self._attributes['isotope mass']
            else:
                base_mass = self._attributes['standard atomic weight']

            if base_mass is None:
                raise MissingAtomicDataError(
                    f"The mass of ion '{self.ionic_symbol}' is not available."
                )

            mass = base_mass - self.integer_charge * const.m_e

            return mass.to(u.kg)

        if self.element:

            if self.isotope:
                mass = self._attributes['isotope mass']
            else:
                mass = self._attributes['standard atomic weight']

            if mass is not None:
                return mass.to(u.kg)

        raise MissingAtomicDataError(f"The mass of {self} is not available.")

    @property
    def atomic_number(self) -> int:
        """
        Return the number of protons in an element, isotope, or ion.

        If the particle is not an element, then this attribute will
        raise an `~plasmapy.utils.InvalidElementError`.

        Examples
        --------
        >>> proton = Particle('p+')
        >>> proton.atomic_number
        1
        >>> curium = Particle('Cm')
        >>> curium.atomic_number
        96

        """
        if not self.element:
            raise InvalidElementError(_category_errmsg(self, 'element'))
        return self._attributes['atomic number']

    @property
    def mass_number(self) -> int:
        """
        Return the number of nucleons in an isotope.

        This attribute will return the number of protons plus the number
        of neutrons in an isotope or nuclide.

        If the particle is not an isotope, then this attribute will
        raise an `~plasmapy.utils.InvalidIsotopeError`.

        Examples
        --------
        >>> alpha = Particle('helium-4 2+')
        >>> alpha.mass_number
        4

        """
        if not self.isotope:
            raise InvalidIsotopeError(_category_errmsg(self, 'isotope'))
        return self._attributes['mass number']

    @property
    def neutron_number(self) -> int:
        """
        Return the number of neutrons in an isotope or nucleon.

        This attribute will return the number of neutrons in an isotope,
        or `1` for a neutron.

        If this particle is not an isotope or neutron, then this
        attribute will raise an `~plasmapy.utils.InvalidIsotopeError`.

        Examples
        --------
        >>> alpha = Particle('He-4++')
        >>> alpha.neutron_number
        2
        >>> Particle('n').neutron_number
        1

        """
        if self.particle == 'n':
            return 1
        elif self.isotope:
            return self.mass_number - self.atomic_number
        else:  # coveralls: ignore
            raise InvalidIsotopeError(_category_errmsg(self, 'isotope'))

    @property
    def electron_number(self) -> int:
        """
        Return the number of electrons in an ion.

        This attribute will return the number of bound electrons in an
        ion, or `1` for an electron.

        If this particle is not an ion or electron, then this attribute
        will raise an `~plasmapy.utils.InvalidIonError`.

        Examples
        --------
        >>> Particle('Li 0+').electron_number
        3
        >>> Particle('e-').electron_number
        1

        """
        if self.particle == 'e-':
            return 1
        elif self.ionic_symbol:
            return self.atomic_number - self.integer_charge
        else:  # coveralls: ignore
            raise InvalidIonError(_category_errmsg(self, 'ion'))

    @property
    def isotopic_abundance(self) -> u.Quantity:
        """
        Return the isotopic abundance of an isotope.

        If the isotopic abundance is not available, this attribute will
        raise a `~plasmapy.utils.MissingAtomicDataError`.  If the
        particle is not an isotope or is an ion of an isotope, then this
        attribute will raise an `~plasmapy.utils.InvalidIsotopeError`.

        Examples
        --------
        >>> D = Particle('deuterium')
        >>> D.isotopic_abundance
        0.000115

        """
        from .atomic import common_isotopes

        if not self.isotope or self.is_ion:  # coveralls: ignore
            raise InvalidIsotopeError(_category_errmsg(self.particle, 'isotope'))

        abundance = self._attributes.get('isotopic abundance', 0.0)

        if not common_isotopes(self.element):
            warnings.warn(
                f'No isotopes of {self.element} have an isotopic abundance. '
                f'The isotopic abundance of {self.isotope} is being returned as 0.0',
                AtomicWarning)

        return abundance

    @property
    def baryon_number(self) -> int:
        """
        Return the number of baryons in a particle.

        This attribute will return the number of protons and neutrons
        minus the number of antiprotons and antineutrons. The baryon
        number is equivalent to the mass number for isotopes.

        If the baryon number is unavailable, then this attribute will
        raise a `~plasmapy.utils.MissingAtomicDataError`.

        Examples
        --------
        >>> alpha = Particle('alpha')
        >>> alpha.baryon_number
        4

        """
        if self._attributes['baryon number'] is None:  # coveralls: ignore
            raise MissingAtomicDataError(
                f"The baryon number for '{self.particle}' is not available.")
        return self._attributes['baryon number']

    @property
    def lepton_number(self) -> int:
        """
        Return `1` for leptons, `-1` for antileptons, and `0` otherwise.

        This attribute returns the number of leptons minus the number of
        antileptons, excluding bound electrons in an atom or ion.

        If the lepton number is unavailable, then this attribute will
        raise a `~plasmapy.utils.MissingAtomicDataError`.

        Examples
        --------
        >>> Particle('e-').lepton_number
        1
        >>> Particle('mu+').lepton_number
        -1
        >>> Particle('He-4 0+').lepton_number
        0

        """
        if self._attributes['lepton number'] is None:  # coveralls: ignore
            raise MissingAtomicDataError(
                f"The lepton number for {self.particle} is not available.")
        return self._attributes['lepton number']

    @property
    def binding_energy(self) -> u.Quantity:
        """
        Return the nuclear binding energy in joules.

        This attribute will raise an
        `~plasmapy.utils.InvalidIsotopeError` if the particle is not a
        nucleon or isotope.

        Examples
        --------
        >>> alpha = Particle('alpha')
        >>> alpha.binding_energy
        <Quantity 4.53346938e-12 J>
        >>> Particle('T').binding_energy.to('MeV')
        <Quantity 8.48179621 MeV>

        """

        if self._attributes['baryon number'] == 1:
            return 0 * u.J

        if not self.isotope:
            raise InvalidIsotopeError(
                f"The nuclear binding energy may only be calculated for nucleons and isotopes.")

        number_of_protons = self.atomic_number
        number_of_neutrons = self.mass_number - self.atomic_number

        mass_of_protons = number_of_protons * const.m_p
        mass_of_neutrons = number_of_neutrons * const.m_n

        mass_of_nucleons = mass_of_protons + mass_of_neutrons

        mass_defect = mass_of_nucleons - self.nuclide_mass
        nuclear_binding_energy = mass_defect * const.c ** 2

        return nuclear_binding_energy.to(u.J)

    @property
    def half_life(self) -> Union[u.Quantity, str]:
        """
        Return the particle's half-life in seconds, or a `str`
        with half-life information.

        Particles that do not have sufficiently well-constrained
        half-lives will return a `str` containing the information
        that is available about the half-life and issue a
        `~plasmapy.utils.MissingAtomicDataWarning`.

        Examples
        --------
        >>> neutron = Particle('n')
        >>> neutron.half_life
        <Quantity 881.5 s>

        """
        if self.element and not self.isotope:
            raise InvalidIsotopeError(_category_errmsg(self.particle, 'isotope'))

        if isinstance(self._attributes['half-life'], str):
            warnings.warn(
                f"The half-life for {self.particle} is not known precisely; "
                "returning string with estimated value.", MissingAtomicDataWarning)

        if self._attributes['half-life'] is None:
            raise MissingAtomicDataError(f"The half-life of '{self.particle}' is not available.")
        return self._attributes['half-life']

    @property
    def spin(self) -> Union[int, float]:
        """
        Return the spin of the particle.

        If the spin is unavailable, then a
        `~plasmapy.utils.MissingAtomicDataError` will be raised.

        Examples
        --------
        >>> positron = Particle('e+')
        >>> positron.spin
        0.5

        """
        if self._attributes['spin'] is None:
            raise MissingAtomicDataError(f"The spin of particle '{self.particle}' is unavailable.")

        return self._attributes['spin']

    @property
    def periodic_table(self) -> collections.namedtuple:
        """
        Return a `~collections.namedtuple` to access category, period,
        group, and block information about an element.

        If the particle is not an element, isotope, or ion, then this
        attribute will raise an `~plasmapy.utils.InvalidElementError`.

        Examples
        --------
        >>> gold = Particle('Au')
        >>> gold.periodic_table.category
        'transition metal'
        >>> gold.periodic_table.period
        6
        >>> gold.periodic_table.group
        11
        >>> gold.periodic_table.block
        'd'

        """
        if self.element:
            return self._attributes['periodic table']
        else:  # coveralls: ignore
            raise InvalidElementError(_category_errmsg(self.particle, 'element'))

    @property
    def categories(self) -> Set[str]:
        """Return the particle's categories.

        Examples
        --------
        >>> gold = Particle('Au')
        >>> 'transition metal' in gold.categories
        True
        >>> 'antilepton' in gold.categories
        False

        """
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

        Examples
        --------
        >>> Particle('e-').is_category('lepton')
        True
        >>> Particle('p+').is_category('baryon', exclude='charged')
        False
        >>> Particle('n').is_category({'matter', 'baryon'}, exclude={'charged'})
        True
        >>> Particle('mu+').is_category('antilepton', exclude='baryon')
        True

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

        if category_tuple != () and require != set():  # coveralls: ignore
            raise AtomicError(
                "No positional arguments are allowed if the `require` keyword "
                "is set in is_category.")

        require = become_set(category_tuple) if category_tuple else become_set(require)

        exclude = become_set(exclude)
        any_of = become_set(any_of)

        if not require and not exclude and not any_of:
            return _valid_categories

        invalid_categories = (require | exclude | any_of) - _valid_categories

        duplicate_categories = require & exclude | exclude & any_of | require & any_of

        categories_and_adjectives = [
            (invalid_categories, 'invalid'),
            (duplicate_categories, 'duplicated'),
        ]

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

    @property
    def is_electron(self) -> bool:
        """
        Return `True` if the particle is an electron, and `False`
        otherwise.

        Examples
        --------
        >>> Particle('e-').is_electron
        True
        >>> Particle('e+').is_electron
        False

        """
        return self == "e-"

    @property
    def is_ion(self) -> bool:
        """
        Return `True` if the particle is an ion, and `False`
        otherwise.

        Examples
        --------
        >>> Particle('D+').is_ion
        True
        >>> Particle('H-1 0+').is_ion
        False
        >>> Particle('e+').is_ion
        False

        """
        return self.is_category('ion')
