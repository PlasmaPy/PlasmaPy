r"""The Particle class and particle_input decorator."""

from typing import (Union, Set, Tuple, List, Optional)
import warnings

from astropy import units as u, constants as const
import numpy as np

from ..utils import (
    AtomicError,
    AtomicWarning,
    InvalidParticleError,
    InvalidElementError,
    InvalidIsotopeError,
    InvalidIonError,
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
from .special import _Particles, ParticleZoo

# TODO: Write a decorator to turn atomic inputs into a Particle.


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

    Attributes
    ----------
    particle : str
        The particle symbol.

    element : str
        The atomic symbol, or None when the particle is not an element.

    isotope : str
        The isotope symbol, or None when the particle is not an isotope.

    ion : str
        The ion symbol, or None when the particle is not an ion.

    element_name : str
        The name of the element.

    integer_charge : int
        The charge in units of the elementary charge.

    charge : Quantity
        The charge in units of coulombs.

    mass : Quantity
        The mass of the particle, element, isotope, or ion.

    standard_atomic_weight : Quantity
        The standard atomic weight of an element, if available.

    nuclide_mass : Quantity
        The mass of a nucleon or of the nucleus of an isotope.

    atomic_number : int
        The atomic number of an element.

    mass_number : int
        The mass number of an isotope.

    baryon_number : int
        The number of baryons (protons and neutrons) minus the number of
        antibaryons (antiprotons and antineutrons) in the particle.

    lepton_number : int
        The number of leptons minus the number of antileptons for special
        particles and nuclides.

    binding_energy : Quantity
        The nuclear binding energy.

    half_life : Quantity
        The half-life of the particle or isotope in seconds.

    spin : int or float
        The spin of the particle, if available.

    reduced_mass
        Returns the reduced mass of the particle and another particle.

    is_category
        Tests whether or not the particle is in or not in certain categories.
        Possible categories include: 'lepton', 'antilepton', 'fermion',
        'boson', 'baryon', 'neutrino', 'antineutrino', 'element', 'isotope',
        'ion', 'matter', 'antimatter', 'stable', and 'unstable'.

    reduced_mass
        Returns the reduced mass of the particle and another particle.

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

        # If the data is missing, then the private attribute should still
        # exist but just be set to None.  This initialization had previously
        # been done in a loop using exec on a string, but this does not play
        # well with static type checkers such as PyCharm.

        self._particle_symbol = None
        self._atomic_symbol = None
        self._isotope_symbol = None
        self._ion_symbol = None
        self._atomic_symbol = None
        self._isotope_symbol = None
        self._ion_symbol = None
        self._unicode_symbol = None
        self._element_name = None
        self._atomic_number = None
        self._mass_number = None
        self._lepton_number = None
        self._baryon_number = None
        self._integer_charge = None
        self._electric_charge = None
        self._standard_atomic_weight = None
        self._mass = None
        self._nuclide_mass = None
        self._half_life = None
        self._spin = None
        self._generation = None
        self._periodic_table_group = None
        self._periodic_table_period = None
        self._charge = None
        self._electric_charge = None

        # Use this set to keep track of particle categories such as 'lepton'
        # for use with the is_category method later on.

        self._categories = set()

        # If the argument corresponds to one of the numerous case-sensitive or
        # case-insensitive aliases for particles, return the standard symbol.
        # Otherwise, just return the original argument.

        particle = _dealias_particle_aliases(argument)

        if particle in _Particles.keys():  # special particles
            self._particle_symbol = particle
            self._name = _Particles[particle]['name']
            self._spin = _Particles[particle]['spin']
            self._class = _Particles[particle]['class']
            self._lepton_number = _Particles[particle]['lepton number']
            self._baryon_number = _Particles[particle]['baryon number']
            self._integer_charge = _Particles[particle]['charge']
            self._half_life = _Particles[particle]['half-life']
            self._mass = _Particles[particle].get('mass', None)

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

                self._atomic_symbol = 'H'
                self._atomic_number = 1
                self._mass_number = 1
                self._element_name = 'hydrogen'
                self._isotope_symbol = 'H-1'
                self._ion_symbol = 'p+'
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
                    f"use:  Particle({repr(self._particle_symbol)})")

        else:  # elements, isotopes, and ions (besides protons)
            try:

                atomic_nomenclature_dict = _parse_and_check_atomic_input(
                    argument, mass_numb=mass_numb, Z=Z)

                self._particle_symbol = atomic_nomenclature_dict['symbol']
                self._atomic_symbol = atomic_nomenclature_dict['element']
                self._isotope_symbol = atomic_nomenclature_dict['isotope']
                self._ion_symbol = atomic_nomenclature_dict['ion']
                self._mass_number = atomic_nomenclature_dict['mass_numb']
                self._integer_charge = atomic_nomenclature_dict['Z']

            except Exception as exc:
                errmsg = _invalid_particle_errmsg(
                    argument, mass_numb=mass_numb, Z=Z)
                raise InvalidParticleError(errmsg) from exc

            element = self._atomic_symbol
            isotope = self._isotope_symbol
            ion = self._ion_symbol

            if element:
                self._categories.add('element')
            if isotope:
                self._categories.add('isotope')
            if ion:
                self._categories.add('ion')

            # Element properties

            self._atomic_number = _Elements[element]['atomic_number']
            self._element_name = _Elements[element]['name']

            self._baryon_number = self._mass_number

            # For the moment, set the lepton number to zero for elements,
            # isotopes, and ions.  The lepton number will probably come up
            # primarily during

            self._lepton_number = 0

            if isotope:
                if _Isotopes[isotope]['is_stable']:
                    self._half_life = np.inf * u.s
                else:
                    self._half_life = _Isotopes[isotope].get('half_life', None)
            elif element and not isotope:
                self._half_life = None

            if ion == 'He-4 2+':
                self._spin = 0
                self._categories.add('boson')

            # Set the masses

            if element and not isotope and not ion:
                try:
                    self._standard_atomic_weight = \
                        _Elements[element]['atomic_mass'].to(u.kg)
                    self._mass = self._standard_atomic_weight
                except KeyError:
                    self._standard_atomic_weight = None
            elif element and isotope and not ion:
                self._isotope_mass = _Isotopes[isotope]['atomic_mass']
            if isotope:
                self._isotope_mass = \
                    _Isotopes[isotope].get('atomic_mass', None).to(u.kg)
                self._standard_atomic_weight = None
                self._mass = self._isotope_mass
            else:
                self._standard_atomic_weight = \
                    _Elements[element].get('atomic_mass', None)
                self._isotope_mass = None

        # Set the charge

        if self._integer_charge is not None:
            self._electric_charge = self._integer_charge * const.e.si

        self._is_element = self._atomic_symbol is not None
        self._is_isotope = self._isotope_symbol is not None
        self._is_ion = self._ion_symbol is not None

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
        r"""Returns the particle symbol.

        Examples
        --------
        >>> electron = Particle('electron')
        >>> electron.particle
        'e-'
        >>> alpha = Particle('alpha')
        >>> alpha.particle
        'He-4 2+'
        """
        return self._particle_symbol

    @property
    def element(self) -> Optional[str]:
        r"""Returns the atomic symbol if the particle corresponds to an
        element, and None otherwise.

        Examples
        --------
        >>> alpha = Particle('He-4 2+')
        >>> alpha.element
        'He'
        >>> proton = Particle('proton')
        >>> electron = Particle('electron')
        >>> print(electron.element)
        None
        >>> if proton.element: print(proton.element)  # can use as conditional
        H
        """
        return self._atomic_symbol

    @property
    def isotope(self) -> Optional[str]:
        r"""Returns the isotope symbol if the particle corresponds to an
        isotope, and None otherwise.

        Example
        -------
        >>> alpha = Particle('alpha')
        >>> alpha.isotope
        'He-4'
        """
        return self._isotope_symbol

    @property
    def ion(self) -> Optional[str]:
        r"""Returns the ion symbol if the particle corresponds to an ion,
        and None otherwise.

        Example
        -------
        >>> alpha = Particle('alpha')
        >>> alpha.ion
        'He-4 2+'
        """
        return self._ion_symbol

    @property
    def element_name(self) -> str:
        r"""Returns the name of the element corresponding to this particle,
        or raises an InvalidElementError if the particle does not correspond
        to an element.

        Example
        -------
        >>> deuterium = Particle('D')
        >>> deuterium.element_name
        'hydrogen'
        """
        if not self.element:
            raise InvalidElementError(self._element_errmsg)
        return self._element_name

    @property
    def integer_charge(self) -> int:
        r"""Returns the integer charge of the partile, or raises a ChargeError
        if the charge has not been specified.

        Example
        -------
        >>> alpha = Particle('alpha')
        >>> alpha.integer_charge
        2
        """
        if self._integer_charge is None:
            raise ChargeError(
                f"The charge of particle {self.particle} has not been "
                f"specified.")
        return self._integer_charge

    @property
    def charge(self) -> u.C:
        r"""Returns the electric charge as a Quantity in units of coulombs,
        or raises a ChargeError if the charge has not been specified.

        Example
        -------
        >>> alpha = Particle('alpha')
        >>> alpha.charge
        <Quantity 3.20435324e-19 C>
        """
        if self._electric_charge is None:
            raise ChargeError(
                f"The charge of particle {self.particle} has not been "
                f"specified.")
        return self._electric_charge

    @property
    def mass(self) -> u.kg:
        r"""Returns the mass of the element, isotope, ion, particle, or
        antiparticle; or raises a MissingAtomicDataError if the mass
        is unavailable.

        Notes
        -----
        For special particles, this attribute will return the standard value
        of the mass of the particle.  The masses of neutrinos are not
        available.

        If the particle is an element and not an isotope or ion, then this
        attribute will return the standard atomic weight if available.

        If the particle is an isotope but not an ion, then this attribute
        will return the isotopic mass.

        If this particle is an ion, then this attribute will return the mass
        of the element or isotope (as described above) minus the integer charge
        times the electron mass.

        Example
        -------
        >>> alpha = Particle('alpha')
        >>> alpha.mass
        <Quantity 6.64647897e-27 kg>
        """
        if self._mass is None:
            raise MissingAtomicDataError(
                f"The mass of particle '{self.particle}' is unavailable.")
        return self._mass.to(u.kg)

    @property
    def standard_atomic_weight(self) -> u.kg:
        r"""Returns the standard atomic weight of an element if available.
        Raises a MissingAtomicDataError if the particle is an element for
        which the standard_atomic_weight is unavailable.  Raises an
        InvalidElementError if the particle is not an element.

        Example
        -------
        >>> H = Particle('H')
        >>> H.standard_atomic_weight
        <Quantity 1.67382335e-27 kg>
        """
        if self.element and not self.isotope and not self.ion:
            if self._standard_atomic_weight is None:
                raise MissingAtomicDataError(
                    f"The standard atomic weight of {self.element} is "
                    f"unavailable.")
            else:
                return self._standard_atomic_weight.to(u.kg)
        else:
            raise InvalidElementError(self._element_errmsg)

    @property
    def nuclide_mass(self) -> u.kg:
        r"""Returns the mass of the nucleus of an isotope, or raises an
        InvalidIsotopeError if the particle is not an isotope or neutron.

        Example
        -------
        >>> isotope = Particle('O-18')
        >>> isotope.nuclide_mass
        <Quantity 2.98810197e-26 kg>
        """
        if self.particle in ['H-1', 'p+']:
            _nuclide_mass = const.m_p
        elif self.particle == 'n':
            _nuclide_mass = const.m_n
        elif self._is_isotope and not self._is_ion:
            try:
                _atomic_number = self._atomic_number
                _isotope_mass = _Isotopes[self.isotope]['atomic_mass'].to(u.kg)
                _nuclide_mass = _isotope_mass - _atomic_number * const.m_e
            except KeyError:  # coveralls: ignore
                raise MissingAtomicDataError(
                    f"The mass of a {self.isotope} nuclide is not available.")
        else:
            raise InvalidIsotopeError(self._isotope_errmsg)

        return _nuclide_mass

    @property
    def atomic_number(self) -> int:
        r"""Returns the atomic number of the element corresponding to this
        particle, or raises an InvalidElementError if the particle does not
        correspond to an element.

        Example
        -------
        >>> iron = Particle('Fe')
        >>> iron.atomic_number
        26
        """
        if not self._is_element:
            raise InvalidElementError(self._element_errmsg)
        return self._atomic_number

    @property
    def mass_number(self) -> int:
        r"""Returns the mass number of the isotope corresponding to this
        particle, or raises an InvalidIsotopeError if the particle does not
        correspond to an isotope.

        Example
        -------
        >>> tritium = Particle('H-3')
        >>> tritium.mass_number
        3
        """
        if not self._is_isotope:
            raise InvalidIsotopeError(self._isotope_errmsg)
        return self._mass_number

    @property
    def baryon_number(self) -> int:
        r"""Returns the number of protons plus neutrons minus the number of
        antiprotons and antineutrons in the particle, or raises an
        AtomicError if the baryon number is unavailable.  The baryon number
        is equivalent to the mass number for isotopes.

        Example
        -------
        >>> antineutron = Particle('antineutron')
        >>> electron = Particle('electron')
        >>> alpha = Particle('He-4 2+')
        >>> antineutron.baryon_number
        -1
        >>> electron.baryon_number
        0
        >>> alpha.baryon_number
        4
        >>> alpha.baryon_number == alpha.mass_number
        True
        """
        if self._baryon_number is None:  # coveralls: ignore
            raise AtomicError(
                f"The baryon number for '{self.particle}' is not "
                f"available.")
        return self._baryon_number

    @property
    def lepton_number(self) -> int:
        r"""Returns 1 for leptons, -1 for antileptons, and 0 for
        nuclides/isotopes; or raises an AtomicError if the lepton number is
        not available.  This attribute does not include the electrons in
        an atom or ion.

        Example
        -------
        >>> positron = Particle('e+')
        >>> electron = Particle('e-')
        >>> proton = Particle('p')
        >>> positron.lepton_number
        -1
        >>> electron.lepton_number
        1
        >>> proton.lepton_number
        0
        """
        if self._lepton_number is None:  # coveralls: ignore
            raise AtomicError(
                f"The lepton number for {self.particle} is not available.")
        return self._lepton_number

    @property
    def binding_energy(self) -> u.J:
        r"""Returns the nuclear binding energy, or raises an
        InvalidIsotopeError if the particle is not a nucleon or isotope.

        Example
        -------
        >>> iron59 = Particle('Fe-59')
        >>> iron59.binding_energy
        <Quantity 8.27574334e-11 J>
        """

        if self._baryon_number == 1:
            return 0 * u.J

        if not self.element:
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
        MissingAtomicDataError if the half-life is unavailable.

        Examples
        --------
        >>> from astropy import units as u
        >>> neutron = Particle('n')
        >>> neutron.half_life
        <Quantity 881.5 s>
        >>> carbon14 = Particle('C-14')
        >>> carbon14.half_life.to(u.yr)
        <Quantity 5730. yr>

        """
        if self._atomic_symbol and not self._isotope_symbol:
            raise InvalidIsotopeError(self._isotope_errmsg)
        if self._half_life is None:
            raise MissingAtomicDataError(
                f"The half-life of '{self.particle}' is not available.")
        return self._half_life

    @property
    def spin(self) -> Union[int, float]:
        r"""Returns the spin of the particle, or raises a
        MissingAtomicDataError if the spin is not available.

        Example
        -------
        >>> electron = Particle('e-')
        >>> electron.spin
        0.5
        """
        if self._spin is None:
            raise MissingAtomicDataError(
                f"The spin of particle '{self.particle}' is unavailable.")
        return self._spin

    def reduced_mass(self, other, Z=None, mass_numb=None) -> u.kg:
        r"""Finds the reduced mass between two particles, or will raise a
        MissingAtomicDataError if either particle's mass is unavailable or
        an AtomicError for any other errors.  The other particle may be
        represented by another Particle object, a Quantity with units of,
        mass, or a string of the other particle's symbol (in conjunction
        with keywords Z and mass_numb)

        Example
        -------
        >>> from plasmapy.atomic import Particle
        >>> electron = Particle('e-')
        >>> proton = Particle('p+')
        >>> proton.reduced_mass(electron)
        <Quantity 9.10442514e-31 kg>

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

    def is_category(self, *categories, any=False,
                    exclude: Union[Set, Tuple, List] = set()) -> bool:
        r"""Returns True if the particle is in all of the inputted categories.

        If any is True, then the particle will return True if any of the listed
        categories are met.

        If the exclude keyword is set, then is_category will return False if
        the particle is in any of the excluded categories, whether or not the
        the particle matches the other criteria.

        The valid categories are: 'lepton', 'antilepton', 'baryon',
        'antibaryon', 'fermion', 'boson', 'neutrino', 'antineutrino', 'matter',
        'antimatter', 'element', 'isotope', and 'ion'.

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
            'antimatter', 'stable', 'unstable'
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


def particle_input():
    pass
