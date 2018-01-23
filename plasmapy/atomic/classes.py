from typing import Union
from astropy import units as u, constants as const
import numpy as np

from ..utils import (
    AtomicError,
    InvalidParticleError,
    InvalidElementError,
    InvalidIsotopeError,
    InvalidIonError,
    ChargeError,
    MissingAtomicDataError,
)

from .particles import (
    _Particles,
    _special_particles,
    _leptons,
    _antileptons,
    _fermions,
    _bosons,
    _neutrinos,
    _antineutrinos,
)

from .elements import _Elements
from .isotopes import _Isotopes

from .parsing import (
    _dealias_particle_aliases,
    _parse_and_check_atomic_input,
    _invalid_particle_errmsg,
)


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
        The atomic symbol.

    isotope : str
        The isotope symbol.

    ion : str
        The ion symbol.

    atomic_number : int
        The atomic number of an element.

    mass_number : int
        The mass number of an isotope.

    mass : Quantity
        The mass of the particle, element, or ion

    integer_charge :

    charge :



    half_life : Quantity

    is_stable : bool
        Returns True if the particle is stable and False if it is unstable.



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

    InvalidIonError
        For when an attribute is being accessed that requires information
        about an ion (excluding charge information), but the particle is not
        an ion.

    ChargeError
        For when either the charge or integer_charge attributes is being
        accessed but the charge information for the particle is not
        available.

    TypeError
        For when any of the arguments or keywords is not of the required
        type.

    """

    # TODO: Write an actual docstring of wonder and amazement
    # TODO: Write a decorator to turn atomic inputs into a Particle.
    # TODO: Add an attribute for nuclide_mass

    def __init__(self,
                 argument: Union[str, int],
                 mass_numb: int = None,
                 Z: int = None):

        self._particle_symbol = None
        self._atomic_symbol = None
        self._isotope_symbol = None
        self._ion_symbol = None
        self._unicode_symbol = None
        self._element_name = None
        self._atomic_number = None
        self._mass_number = None
        self._lepton_number = None
        self._baryon_number = None
        self._is_stable = None
        self._is_antimatter = None
        self._integer_charge = None
        self._electric_charge = None
        self._standard_atomic_weight = None
        self._mass = None
        self._half_life = None
        self._spin = None
        self._lepton = None
        self._antilepton = None
        self._baryon = None
        self._antibaryon = None
        self._fermion = None
        self._generation = None
        self._periodic_table_group = None
        self._periodic_table_period = None

        particle = _dealias_particle_aliases(argument)

        is_special_particle = particle in _special_particles

        if is_special_particle:
            self._particle_symbol = particle
            self._name = _Particles[particle]['name']
            self._spin = _Particles[particle]['spin']
            self._class = _Particles[particle]['class']
            self._lepton_number = _Particles[particle]['lepton number']
            self._baryon_number = _Particles[particle]['baryon number']
            self._integer_charge = _Particles[particle]['charge']
            self._mass = _Particles[particle].get('mass', None)
            self._is_stable = _Particles[particle]['half-life'] == np.inf * u.s
            self._half_life = _Particles[particle]['half-life']
            self._is_antimatter = _Particles[particle]['antimatter']
            self._is_lepton = particle in _leptons
            self._is_antilepton = particle in _antileptons
            self._is_fermion = particle in _fermions
            self._is_boson = particle in _bosons

            self._is_neutrino_or_antineutrino = \
                particle in _neutrinos + _antineutrinos

            if particle in _leptons + _antileptons:
                self._generation = _Particles[particle]['generation']

            if mass_numb is not None or Z is not None:
                raise InvalidParticleError(
                    "The keywords 'mass_numb' and 'Z' cannot be used when "
                    "creating Particle objects for special particles. To "
                    f"create a Particle object for {self._name}s, "
                    f"use:  Particle({repr(self._particle_symbol)})")
        else:
            try:
                particle_dict = _parse_and_check_atomic_input(
                    argument, mass_numb=mass_numb, Z=Z)
            except Exception as exc:
                errmsg = _invalid_particle_errmsg(
                    argument, mass_numb=mass_numb, Z=Z)
                raise InvalidParticleError(errmsg) from exc

            self._particle_symbol = particle_dict['symbol']

            element = particle_dict['element']
            isotope = particle_dict['isotope']
            ion = particle_dict['ion']

            self._atomic_symbol = element
            self._isotope_symbol = isotope
            self._ion_symbol = ion

            self._integer_charge = particle_dict['Z']
            self._mass_number = particle_dict['mass_numb']
            self._baryon_number = self._mass_number
            self._atomic_number = _Elements[element]['atomic_number']
            self._element_name = _Elements[element]['name']

            if element:
                self._lepton_number = 0
            if isotope:
                try:
                    self._half_life = _Isotopes[isotope]['half-life']
                except KeyError:
                    self._half_life = np.inf * u.s
            elif element and not isotope:
                self._half_life = None

            if ion == 'p+':
                self._mass = const.m_p
                self._spin = 1/2
                self._lepton_number = 0
            elif ion == 'He-4 2+':
                self._spin = 0

            elif element and not isotope and not ion:
                try:
                    self._standard_atomic_weight = \
                        _Elements[element]['atomic_mass'].to(u.kg)
                    self._mass = self._standard_atomic_weight
                except KeyError:
                    self._standard_atomic_weight = None
            elif element and isotope and not ion:
                self._isotope_mass = _Isotopes[isotope]['atomic_mass']

            self._is_antimatter = False

        if self._integer_charge is not None:
            self._electric_charge = self._integer_charge * const.e.si
        else:
            self._electric_charge = None

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

    def __repr__(self):
        r"""Returns a string of the call that would recreate this object."""
        return f'Particle("{self.particle}")'

    def __str__(self):
        r"""Returns a string of the particle symbol."""
        return f"{self.particle}"

    def __eq__(self, other):
        r"""Returns True when comparing two Particle objects that correspond
        to the same particle, and False when the two objects differ."""
        try:
            if self.__dict__ == other.__dict__:
                return True
            else:
                return False
        except Exception:
            return False

    def __ne__(self, other):
        r"""Returns True when the two objects differ, and False when
        comparing two Particle objects that correspond to the same particle."""
        return not self.__eq__(other)

    @property
    def particle(self) -> str:
        r"""Returns the particle symbol."""
        return self._particle_symbol

    @property
    def element(self) -> str:
        r"""Returns the atomic symbol, or raises an InvalidElementError
        if the particle is not an element, isotope, or ion."""
        if not self._atomic_symbol:
            raise InvalidElementError(self._element_errmsg)
        return self._atomic_symbol

    @property
    def isotope(self) -> str:
        r"""Returns the isotope symbol, or raises an InvalidIsotopeError
        if the particle does not correspond to an isotope."""
        if not self._is_isotope:
            raise InvalidIsotopeError(self._isotope_errmsg)
        return self._isotope_symbol

    @property
    def ion(self) -> str:
        r"""Returns the ion symbol, or raises an InvalidIonError if the
        particle is not an ion."""
        if not self._ion_symbol:
            raise InvalidIonError(self._ion_errmsg)
        return self._ion_symbol

    @property
    def element_name(self) -> str:
        r"""Returns the name of the element, or raises an InvalidElementError
        if the particle is not an element, isotope, or ion."""
        if not self._is_element:
            raise InvalidElementError(self._element_errmsg)
        return self._element_name

    @property
    def atomic_number(self) -> int:
        r"""Returns the atomic number of the particle if it is an element,
        ion, or isotope."""
        if not self._is_element:
            raise InvalidElementError(self._element_errmsg)
        return self._atomic_number

    @property
    def mass_number(self) -> int:
        r"""Returns the mass number of an isotope, or raises an
        InvalidIsotopeError if the particle does not have a mass number."""
        if not self._is_isotope:
            raise InvalidIsotopeError(self._isotope_errmsg)
        return self._mass_number

    @property
    def baryon_number(self) -> int:
        r"""Returns the baryon number of a subatomic particle or the nuclide
        of an element or isotope, or raises a MissingAtomicDataError if
        the baryon number is unavailable."""
        if self._baryon_number is None:  # coveralls: ignore
            raise AtomicError(
                f"The baryon number for '{self.particle}' is not "
                f"available.")
        return self._baryon_number

    @property
    def lepton_number(self) -> int:
        r"""Returns the lepton number of a subatomic particle or the nuclide
        of an element or isotope, or raises an AtomicError if the lepton
        number is not available."""
        if self._lepton_number is None:  # coveralls: ignore
            raise AtomicError(
                f"The lepton number for {self.particle} is not available.")
        return self._lepton_number

    @property
    def half_life(self) -> u.s:
        r"""Returns the half-life of the particle, or raises a
        MissingAtomicDataError if the half-life is unavailable."""
        if self._atomic_symbol and not self._isotope_symbol:
            raise InvalidIsotopeError(self._isotope_errmsg)
        if not self._half_life:
            raise MissingAtomicDataError(
                f"The half-life of '{self.particle}' is not available.")
        return self._half_life

    @property
    def is_stable(self) -> bool:
        r"""Returns True if the particle is stable and False if the
        particle is unstable, or raises a MissingAtomicDataError if
        stability information is not available."""
        if self._atomic_symbol and not self._isotope_symbol:
            raise InvalidIsotopeError(self._isotope_errmsg)
        if not self._half_life:
            raise MissingAtomicDataError(
                f"The stability of '{self.particle}' is not available.")
        return self._half_life == np.inf * u.s

    @property
    def is_antimatter(self) -> bool:
        r"""Returns True is the particle is antimatter and False if the
        particle is matter."""
        return self._is_antimatter

    @property
    def integer_charge(self) -> int:
        r"""Returns the integer charge, or raises a ChargeError if the
        charge has not been specified."""
        if self._integer_charge is None:
            raise ChargeError(
                f"The charge of particle {self.particle} has not been "
                f"specified.")
        return self._integer_charge

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
        is unavailable."""
        # TODO: Take care of the cases in isotope_mass and ion_mass
        if self._mass is None:
            raise MissingAtomicDataError(
                f"The mass of particle '{self.particle}' is unavailable.")
        return self._mass

    @property
    def standard_atomic_weight(self) -> u.kg:
        r"""Returns the standard atomic weight of an element if available.
        Raises a MissingAtomicDataError if the particle is an element for
        which the standard_atomic_weight is unavailable.  Raises an
        InvalidElementError if the particle is not an element."""
        if (self._is_element and not self._is_isotope and not self._is_ion):
            if self._standard_atomic_weight is None:
                raise MissingAtomicDataError(
                    f"The standard atomic weight of {self.element} is "
                    f"unavailable.")
            else:
                return self._standard_atomic_weight.to(u.kg)
        else:
            raise InvalidElementError(self._element_errmsg)

    @property
    def spin(self) -> Union[int, float]:
        r"""Returns the spin of the Particle, or raises a
        MissingAtomicDataError if the spin is not available."""
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
        except MissingAtomicDataError:  # coveralls: ignore
            raise MissingAtomicDataError(
                f"Unable to find the reduced mass because the mass of "
                f"{self.particle} is not available.") from None

        if isinstance(other, (str, int)):
                other = Particle(other, Z=Z, mass_numb=mass_numb)

        if isinstance(other, Particle):
            try:
                mass_that = other.mass.to(u.kg)
            except MissingAtomicDataError:  # coveralls: ignore
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
