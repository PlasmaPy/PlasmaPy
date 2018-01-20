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
    _baryons,
    _antibaryons,
    _particles,
    _antiparticles,
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
)


class Particle():
    r"""A class for individual particles or antiparticles."""

    def __repr__(self):
        return f'<Particle "{self._particle_symbol}">'

    def __str__(self):
        return f"{self._particle_symbol}"

    def __init__(self,
                 argument: Union[str, int],
                 mass_numb: int = None,
                 Z: int = None):

        self._original_argument = argument
        self._original_mass_number = mass_numb
        self._original_integer_charge = Z
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

        if is_special_particle and (mass_numb is not None or Z is not None):
            raise InvalidParticleError(
                f"Particle '{argument}' with mass_numb = {mass_numb} "
                f"and Z = {Z} is not a valid particle.")

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

        else:
            try:
                particle_dict = _parse_and_check_atomic_input(
                    argument, mass_numb=mass_numb, Z=Z)
            except Exception as exc:
                raise InvalidParticleError(
                    f"Particle '{argument}' with mass_numb = {mass_numb} "
                    f"and Z = {Z} is not a valid particle.") from exc

            self._particle_symbol = particle_dict['symbol']

            element = particle_dict['element']
            isotope = particle_dict['isotope']
            ion = particle_dict['ion']

            self._atomic_symbol = element
            self._isotope_symbol = isotope
            self._ion_symbol = ion

            self._mass_number = particle_dict['mass_numb']

            self._integer_charge = particle_dict['Z']

            self._atomic_number = _Elements[element]['atomic_number']
            self._element_name = _Elements[element]['name']

            if element and not isotope and not ion:
                self._standard_atomic_weight = \
                    _Elements[element]['atomic_mass'].to(u.kg)
                self._mass = self._standard_atomic_weight
            elif element and isotope and not ion:
                self._isotope_mass = _Isotopes[isotope]['atomic_mass']

            self._is_antimatter = False

        if self._integer_charge is not None:
            self._electric_charge = self._integer_charge * const.e.si
        else:
            self._electric_charge = None

        self._element_errmsg = (
            f"The particle '{self.particle}' is not an element, so "
            f"this attribute is not available.")

        self._isotope_errmsg = (
            f"The particle '{self.particle}' does not have an "
            f"isotope specified, so this attribute is not available.")

        self._ion_errmsg = (
            f"The particle '{self.particle}' is not an ion, so this"
            f"attribute is not available.")

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
        if not self._isotope_symbol:
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
        if not self.element:
            raise InvalidElementError(self._element_errmsg)
        return self._element_name

    @property
    def atomic_number(self) -> int:
        if not self.element:
            raise InvalidElementError(self._element_errmsg)
        return self._atomic_number

    @property
    def mass_number(self) -> int:
        if not self.isotope:
            raise InvalidIsotopeError(self._isotope_errmsg)
        return self._mass_number

    @property
    def lepton_number(self) -> int:
        if self._lepton_number is None:
            raise AtomicError(
                f"The lepton number for {self.particle} is not available.")
        return self._lepton_number

    @property
    def baryon_number(self) -> int:
        r"""Returns the baryon number, or raises a MissingAtomicDataError if
        the baryon number is unavailable."""
        if self._baryon_number is None:
            raise AtomicError(
                f"The baryon number for '{self.particle}' is not "
                f"available.")
        return self._baryon_number

    @property
    def half_life(self) -> u.Quantity:
        r"""Returns the half-life of the particle, or raises a
        MissingAtomicDataError if the half-life is unavailable."""
        if not self._half_life:
            raise MissingAtomicDataError(
                f"The half-life of '{self.particle}' is not available.")
        return self._half_life

    @property
    def is_stable(self) -> bool:
        r"""Returns True if the particle is stable and False if the
        particle is unstable, or raises a MissingAtomicDataError if
        stability information is not available."""
        if not self._half_life:
            raise MissingAtomicDataError(
                f"The stability of '{self.particle}' is not available.")
        return self._half_life

    @property
    def is_antimatter(self) -> bool:
        r"""Returns True is the particle is antimatter and False if the
        particle is matter."""
        return self._is_antimatter

    @property
    def Z(self) -> int:
        r"""Returns the integer charge, or raises a ChargeError if the
        charge has not been specified."""
        if self._integer_charge is None:
            raise ChargeError(
                f"The charge of particle {self.particle} has not been "
                f"specified.")
        return self._integer_charge

    @property
    def q(self) -> u.Quantity:
        r"""Returns the electric charge as a Quantity in units of coulombs,
        or raises a ChargeError if the charge has not been specified."""
        if self._electric_charge is None:
            raise ChargeError(
                f"The charge of particle {self.particle} has not been "
                f"specified.")
        return self._electric_charge

    @property
    def m(self) -> u.Quantity:
        r"""Returns the mass of the element, isotope, ion, particle, or
        antiparticle; or raises a MissingAtomicDataError if the mass
        is unavailable."""
        # TODO: Take care of the cases in isotope_mass and ion_mass
        if self._mass is None:
            raise MissingAtomicDataError(
                f"The mass of particle '{self.particle}' is unavailable.")
        return self._mass

    @property
    def spin(self) -> Union[int, float]:
        r"""Returns the spin of the Particle, or raises a
        MissingAtomicDataError if the spin is not available."""
        if self._spin is None:
            raise MissingAtomicDataError(
                f"The spin of particle '{self.particle}' is unavailable.")
        return self._spin
