from typing import Union
from astropy import units as u, constants as const

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

from .elements import (_Elements, _atomic_symbols, _atomic_symbols_dict)

from .parsing import (
    _dealias_particle_aliases,
    _parse_and_check_atomic_input,
)


class Particle():
    r"""A class for individual particles or antiparticles."""

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
            raise InvalidParticleError()

        if is_special_particle:
            self._particle_symbol = particle
            self._name = _Particles[particle]['name']
            self._spin = _Particles[particle]['spin']
            self._class = _Particles[particle]['class']
            self._lepton_number = _Particles[particle]['lepton number']
            self._baryon_number = _Particles[particle]['baryon number']
            self._integer_charge = _Particles[particle]['charge']
            self._mass = _Particles[particle]['mass']

            self._is_stable = _Particles[particle]['is_stable']

            self._half_life = _Particles[particle]['half-life']
            self._is_antimatter = _Particles[particle]['antimatter']

            if particle in _leptons + _antileptons:
                self._is_lepton = True
                self._generation = _Particles[particle]['generation']
            else:
                self._is_lepton = False

            if particle in _fermions:
                self._is_fermion = True
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
                    _Elements[element]['atomic_mass']
                self._mass = self._standard_atomic_weight
            elif element and isotope and not ion:
                self._isotope_mass = _Isotopes[isotope]['atomic_mass']

            self._is_antimatter = False

        if self._integer_charge is not None:
            self._electric_charge = self._integer_charge * const.e.si
        else:
            self._electric_charge = None

    @property
    def particle(self):
        return self._particle_symbol

    @property
    def element(self):
        return self._atomic_symbol

    @property
    def isotope(self):
        return self._isotope_symbol

    @property
    def ion(self):
        return self._ion_symbol

    @property
    def atomic_number(self):
        if self._atomic_number is not None:
            return self._atomic_number
        else:
            raise InvalidElementError(
                f"Particle '{self._particle_symbol}' is not an element, "
                f"so no atomic number is available.")

    @property
    def mass_number(self):
        if self._mass_number is not None:
            return self._mass_number
        else:
            raise InvalidIsotopeError(
                f"Particle '{self._particle_symbol}' is not an isotope, so "
                f"no mass number is available.")

    @property
    def lepton_number(self):
        if self._lepton_number is not None:
            return self._lepton_number
        else:
            raise AtomicError(f"The leption number for "
                              f"'{self._particle_symbol}' is not available.")

    @property
    def baryon_number(self):
        if self._baryon_number is not None:
            return self._baryon_number
        elif self._mass_number is not None:
            return self._mass_number
        else:
            raise AtomicError("The baryon number for "
                              f"'{self._particle_symbol}' is not available.")

    @property
    def half_life(self):
        if self._half_life is not None:
            return self._half_life
        else:
            raise MissingAtomicDataError(
                f"Half-life data for '{self._particle_symbol}' is not "
                "available.")

    @property
    def is_stable(self):
        if self._half_life is not None:
            return self._half_life == np.inf * u.s
        else:
            raise MissingAtomicDataError(
                f"Stability data for '{self._particle_symbol}' is not "
                "available.")

    @property
    def is_antimatter(self):
        return self._is_antimatter

    @property
    def Z(self):
        if self._integer_charge is not None:
            return self._integer_charge
        else:
            raise ChargeError(
                f"No charge information is available for particle "
                f"'{self._particle_symbol}'.")

    @property
    def q(self):
        if self._electric_charge is not None:
            return self._electric_charge
        else:
            raise ChargeError(
                f"No charge information is available for particle "
                f"'{self._particle_symbol}'.")

    @property
    def m(self):
        if self._mass is not None:
            return self._mass
        else:
            raise MissingAtomicDataError(
                f"The mass of particle '{self._particle_symbol}' is "
                f"not available.")

    def __repr__(self):
        return f'<Particle "{self._particle_symbol}">'

    def __str__(self):
        return f"{self._particle_symbol}"
