"""Classes for storing ionization state data."""

import numpy as np
import astropy.units as u
from typing import Union, Dict, List, Tuple
import collections

from .particle_class import Particle
from .particle_input import particle_input
from ..utils import (
    AtomicError,
    check_quantity,
    _check_quantity,
    InvalidParticleError,
    ChargeError,
)

State = collections.namedtuple('State', ['integer_charge', 'ionic_fraction', 'ionic_symbol'])

_number_density_errmsg = (
    "Number densities must be Quantity objects with units of inverse "
    "volume."
)


class IonizationState:
    """
    Describe the ionization state distribution of a single element.

    Parameters
    ----------
    particle: str, int, np.integer, or ~plasmapy.atomic.Particle
        A `str` or `~plasmapy.atomic.Particle` instance representing
        an element or isotope, or an `int` representing the atomic
        number of an element.

    ionic_fractions: ~numpy.ndarray, list, tuple, or ~astropy.units.Quantity
        The ionization fractions of an element, where the indices
        correspond to integer charge.  This argument should contain the
        atomic number plus one items, and must sum to one within an
        absolute tolerance of `tol`.  Alternatively, this argument may
        be a `~astropy.units.Quantity` that represents the number
        densities of each neutral/ion.

    T_e: ~astropy.units.Quantity, keyword-only, optional
        The electron temperature or thermal energy per particle.

    n_e: ~astropy.units.Quantity, keyword-only, optional
        The electron number density.

    n_elem: ~astropy.units.Quantity, keyword-only, optional
        The number density of ions plus atoms.

    tol: float or int, keyword-only, optional
        The absolute tolerance used by `~numpy.isclose` when testing
        normalizations and making comparisons.  Defaults to `1e-15`.

    Raises
    ------
    ~plasmapy.utils.AtomicError
        If the ionic fractions are not normalized or contain invalid
        values

    ~plasmapy.utils.InvalidParticleError
        If the particle is invalid.

    Notes
    -----
    Only one of `n_e` and `n_elem` may be set.

    """

    @check_quantity({
        "T_e": {"units": u.K, "none_shall_pass": True},
        "n_e": {"units": u.m ** -3, "none_shall_pass": True},
        "n_elem": {"units": u.m ** -3, "none_shall_pass": True},
    })
    @particle_input(require='element', exclude='ion')
    def __init__(self,
                 particle: Particle,
                 ionic_fractions=None,
                 *,
                 T_e=None,
                 n_e=None,
                 n_elem=None,
                 tol: Union[float, int] = 1e-15):
        """Initialize a `~plasmapy.atomic.IonizationState` instance."""

        try:
            self._n_e = None
            self._n_elem = None
            self.tol = tol
            self.T_e = T_e
            self.n_elem = n_elem
            self.n_e = n_e
            self.ionic_fractions = ionic_fractions
            self._particle = particle  # Store Particle class instance.
            if self._ionic_fractions is None and self.T_e is not None:
                self.equilibrate()
        except Exception as exc:
            raise AtomicError(
                f"Unable to create IonizationState instance for "
                f"{particle.particle}.")

    def __getitem__(self, value):
        """Return the ionic fraction(s)."""
        if isinstance(value, slice):
            result = State(
                np.arange(value.start, value.stop, value.step),
                self.ionic_fractions[value.start:value.stop:value.step],
                self.ionic_symbols[value.start:value.stop:value.step],
            )
        elif isinstance(value, (int, np.integer)) and 0 <= value <= self.atomic_number:
            result = State(value, self.ionic_fractions[value], self.ionic_symbols[value])
        else:
            if not isinstance(value, Particle):
                try:
                    value = Particle(value)
                except InvalidParticleError as exc:
                    raise InvalidParticleError(
                        f"{value} is not a valid integer charge or "
                        f"particle.") from exc

            same_element = value.element == self.element
            same_isotope = value.isotope == self.isotope
            has_charge_info = value.is_category(any_of=["charged", "uncharged"])

            if same_element and same_isotope and has_charge_info:
                Z = value.integer_charge
                result = State(Z, self.ionic_fractions[Z], self.ionic_symbols[Z])
            else:
                if not same_element or not same_isotope:
                    raise AtomicError("Inconsistent element or isotope.")
                elif not has_charge_info:
                    raise ChargeError("No integer charge provided.")
        return result

    def __iter__(self):
        """Initialize the instance prior to an iteration."""
        self._charge_index = 0
        return self

    def __next__(self):
        """
        Return a `~State` `~collections.namedtuple` which contains
        ``integer_charge``, ``ionic_fraction``, and ``ionic_symbol``.
        """
        if self._charge_index <= self.atomic_number:
            result = State(
                self._charge_index,
                self._ionic_fractions[self._charge_index],
                self.ionic_symbols[self._charge_index],
            )
            self._charge_index += 1
            return result
        else:
            del self._charge_index
            raise StopIteration

    def __eq__(self, other):
        """
        Return `True` if the ionic fractions for two `IonizationState`
        instances are approximately equal to within the minimum `tol`
        specified by either, and `False` otherwise.

        Raises
        ------
        AtomicError
            If `other` is not an `~plasmapy.atomic.IonizationState`
            instance, or if `other` corresponds to a different element.

        Examples
        --------
        >>> IonizationState('H', [1, 0], tol=1e-6) == IonizationState('H', [1, 1e-6], tol=1e-6)
        True
        >>> IonizationState('H', [1, 0], tol=1e-8) == IonizationState('H', [1, 1e-6], tol=1e-5)
        False

        """
        if not isinstance(other, IonizationState):
            raise AtomicError(
                "Instances of the IonizationState class may only be "
                "compared with other IonizationState instances.")

        if self.element != other.element:
            raise AtomicError("Only ionization states of the same element may be compared.")

        # Use the tightest of the two absolute tolerances
        min_tol = np.min([self.tol, other.tol])

        return np.allclose(self.ionic_fractions, other.ionic_fractions, atol=min_tol)

    @property
    def ionic_fractions(self) -> np.ndarray:
        """
        Return the ionic fractions, where the index corresponds to
        the integer charge.

        Examples
        --------
        >>> hydrogen_states = IonizationState('H', [0.9, 0.1])
        >>> hydrogen_states.ionic_fractions
        array([0.9, 0.1])

        """
        return self._ionic_fractions

    @ionic_fractions.setter
    def ionic_fractions(self, fractions):
        """
        Set the ionic fractions, while checking that the new values are
        valid and normalized to one.
        """
        try:
            if not isinstance(fractions, np.ndarray) or 'float' not in str(fractions.dtype):
                fractions = np.array(fractions, dtype=np.float)
        except Exception as exc:
            raise AtomicError(
                f"Unable to set ionic fractions of {self.element} "
                f"to {fractions}.") from exc

        if np.min(fractions) < 0:
            raise AtomicError("Ionic fractions cannot be negative.")

        if isinstance(fractions, u.Quantity):
            if self._n_e or self._n_elem:
                raise AtomicError(
                    "The ionization state may be set using number "
                    "densities for each ion only if neither of the "
                    "electron density and element density has already "
                    "been set.")
            self.number_densities = fractions
        else:

            total = np.sum(fractions)
            if not np.isclose(total, 1, atol=self.tol, rtol=0):
                raise AtomicError(
                    f"The sum of the ionic fractions of {self.element} "
                    f"equals {total}, which is not approximately one.")

            self._ionic_fractions = fractions

    @property
    def n_e(self):
        """
        Return the electron number density assuming a single species
        plasma.
        """
        if self._n_e is not None:
            return self._n_e
        elif self._n_elem is not None:
            return np.sum(self._n_elem * self.ionic_fractions * self.integer_charges)
        else:
            raise AtomicError(
                "Insufficient information to calculate the electron "
                "density.")

    @n_e.setter
    def n_e(self, value):
        """
        Return the electron density assuming a single-species
        plasma.
        """
        if value is None:
            self._n_e = None
            return
        elif self._n_elem is not None:
            raise AtomicError(
                "Only one of n_e and n_elem may be set for a "
                "single element, quasineutral plasma.")
        try:
            self._n_e = value.to(u.m ** -3)
            self._n_elem = None
        except (AttributeError, u.UnitConversionError):
            raise AtomicError(_number_density_errmsg)

    @property
    def n_elem(self):
        """
        Return the number density of atoms plus ions for this
        species.
        """
        if self._n_elem is not None:
            return self._n_elem
        elif self._n_e is not None:
            return self._n_e / (self.ionic_fractions * self.integer_charges)

    @n_elem.setter
    def n_elem(self, value):
        """The number density of atoms plus ions of this species."""
        if value is None:
            self._n_elem = None
        else:
            if '_n_e' in dir(self) and self._n_e is not None:
                raise AtomicError(
                    "Only one of n_e and n_elem may be set for a "
                    "single element, quasineutral plasma.")
            try:
                self._n_elem = value.to(u.m ** -3)
            except (AttributeError, u.UnitConversionError):
                raise AtomicError(_number_density_errmsg)

    @property
    def number_densities(self):
        """Return the number densities for each state."""
        if self._n_e is not None or self._n_elem is not None:
            return (self.n_elem * self.ionic_fractions).to(u.m ** -3)
        else:
            raise AtomicError(
                "Insufficient information to return number densities.")

    @number_densities.setter
    def number_densities(self, value):
        if self._n_elem or self._n_e:
            raise AtomicError
        if not isinstance(value, u.Quantity):
            raise AtomicError
        if np.any(value < 0):
            raise AtomicError("Number densities cannot be negative.")
        try:
            value = value.to(u.m ** -3)
        except (AttributeError, u.UnitsError):
            raise AtomicError

        self._n_elem = value.sum()
        self._ionic_fractions = number_densities / self._n_elem

    @property
    def T_e(self):
        """Return the electron temperature."""
        if not self._T_e:
            raise AtomicError("No electron temperature has been specified.")
        return self._T_e.to(u.K, equivalencies=u.temperature_energy())


    @T_e.setter
    def T_e(self, value):

        if value is None:
            self._T_e = None
        else:
            try:
                value = value.to(u.K, equivalencies=u.temperature_energy())
            except (AttributeError, u.UnitsError):
                raise AtomicError("Invalid temperature.")
            self._T_e = value

    @property
    def equil_ionic_fractions(self, T_e=None):
        """
        Return the equilibrium ionic fractions for temperature `T_e` or
        the temperature set in the IonizationState instance.  Not
        implemented.
        """
        raise NotImplementedError

    def equilibrate(self, T_e=None):
        """
        Set the ionic fractions to collisional ionization equilibrium
        for temperature `T_e`.  Not implemented.
        """
        self.ionic_fractions = self.equil_ionic_fractions

    @property
    def atomic_number(self) -> int:
        """Return the atomic number of the element."""
        return self._particle.atomic_number

    @property
    def element(self) -> str:
        """Return the atomic symbol of the element."""
        return self._particle.element

    @property
    def isotope(self) -> str:
        """
        Return the isotope symbol for an isotope, or `None` if the
        particle is not an isotope.
        """
        return self._particle.isotope

    @property
    def base_particle(self) -> str:
        """
        Return the element or isotope corresponding to this
        `~plasmapy.atomic.IonizationState` instance.
        """
        return self._particle.particle

    @property
    def particles(self) -> List[Particle]:
        """
        Return a list of the `~plasmapy.atomic.Particle` class
        instances.
        """
        return [Particle(self._particle.particle, Z=i) for i in range(self.atomic_number + 1)]

    @property
    def ionic_symbols(self) -> List[str]:
        """Return the ionic symbols for all charge states."""
        return [p.ionic_symbol for p in self.particles]

    @property
    def integer_charges(self) -> np.ndarray:
        """Return an array with the integer charges."""
        return np.arange(0, self.atomic_number + 1, dtype=np.int)

    def is_normalized(self, tol=None) -> bool:
        """
        Return `True` if the sum of the ionization fractions is equal to
        one within the allowed tolerance, and `False` otherwise.
        """

        tol = tol if tol is not None else self.tol

        if not isinstance(tol, (float, int, np.integer)):
            raise TypeError("tol must be an int or float.")

        if tol < 0 or tol > 1:
            raise ValueError("tol must be between 0 and 1, inclusive.")

        total = np.sum(self._ionic_fractions)

        return np.isclose(total, 1, atol=tol, rtol=0)

    def normalize(self):
        """
        Normalize the ionization state distribution so that the sum
        becomes equal to one.
        """
        self._ionic_fractions = self._ionic_fractions / np.sum(self._ionic_fractions)

    @property
    def tol(self) -> float:
        """Return the absolute tolerance for comparisons."""
        return self._tol

    @tol.setter
    def tol(self, atol):
        """Set the absolute tolerance for comparisons."""
        if not isinstance(atol, (float, int, np.integer)):
            raise TypeError("The attribute tol must be a float or an int.")
        if 0 <= atol <= 1.0:
            self._tol = atol
        else:
            raise ValueError("Need 0 <= tol <= 1.")




class IonizationStates:
    """
    Describe the ionization state distributions of multiple elements.

    Parameters
    ----------
    inputs: `list`, `tuple`, or `dict`
        A `list` or `tuple` of elements or isotopes (if `T_e` is
        provided); a `list` of `~plasmapy.atomic.IonizationState`
        instances; or a `dict` with elements or isotopes as keys and
        a `~numpy.ndarray` of ionic fractions as the values.

    T_e: `~astropy.units.Quantity`
        The electron temperature in units of temperature or thermal
        energy per particle.

    abundances: `dict` or `str`, optional
        The relative abundances of each element in the plasma.

    log_abundances: `dict`, optional
        The base 10 logarithm of the relative abundances of each element
        in the plasma.

    number_densities: `dict`, optional
        The number densities of elements (including both neutral atoms
        and ions) in units of inverse volume.

    n_H: ~astropy.units.Quantity, optional
        The number density of neutral and ionized hydrogen atoms.  May
        only be used if the ionization state of hydrogen is included.

    Raises
    ------
    AtomicError

    Notes
    -----
    Exactly one of `abundances`, `log_abundances`, and
    `number_densities` must be specified.



    """

    @check_quantity({
        "T_e": {"units": u.K, "none_shall_pass": True},
        "n_H": {"units": u.m ** -3, "none_shall_pass": True},
    })
    def __init__(
            self,
            inputs,
            *,
            T_e=None,
            abundances=None,
            log_abundances=None,
            n_H=None,
            tol=1e-15,
        ):

        self.ionic_fractions = inputs

        self._pars = collections.defaultdict(lambda: None)
        self.T_e = T_e

        self.abundances = abundances
        self.log_abundances = log_abundances


    @property
    def ionic_fractions(self):
        return self._ionic_fractions

    @ionic_fractions.setter
    def ionic_fractions(self, inputs: Union[Dict, List, Tuple]):
        """Set the ionic fractions"""
        if isinstance(inputs, dict):
            original_keys = inputs.keys()

            ionfrac_types = {type(inputs[key]) for key in original_keys}
            if u.Quantity in ionfrac_types and len(ionfrac_types != 1):
                raise TypeError(
                    "Ionic fraction information may only be inputted "
                    "as a Quantity object if all ionic fractions are "
                    "Quantity objects with units of inverse volume.")

            # Create a dictionary of Particle instances
            particles = dict()
            for key in original_keys:
                try:
                    particles[key] = key if isinstance(key, Particle) else Particle(key)
                except (InvalidParticleError, TypeError) as exc:
                    raise AtomicError(
                        f"Unable to create IonizationStates instance "
                        f"because {key} is not a valid particle") from exc

            # The particles whose ionization states are to be recorded
            # should be elements or isotopes but not ions or neutrals.
            is_element = particles[key].is_category('element')
            has_charge_info = particles[key].is_category(any_of=["charged", "uncharged"])
            if not is_element or has_charge_info:
                raise AtomicError(
                    f"{key} is not an element or isotope without "
                    f"charge information.")

            # Sort the original keys by atomic number (and if needed, by mass number)
            sorted_tuples = sorted([(
            particles[key].atomic_number,
            particles[key].mass_number if particles[key].isotope else 0,
            key) for key in original_keys])
            sorted_keys = [sorted_tuple[2] for sorted_tuple in sorted_tuples]

            self._ionic_fractions = {}
            for key in sorted_keys:
                new_key = particles[key].particle
                if isinstance(inputs[key], u.Quantity):
                    try:
                        number_densities = inputs[key].to(u.m ** -3)
                        n_elem = np.sum(number_densities)
                        self._ionic_fractions[new_key] = np.array(number_densities / n_elem)
                    except u.UnitConversionError as exc:
                        raise AtomicError("Units are not inverse volume.") from exc
                elif isinstance(inputs[key], np.ndarray) and inputs[key].dtype.kind == 'f':
                    self._ionic_fractions[particles[key].particle] = inputs[key]
                else:
                    try:
                        self._ionic_fractions[particles[key].particle] = \
                            np.array(inputs[key], dtype=np.float)
                    except ValueError:
                        raise AtomicError(f"Inappropriate ionic fractions for {key}.")

            print(self._ionic_fractions)

    def __getitem__(self, value):
        ...

    def __iter__(self):
        ...

    def __next__(self):
        ...

    def __eq__(self, other):
        ...

    @property
    def elements(self):
        return self._elements

    @property
    def abundances(self):
        return self._pars['abundances']

    @abundances.setter
    def abundances(self, value):
        self._pars['abundances'] = value

    @property
    def log_abundances(self):
        return np.log10(self.abundances)

    @log_abundances.setter
    def log_abundances(self, value):
        if value is None:
            self._pars['abundances'] = None
        else:
            # Add checks
            self._pars['abundances'] = 10 ** value

    @property
    def T_e(self):
        return self._pars['T_e']

    @T_e.setter
    def T_e(self, electron_temperature):
        if electron_temperature is None:
            self._pars['T_e'] = None
        else:
            try:
                self._pars['T_e'] = electron_temperature.to(u.K)
            except (AttributeError, u.UnitsError):
                raise AtomicError("Invalid electron temperature.")

    def equilibrate(self, T_e=None):
        if T_e is not None:
            raise NotImplementedError
        elif self._T_e is not None:
            raise NotImplementedError
        else:
            raise AtomicError

    @property
    def tol(self) -> float:
        """Return the absolute tolerance for comparisons."""
        return self._tol

    @tol.setter
    def tol(self, atol):
        """Set the absolute tolerance for comparisons."""
        if not isinstance(atol, (float, int, np.integer)):
            raise TypeError("The attribute tol must be a float or an int.")
        if 0 <= atol <= 1.0:
            self._tol = atol
        else:
            raise ValueError("Need 0 <= tol <= 1.")
