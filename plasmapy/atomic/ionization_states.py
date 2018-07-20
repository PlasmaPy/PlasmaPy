"""Classes for storing ionization state data."""

from typing import Dict, List, Optional, Tuple, Union

import astropy.units as u
import collections
import numpy as np

from plasmapy.atomic.atomic import atomic_number
from plasmapy.atomic.particle_class import Particle
from plasmapy.atomic.particle_input import particle_input
from plasmapy.atomic.symbols import particle_symbol
from plasmapy.utils import (AtomicError, ChargeError, InvalidParticleError, check_quantity)

State = collections.namedtuple(
    'State', [
        'integer_charge',
        'ionic_fraction',
        'ionic_symbol',
    ])

_number_density_errmsg = (
    "Number densities must be Quantity objects with units of inverse "
    "volume."
)


class IonizationState:
    """
    Representation of the ionization state distribution of a single
    element or isotope.

    Parameters
    ----------
    particle: str, integer, or ~plasmapy.atomic.Particle
        A `str` or `~plasmapy.atomic.Particle` instance representing
        an element or isotope, or an `int` representing the atomic
        number of an element.

    ionic_fractions: ~numpy.ndarray, list, tuple, or ~astropy.units.Quantity
        The ionization fractions of an element, where the indices
        correspond to integer charge.  This argument should contain the
        atomic number plus one items, and must sum to one within an
        absolute tolerance of `tol` if dimensionless.  Alternatively,
        this argument may be a `~astropy.units.Quantity` that represents
        the number densities of each neutral/ion.

    T_e: ~astropy.units.Quantity, keyword-only, optional
        The electron temperature or thermal energy per particle.

    n_e: ~astropy.units.Quantity, keyword-only, optional
        The electron number density.

    n_elem: ~astropy.units.Quantity, keyword-only, optional
        The number density of the element, including neutrals and all
        ions.

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

    Examples
    --------
    >>> states = IonizationState('H', [0.6, 0.4], n_elem=1*u.cm**-3, T_e=11000*u.K)
    >>> states.ionic_fractions
    array([0.6, 0.4])
    >>> states.n_e  # electron number density
    <Quantity 400000. 1 / m3>
    >>> states.n_elem  # element number density
    <Quantity 1000000. 1 / m3>

    Notes
    -----
    Only one of `n_e` and `n_elem` may be set.

    """

    # TODO: Allow this class to (optionally?) handle negatively charged
    # ions.  There are instances where singly negatively charged ions
    # are important in astrophysical plasmas, such as H- in the
    # atmospheres of relatively cool stars.  There may be some rare
    # situations where doubly negatively charged ions show up too,
    # though triply negatively charged ions are very unlikely.

    # TODO: Add in functionality to find equilibrium ionization states.
    # How much data will this require?

    @check_quantity(
        T_e={"units": u.K, "none_shall_pass": True},
        n_e={"units": u.m ** -3, "none_shall_pass": True},
        n_elem={"units": u.m ** -3, "none_shall_pass": True},
    )
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

        self._particle = particle

        try:
            self.tol = tol
            self.T_e = T_e
            self.n_elem = n_elem
            self.n_e = n_e
            self.ionic_fractions = ionic_fractions

            # This functionality has not yet been implemented:

            # if self._ionic_fractions is None and self.T_e is not None:
            #     self.equilibrate()

        except Exception as exc:
            raise AtomicError(
                f"Unable to create IonizationState instance for "
                f"{particle.particle}.") from exc

    def __getitem__(self, value) -> State:
        """Return the ionic fraction(s)."""
        if isinstance(value, slice):
            raise TypeError("IonizationState instances cannot be sliced.")

        if isinstance(value, (int, np.integer)) and 0 <= value <= self.atomic_number:
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

        if fractions is None:
            self._ionic_fractions = np.full(self.atomic_number + 1, np.nan, dtype=np.float64)
        else:
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
                if self._n_e is not None or self._n_elem is not None:
                    raise AtomicError(
                        "The ionization state may be set using number "
                        "densities for each ion only if neither of the "
                        "electron density and element density has "
                        "already been set.")
                self.number_densities = fractions
            else:

                if not np.any(np.isnan(fractions)):

                    total = np.sum(fractions)
                    if not np.isclose(total, 1, atol=self.tol, rtol=0):
                        raise AtomicError(
                            f"The sum of the ionic fractions of {self.element} "
                            f"equals {total}, which is not approximately one.")
                    if not len(fractions) == self.atomic_number + 1:
                        raise AtomicError(
                            f"len(fractions) equals {len(fractions)}, but "
                            f"should equal {self.atomic_number + 1} which "
                            f"is the atomic number of {self.element} + 1."
                        )

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
                raise AtomicError(_number_density_errmsg) from None

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
        if self._n_elem is not None or self._n_e is not None:
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
        self._ionic_fractions = value / self._n_elem

    @property
    def T_e(self):
        """Return the electron temperature."""
        if self._T_e is None:
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
                raise AtomicError("Invalid temperature.") from None
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

    @property
    def Z_mean(self) -> np.float64:
        """Return the mean integer charge"""
        if np.nan in self.ionic_fractions:
            raise ChargeError(
                "Z_mean cannot be found because no ionic fraction "
                f"information is available for {self.particle}.")
        return np.sum(self.ionic_fractions * np.arange(self.atomic_number + 1))

    @property
    def Z_rms(self) -> np.float64:
        return np.sqrt(np.sum(self.ionic_fractions * np.arange(self.atomic_number + 1) ** 2))


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

    n: ~astropy.units.Quantity, optional
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

    @check_quantity(
        T_e={"units": u.K, "none_shall_pass": True},
        n={"units": u.m ** -3, "none_shall_pass": True},
    )
    def __init__(
            self,
            inputs,
            *,
            T_e=None,
            abundances=None,
            log_abundances=None,
            n=None,
            tol=1e-15,
        ):

        self._pars = collections.defaultdict(lambda: None)
        self.T_e = T_e
        self.n = n

        self.tol = tol

        if isinstance(inputs, dict):
            self.ionic_fractions = inputs
        elif isinstance(inputs, (list, tuple)):
            self.ionic_fractions = inputs
        else:
            raise TypeError(f"{inputs} are invalid inputs.")

        self.abundances = abundances
        self.log_abundances = log_abundances

    def __str__(self) -> str:
        join_str = ", " if len(self.elements) <= 5 else ","
        return f"<IonizationStates for: {join_str.join(self.elements)}>"

    def __repr__(self) -> str:
        return self.__str__()

    @property
    def ionic_fractions(self):
        try:
            return self._ionic_fractions
        except AttributeError as exc:
            raise AttributeError("Ionic fractions were not set.") from exc

    @ionic_fractions.setter
    def ionic_fractions(self, inputs: Union[Dict, List, Tuple]):
        """Set the ionic fractions"""
        if isinstance(inputs, dict):
            original_keys = inputs.keys()

            ionfrac_types = {type(inputs[key]) for key in original_keys}
            if u.Quantity in ionfrac_types and len(ionfrac_types) != 1:
                raise TypeError(
                    "Ionic fraction information may only be inputted "
                    "as a Quantity object if all ionic fractions are "
                    "Quantity arrays with units of inverse volume.")

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

            # We are sorting the elements/isotopes by atomic number and
            # mass number since we will often want to plot and analyze
            # things and this is the most sensible order.

            sorted_keys = sorted(original_keys, key=lambda k: (
                particles[k].atomic_number,
                particles[k].mass_number if particles[k].isotope else 0,
            ))

            _elements = []
            _particles = []
            new_ionic_fractions = {}
            for key in sorted_keys:
                new_key = particles[key].particle
                _particles.append(particles[key])
                if new_key in _elements:
                    raise AtomicError("Repeated particles in IonizationStates.")

                _elements.append(new_key)
                if isinstance(inputs[key], u.Quantity):
                    try:
                        number_densities = inputs[key].to(u.m ** -3)
                        n_elem = np.sum(number_densities)
                        new_ionic_fractions[new_key] = np.array(number_densities / n_elem)
                    except u.UnitConversionError as exc:
                        raise AtomicError("Units are not inverse volume.") from exc
                elif isinstance(inputs[key], np.ndarray) and inputs[key].dtype.kind == 'f':
                    new_ionic_fractions[particles[key].particle] = inputs[key]
                else:
                    try:
                        new_ionic_fractions[particles[key].particle] = \
                            np.array(inputs[key], dtype=np.float)
                    except ValueError:
                        raise AtomicError(f"Inappropriate ionic fractions for {key}.")

            for key in _elements:
                if np.min(new_ionic_fractions[key]) < 0 or np.max(new_ionic_fractions[key]) > 1:
                    raise AtomicError(f"Ionic fractions for {key} are not between 0 and 1.")
                if not np.isclose(np.sum(new_ionic_fractions[key]), 1, atol=self.tol, rtol=0):
                    raise AtomicError(f"Ionic fractions for {key} are not normalized to 1.")

        elif isinstance(inputs, (list, tuple)):

            try:
                _particles = [Particle(particle) for particle in inputs]
            except (InvalidParticleError, TypeError):
                raise AtomicError("Invalid inputs to IonizationStates")

            _particles.sort(key=lambda p: (p.atomic_number, p.mass_number if p.isotope else 0))
            _elements = [particle.particle for particle in _particles]
            new_ionic_fractions = {
                particle.particle: np.full(
                    particle.atomic_number + 1,
                    fill_value=np.nan,
                    dtype=np.float64
                ) for particle in _particles
            }
        else:
            raise TypeError

        # Because this depends on _particles being sorted, we add in an
        # easy check that atomic numbers do not decrease.
        for i in range(1, len(_particles)):
            if _particles[i - 1].element == _particles[i].element:
                if not _particles[i - 1].isotope and _particles[i].isotope:
                    raise AtomicError("Cannot have an element and isotopes of that element.")
            elif _particles[i - 1].atomic_number > _particles[i].atomic_number:
                raise AtomicError("_particles has not been sorted.")

        self._particles = _particles
        self._elements = _elements
        self._ionic_fractions = new_ionic_fractions

    def __getitem__(self, *values):

        errmsg = f"Invalid indexing for IonizationStates instance: {values[0]}"

        one_input = not isinstance(values[0], tuple)
        two_inputs = len(values[0]) == 2

        if not one_input and not two_inputs:
            raise TypeError(errmsg)

        try:
            arg1 = values[0] if one_input else values[0][0]
            int_charge = None if one_input else values[0][1]
            particle = arg1 if arg1 in self.elements else particle_symbol(arg1)

            if int_charge is None:
                return IonizationState(
                    particle=particle,
                    ionic_fractions=self.ionic_fractions[particle],
                    T_e=self._pars["T_e"],
                    tol=self.tol,
                )
            else:
                if not isinstance(int_charge, (int, np.integer)):
                    raise TypeError(f"{int_charge} is not a valid charge for {particle}.")
                elif not 0 <= int_charge <= atomic_number(particle):
                    raise ChargeError(f"{int_charge} is not a valid charge for {particle}.")
                return State(
                    integer_charge=int_charge,
                    ionic_fraction=self.ionic_fractions[particle][int_charge],
                    ionic_symbol=particle,
                )
        except Exception as exc:
            raise AtomicError(errmsg) from exc

    def __setitem__(self, key, value):
        if isinstance(value, dict):
            raise NotImplementedError("Dictionary assignment not implemented.")
        else:
            try:
                particle = particle_symbol(key)
                if particle not in self.elements:
                    raise AtomicError(
                        f"{key} is not one of the particles kept track of "
                        f"by this IonizationStates instance.")
                new_fractions = np.array(value, dtype=np.float64)
                if new_fractions.min() < 0 or new_fractions.max() > 1:
                    raise ValueError("Ionic fractions must be between 0 and 1.")
                if not np.isclose(np.sum(new_fractions), 1):
                    raise ValueError("Ionic fractions are not normalized.")
                if len(new_fractions) != atomic_number(particle) + 1:
                    raise ValueError(f"Incorrect size of ionic fraction array for {key}.")
                self._ionic_fractions[particle][:] = new_fractions[:]
            except Exception as exc:
                raise AtomicError(
                    f"Cannot set item for this IonizationStates "
                    f"instance for key = {repr(key)} and value = "
                    f"{repr(value)}")

    def __iter__(self):
        self._element_index = 0
        return self

    def __next__(self):
        if self._element_index < len(self.elements):
            particle = self.elements[self._element_index]
            result = IonizationState(
                particle,
                self.ionic_fractions[particle],
                T_e=self.T_e,
                tol=self.tol,
            )
            self._element_index += 1
            return result
        else:
            del self._element_index
            raise StopIteration

    def __eq__(self, other):
        """
        Test that the ionic fractions are approximately equal to
        another `~plasmapy.atomic.IonizationStates` instance.
        """
        # TODO: Should we check temperatures, densities, and other
        # parameters too?

        if not isinstance(other, IonizationStates):
            raise TypeError(
                "IonizationStates instance can only be compared with "
                "other IonizationStates instances.")

        if self.elements != other.elements:
            raise AtomicError

        tol = np.min([self.tol, other.tol])

        for element in self.elements:

            are_all_close = np.allclose(
                self.ionic_fractions[element],
                other.ionic_fractions[element],
                atol=tol,
                rtol=0,
            )

            if not are_all_close:
                return False

        return True

    @property
    def elements(self) -> List[str]:
        """
        Return a list of the elements whose ionization states are being
        kept track of.
        """
        return self._elements

    @elements.setter
    def elements(self, elems):
        if hasattr(self, "_elements"):
            raise AtomicError("Cannot change elements once they have been set.")
        else:
            self._elements = elems

    @property
    def abundances(self) -> Dict:
        if self._pars['abundances'] is None:
            raise AtomicError("No abundances are available.")
        return self._pars['abundances']

    @abundances.setter
    def abundances(self, abundances_dict: Optional[Dict]):
        """
        Set the elemental (or isotopic) abundances.  The elements and
        isotopes must be the same as or a superset of the elements whose
        ionization states are being tracked.
        """
        if abundances_dict is None:
            self._pars['abundances'] = None
        elif not isinstance(abundances_dict, dict):
            raise TypeError(
                f"The abundances argument {abundances_dict} must be a dict with elements "
                "or isotopes as keys and ")
        else:
            old_keys = abundances_dict.keys()
            try:
                new_keys_dict = {particle_symbol(old_key): old_key for old_key in old_keys}
            except Exception:
                raise AtomicError(
                    "The key {repr(old_key)} in the abundances "
                    "dictionary is not a valid element or isotope.")

            new_elements = new_keys_dict.keys()

            old_elements_set = set(self.elements)
            new_elements_set = set(new_elements)

            if old_elements_set > new_elements_set:
                raise AtomicError(
                    f"The abundances of the following particles are "
                    f"missing: {old_elements_set - new_elements_set}")

            new_abundances_dict = {}

            for element in new_elements:
                inputted_abundance = abundances_dict[new_keys_dict[element]]
                try:
                    inputted_abundance = float(inputted_abundance)
                except Exception:
                    raise TypeError

                if inputted_abundance < 0:
                    raise AtomicError(f"The abundance of {element} is negative.")
                new_abundances_dict[element] = inputted_abundance

            self._pars['abundances'] = new_abundances_dict

    @property
    def log_abundances(self):
        if self._pars['abundances'] is not None:
            log_abundances_dict = {}
            for key in self.abundances.keys():
                log_abundances_dict[key] = np.log10(self.abundances[key])
            return log_abundances_dict
        else:
            raise AtomicError("No abundances are available.")

    @log_abundances.setter
    def log_abundances(self, value):

        if value is not None:
            try:
                new_abundances_input = {}
                for key in value.keys():
                    new_abundances_input[key] = 10 ** value[key]
                self.abundances = new_abundances_input
            except Exception as exc:
                raise AtomicError("Invalid log_abundances.")

    @property
    def T_e(self):
        """Return the electron temperature."""
        return self._pars['T_e']

    @T_e.setter
    def T_e(self, electron_temperature):
        if electron_temperature is None:
            self._pars['T_e'] = None
        else:
            try:
                temp = electron_temperature.to(u.K, equivalencies=u.temperature_energy())
            except (AttributeError, u.UnitsError):
                raise AtomicError("Invalid electron temperature.")
            else:
                if temp < 0 * u.K:
                    raise AtomicError("The electron temperature cannot be negative.")
                self._pars['T_e'] = temp

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

    def normalize(self):
        for particle in self.elements:
            tot = np.sum(self.ionic_fractions[particle])
            self.ionic_fractions[particle] = self.ionic_fractions[particle] / tot

    @property
    def number_densities(self):
        ...

    @property
    def n_e(self):
        raise NotImplementedError

    @property
    def n(self):
        """
        The number density scaling factor, if defined.
        """
        if 'H' not in self.elements or self._pars['n'] is None:
            raise AtomicError("The number density of hydrogen is not ")
        return self._pars['n']

    @n.setter
    def n(self, n):
        if n is None:
            self._pars['n'] = n
        else:
            try:
                self._pars['n'] = n.to(u.m ** -3)
            except u.UnitConversionError:
                raise AtomicError("Units cannot be converted to u.m**-3.")
            except Exception:
                raise AtomicError(f"{n} is not a valid number density.") from None
