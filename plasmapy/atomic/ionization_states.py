import numpy as np
import astropy.units as u
import typing
import collections

from .particle_class import Particle
from .particle_input import particle_input
from ..utils import AtomicError, check_quantity

State = collections.namedtuple('State', ['integer_charge', 'ionic_fraction', 'ionic_symbol'])


class IonizationState:
    """
    Describe the ionization state distribution of a single element.

    Parameters
    ----------
    particle: str, int, or ~plasmapy.atomic.Particle
        A `str` or `~plasmapy.atomic.Particle` instance representing
        an element or isotope, or an `int` representing the atomic
        number of an element.

    ionic_fractions: ~numpy.ndarray, list, or tuple
        The ionization fractions of an element, where the index
        corresponds to the integer charge.  This parameter must have a
        length equal to one more than the atomic number, and must sum to
        one within an absolute tolerance of `tol`.

    T_e: ~astropy.units.Quantity, optional
        The electron temperature or thermal energy per particle.  Being
        implemented.

    n_e: ~astropy.units.Quantity, optional
        The electron number density.  Being implemented.

    element_density: ~astropy.units.Quantity, optional
        The number density of ions and atoms.  Being implemented.

    number_densities: ~astropy.units.Quantity, optional
        The number densities of each charge state.  Being implemented.

    tol: float or int, optional
        The absolute tolerance used by `~numpy.isclose` when testing
        normalizations and making comparisons.  Defaults to `1e-15`.

    Raises
    ------
    ~plasmapy.utils.AtomicError

    """

    @check_quantity({
        "T_e": {"units": u.K, "none_shall_pass": True},
        "n_e": {"units": u.m ** -3, "none_shall_pass": True},
        "element_density": {"units": u.m ** -3, "none_shall_pass": True},
        "number_densities": {"units": u.m ** -3, "none_shall_pass": True},
    })
    @particle_input(require='element', exclude='ion')
    def __init__(self,
                 particle: Particle,
                 ionic_fractions=None,
                 T_e=None,
                 n_e=None,
                 element_density=None,
                 number_densities=None,
                 tol: typing.Union[float, int] = 1e-15):
        """
        Initialize a `~plasmapy.atomic.IonizationState` instance.
        """

        self.tol = tol

        # Thermodynamic variables are not yet implemented.

        self.T_e = T_e
        self.n_e = n_e
        self.element_density = element_density
        self.number_densities = number_densities

        # Store the Particle class instance so that we can access the
        # attributes of that as needed.

        self._particle = particle

        if ionic_fractions is None and T_e is not None:
            self.equilibrate(T_e)  # Not implemented.

        self.ionic_fractions = ionic_fractions

        if not self.is_normalized:
            raise AtomicError

    def __getitem__(self, value):
        """Return the ionic fraction(s)."""
        if isinstance(value, slice):
            return self.ionic_fractions[value.start:value.stop:value.step]
        elif isinstance(value, int) and 0 <= value <= self.atomic_number:
            return self.ionic_fractions[value]

        if not isinstance(value, Particle):
            try:
                value = Particle(value)
            except InvalidParticleError as exc:
                raise InvalidParticleError from exc(
                    f"{value} is not a valid integer charge or particle.")

        if value.element == self.element:
            return self.ionic_fractions[value.integer_charge]

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
        >>> IonizationState('H', [1, 0], tol=1e-6) == \
                IonizationState('H', [1, 1e-6], tol=1e-6)
        True
        >>> IonizationState('H', [1, 0], tol=1e-8) == \
                IonizationState('H', [1, 1e-6], tol=1e-5)
        False

        """
        if not isinstance(other, IonizationState):
            raise AtomicError(
                "Instances of the IonizationState class may only be "
                "compared with other IonizationState instances.")

        if self.element != other.element:
            raise AtomicError(
                "Only ionization states of the same element may be "
                "compared.")

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
    def ionic_fractions(self, value):
        """
        Set the ionic fractions, while checking that the new values are
        valid and normalized to one.
        """
        try:
            if not isinstance(value, np.ndarray) or 'float' not in str(value.dtype):
                value = np.array(value, dtype=np.float)
        except Exception as exc:
            raise AtomicError(
                f"Unable to set ionic fractions of {self.element} to "
                f"{value}.") from exc

        if np.any(value < 0) or np.any(value > 1):
            raise AtomicError(f"The values of ")

        total = np.sum(value)
        if not np.isclose(total, 1, atol=self.tol, rtol=0):
            raise AtomicError(
                f"The sum of the ionic fractions of {self.element} "
                f"equals {total}, which is not approximately one.")

        self._ionic_fractions = value

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
        self.ionic_fractions = self.equilibrium_ionic_fractions

    @property
    def atomic_number(self) -> int:
        """Return the atomic number of the element."""
        return self._particle.atomic_number

    @property
    def element(self) -> str:
        """Return the atomic symbol of the element."""
        return self._particle.element

    @property
    def particles(self) -> typing.List[Particle]:
        """
        Return a list of the `~plasmapy.atomic.Particle` class
        instances.
        """
        return [Particle(self._particle.particle, Z=i) for i in range(self.atomic_number + 1)]

    @property
    def ionic_symbols(self) -> typing.List[str]:
        """
        Return the ionic symbols for all charge states.
        """
        return [p.ionic_symbol for p in self.particles]

    def is_normalized(self, tol=None) -> bool:
        """
        Return `True` if the sum of the ionization fractions is
        """

        tol = tol if tol is not None else self.tol

        if not isinstance(tol, (int, float)):
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
        if not isinstance(atol, (float, int)):
            raise TypeError("The attribute tol must be a float or an int.")
        if 0 <= atol <= 1.0:
            self._tol = atol
        else:
            raise ValueError("Need 0 <= tol <= 1.")
