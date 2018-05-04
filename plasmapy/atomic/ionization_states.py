import numpy as np
import astropy.units as u
import typing
import collections

from .particle_class import Particle
from .particle_input import particle_input
from ..utils import AtomicError, check_quantity

State = collections.namedtuple('State', ['Z', 'f', 'ionic_symbol'])


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
        The ionization fractions of an element.  This parameter must
        have a length equal to one more than the atomic number, and
        must sum to one within an absolute tolerance of `atol`.

    T_e: ~astropy.units.Quantity, optional
        The electron temperature or thermal energy per particle.  Being
        implemented.

    n_e: ~astropy.units.Quantity, optional
        The electron number density.  Being implemented.

    element_density: ~astropy.units.Quantity, optional
        The number density of ions and atoms.  Being implemented.

    number_densities: ~astropy.units.Quantity, optional
        The number densities of each charge state.  Being implemented.

    atol: float or int, optional
        The absolute tolerance used by `~numpy.isclose` when making
        comparisons and testing normalizations.  Defaults to `1e-15`.

    rtol: float or int, optional
        The relative tolerance used by `~numpy.isclose` when making
        comparisons and testing normalizations.  Defaults to `0`.

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
                 atol: typing.Union[float, int] = 1e-15,
                 rtol: typing.Union[float, int] = 0.0, ):
        """
        Initialize a `~plasmapy.atomic.IonizationState` instance.
        """

        self.atol = atol
        self.rtol = rtol

        self.n_e = n_e
        self.element_density = element_density
        self.number_densities = number_densities

        self.T_e = T_e

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
        self.charge_index = -1
        return self

    def __next__(self):
        if self.charge_index < self.atomic_number:
            self.charge_index += 1
            return State(self.charge_index, self._ionization_state[self.charge_index])
        else:
            del self.charge_index
            raise StopIteration

    def __eq__(self, other):
        if not isinstance(other, IonizationState):
            raise AtomicError(
                "Instances of the IonizationState class may only be "
                "compared with other IonizationState instances.")
        if self.element != other.element:
            raise AtomicError(
                "Only ionization states of the same element may be "
                "compared.")

        same_fractions = np.allclose(
            self.ionic_fractions,
            other.ionic_fractions,
            atol=self.atol,
            rtol=self.rtol,
        )

        return same_fractions

    @property
    def ionic_fractions(self) -> np.ndarray:
        """Return the ionic fractions."""
        return self._ionic_fractions

    @ionic_fractions.setter
    def ionic_fractions(self, value):
        """Set the ionic fractions."""
        try:
            if not isinstance(value, np.ndarray) or 'float' not in str(value.dtype):
                value = np.array(value, dtype=np.float)
        except Exception as exc:
            raise AtomicError(
                f"Unable to set ionic fractions of {self.element} to "
                f"{value}.") from exc
        total = np.sum(value)
        if not np.isclose(total, 1, atol=self.atol, rtol=self.rtol):
            raise AtomicError(
                f"The sum of the ionic fractions of {self.element} "
                f"equals {total}, which is not approximately one.")
        self._ionic_fractions = value

    @property
    def equilibrium_ionic_fractions(self, T_e=None):
        """
        Return the equilibrium ionic fractions for temperature `T_e` or
        the temperature set in the IonizationState instance.
        """
        raise NotImplementedError

    def equilibrate(self, T_e=None):
        """
        Set the ionic fractions to collisional ionization equilibrium
        for temperature `T_e`.
        """
        self.ionic_fractions = self.equilibrium_ionic_fractions

    @property
    def atomic_number(self) -> int:
        return self._particle.atomic_number

    @property
    def element(self) -> str:
        return self._particle.element

    @property
    def particles(self) -> typing.List[Particle]:
        """Return a list of the """
        return [Particle(self._particle.particle, Z=i) for i in range(self.atomic_number + 1)]

    @property
    def ionic_symbols(self) -> typing.List[str]:
        """
        Return the ionic symbols for all charge states.
        """
        return [p.ionic_symbol for p in self.particles]

    def is_normalized(self, atol=None, rtol=None) -> bool:
        """
        Return `True` if the sum of the ionization fractions is
        """
        atol = atol if atol is not None else self.atol
        rtol = rtol if rtol is not None else self.rtol

        if not isinstance(atol, (int, float)) or not isinstance(rtol, (int, float)):
            raise TypeError

        total = np.sum(self._ionic_fractions)

        return np.isclose(total, 1, atol=atol, rtol=rtol)

    def normalize(self):
        """
        Normalize the ionization state distribution so that the sum
        becomes equal to one.
        """
        self._ionic_fractions = self._ionic_fractions / np.sum(self._ionic_fractions)

    @property
    def atol(self) -> float:
        """Return the absolute tolerance for comparisons."""
        return self._atol

    @atol.setter
    def atol(self, value):
        """Set the absolute tolerance for comparisons."""
        if not isinstance(value, (float, int)):
            raise TypeError("The attribute atol must be a float or an int")
        if 0 <= value <= 1.0:
            self._atol = value
        else:
            raise ValueError("Need 0 <= atol <= 1.")

    @property
    def rtol(self) -> float:
        """Return the relative tolerance for comparisons."""
        return self._rtol

    @rtol.setter
    def rtol(self, value):
        """Set the relative tolerance for comparisons."""
        if not isinstance(value, (float, int)):
            raise TypeError("The attribute rtol must be a float or an int")
        if 0 <= value <= 1.0:
            self._rtol = value
        else:
            raise ValueError("Need 0 <= rtol <= 1.")


class IonizationStates:
    ...