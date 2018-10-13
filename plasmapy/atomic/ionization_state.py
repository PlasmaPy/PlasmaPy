import typing
import astropy.units as u
import collections
import numpy as np
import numbers
import warnings

from plasmapy.atomic import Particle, particle_input
from plasmapy.utils.exceptions import AtomicError, ChargeError, InvalidParticleError
from plasmapy.utils.checks import check_quantity

State = collections.namedtuple(
    'State', [
        'integer_charge',
        'ionic_fraction',
        'ionic_symbol',
        'number_density',
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
        an element or isotope, or an integer representing the atomic
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

    n_elem: ~astropy.units.Quantity, keyword-only, optional
        The number density of the element, including neutrals and all
        ions.

    tol: float or integer, keyword-only, optional
        The absolute tolerance used by `~numpy.isclose` when testing
        normalizations and making comparisons.  Defaults to `1e-15`.

    Raises
    ------
    ~plasmapy.utils.AtomicError
        If the ionic fractions are not normalized or contain invalid
        values, or if number density information is provided through
        both `ionic_fractions` and `n_elem`.

    ~plasmapy.utils.InvalidParticleError
        If the particle is invalid.

    Examples
    --------
    >>> states = IonizationState('H', [0.6, 0.4], n_elem=1*u.cm**-3, T_e=11000*u.K)
    >>> states.ionic_fractions[0]  # fraction of hydrogen that is neutral
    0.6
    >>> states.ionic_fractions[1]  # fraction of hydrogen that is ionized
    0.4
    >>> states.n_e  # electron number density
    <Quantity 400000. 1 / m3>
    >>> states.n_elem  # element number density
    <Quantity 1000000. 1 / m3>

    Notes
    -----
    Calculation of collisional ionization equilibrium has not yet been
    implemented.

    """

    # TODO: Allow this class to (optionally?) handle negatively charged
    # ions.  There are instances where singly negatively charged ions
    # are important in astrophysical plasmas, such as H- in the
    # atmospheres of relatively cool stars.  There may be some rare
    # situations where doubly negatively charged ions show up too,
    # though triply negatively charged ions are very unlikely.

    # TODO: Add in functionality to find equilibrium ionization states.

    @check_quantity(
        T_e={"units": u.K, "none_shall_pass": True},
        n_elem={"units": u.m ** -3, "none_shall_pass": True},
    )
    @particle_input(require='element', exclude='ion')
    def __init__(self,
                 particle: Particle,
                 ionic_fractions=None,
                 *,
                 T_e: u.K = np.nan * u.K,
                 n_elem: u.m ** -3 = np.nan * u.m ** -3,
                 tol: typing.Union[float, int] = 1e-15):
        """Initialize a `~plasmapy.atomic.IonizationState` instance."""

        self._particle = particle

        try:
            self.tol = tol
            self.T_e = T_e

            if not np.isnan(n_elem) and isinstance(ionic_fractions, u.Quantity) and \
                    ionic_fractions.si.unit == u.m ** -3:
                raise AtomicError("Cannot simultaneously provide number density "
                                  "through both n_elem and ionic_fractions.")

            self.n_elem = n_elem
            self.ionic_fractions = ionic_fractions

            if ionic_fractions is None and not np.isnan(self.T_e):
                warnings.warn(
                    "Collisional ionization equilibration has not yet "
                    "been implemented in IonizationState; cannot set "
                    "ionic fractions.")

        except Exception as exc:
            raise AtomicError(
                f"Unable to create IonizationState instance for "
                f"{particle.particle}.") from exc

    def __str__(self) -> str:
        symbol = self.isotope if self._particle.isotope else self.element
        return f"<IonizationState of {symbol}>"

    def __repr__(self):
        """Show diagnostic information of an IonizationState instance."""

        minimum_ionic_fraction_to_show = 1e-3

        output = f"IonizationState instance of {self.particle}"

        if np.any(np.isnan(self.ionic_fractions)):
            return output

        Z_mean = "{:.2f}".format(self.Z_mean)
        output += f'with Z_mean = {Z_mean}\n\n'

        for state in self:
            if state.ionic_fraction > minimum_ionic_fraction_to_show:

                symbol = state.ionic_symbol
                if state.integer_charge < 10:
                    symbol = symbol[:-2] + ' ' + symbol[-2:]

                fraction = "{:.3f}".format(state.ionic_fraction)

                output += (f'  {symbol}: {fraction}')

                if np.isfinite(self.n_elem):
                    value = "{:.2e}".format(state.number_density.si.value)
                    output += f"     n_i = {value} m**-3"

                output += '\n'

        if np.isfinite(self.T_e) or np.isfinite(self.n_elem):
            output += '\n'
        if np.isfinite(self.T_e):
            T_e = "{:.2e}".format(self.T_e.value) + " K"
            output += f'  T_e    = {T_e}\n'
        if np.isfinite(self.n_elem):
            n_elem = "{:.2e}".format(self.n_elem.value) + " m**-3"
            n_e = "{:.2e}".format(self.n_e.value) + " m**-3"
            output += f'  n_elem = {n_elem}\n'
            output += f'  n_e    = {n_e}'

        return output

    def __getitem__(self, value) -> State:
        """Return the ionic fraction(s)."""
        if isinstance(value, slice):
            raise TypeError("IonizationState instances cannot be sliced.")

        if isinstance(value, numbers.Integral) and 0 <= value <= self.atomic_number:
            result = State(
                value,
                self.ionic_fractions[value],
                self.ionic_symbols[value],
                self.number_densities[value],
            )
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
                result = State(
                    Z,
                    self.ionic_fractions[Z],
                    self.ionic_symbols[Z],
                    self.number_densities[Z],
                )
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
                self.number_densities[self._charge_index],
            )
            self._charge_index += 1
            return result
        else:
            del self._charge_index
            raise StopIteration

    def __eq__(self, other):
        """
        Return `True` if the ionic fractions, number density (if set),
        and electron temperature (if set) are all equal, and `False`
        otherwise.

        Raises
        ------
        TypeError
            If `other` is not an `~plasmapy.atomic.IonizationState`
            instance.

        AtomicError
            If `other` corresponds to a different element or isotope.

        Examples
        --------
        >>> IonizationState('H', [1, 0], tol=1e-6) == IonizationState('H', [1, 1e-6], tol=1e-6)
        True
        >>> IonizationState('H', [1, 0], tol=1e-8) == IonizationState('H', [1, 1e-6], tol=1e-5)
        False

        """
        if not isinstance(other, IonizationState):
            raise TypeError(
                "An instance of the IonizationState class may only be "
                "compared with another IonizationState instance.")

        same_element = self.element == other.element
        same_isotope = self.isotope == other.isotope

        if not same_element or not same_isotope:
            raise AtomicError(
                "An instance of the IonizationState class may only be "
                "compared with another IonizationState instance if "
                "both correspond to the same element and/or isotope.")

        # Use the tighter of the two tolerances. For thermodynamic
        # quantities, use it as a relative tolerance because the values
        # may substantially depart from order unity.

        min_tol = np.min([self.tol, other.tol])

        same_T_e = np.isnan(self.T_e) and np.isnan(other.T_e) or \
            u.allclose(self.T_e, other.T_e, rtol=min_tol*u.K, atol=0*u.K)

        same_n_elem = np.isnan(self.n_elem) and np.isnan(other.n_elem) or \
            u.allclose(self.n_elem, other.n_elem, rtol=min_tol*u.m**-3, atol=0*u.m**-3)

        same_fracs = np.allclose(self.ionic_fractions, other.ionic_fractions, rtol=0, atol=min_tol)

        return np.all([same_element, same_isotope, same_T_e, same_n_elem, same_fracs])

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
            return

        try:
            if np.min(fractions) < 0:
                raise AtomicError("Cannot have negative ionic fractions.")

            if len(fractions) != self.atomic_number + 1:
                raise AtomicError(
                    "The length of ionic_fractions must be "
                    f"{self.atomic_number + 1}.")

            elif isinstance(fractions, u.Quantity):
                fractions = fractions.to(u.m ** -3)
                self.n_elem = np.sum(fractions)
                self._ionic_fractions = np.array(fractions/self.n_elem)
            else:
                fractions = np.array(fractions, dtype=np.float64)
                if np.any(fractions < 0) or np.any(fractions > 1):
                    raise AtomicError("Ionic fractions must be between 0 and 1.")
                if not np.isclose(np.sum(fractions), 1, rtol=0, atol=self.tol):
                    raise AtomicError("Ionic fractions must sum to one.")
                self._ionic_fractions = fractions

        except Exception as exc:
            raise AtomicError(
                f"Unable to set ionic fractions of {self.element} "
                f"to {fractions}.") from exc

    @property
    @u.quantity_input
    def n_e(self) -> u.m ** -3:
        """
        Return the electron number density assuming a single species
        plasma.
        """
        return np.sum(self._n_elem * self.ionic_fractions * self.integer_charges)

    @property
    @u.quantity_input
    def n_elem(self) -> u.m ** -3:
        """
        Return the number density of atoms plus ions for this
        species.
        """
        return self._n_elem.to(u.m ** -3)

    @n_elem.setter
    @u.quantity_input
    def n_elem(self, value: u.m ** -3):
        """The number density of atoms plus ions of this species."""
        if value < 0 * u.m ** -3:
            raise AtomicError
        if 0 * u.m ** -3 < value <= np.inf * u.m ** -3:
            self._n_elem = value.to(u.m ** -3)
        elif np.isnan(value):
            self._n_elem = np.nan * u.m ** -3

    @property
    @u.quantity_input
    def number_densities(self) -> u.m ** -3:
        """Return the number densities for each state."""
        try:
            return (self.n_elem * self.ionic_fractions).to(u.m ** -3)
        except Exception:
            return np.full(self.atomic_number + 1, np.nan) * u.m ** -3

    @number_densities.setter
    @u.quantity_input
    def number_densities(self, value: u.m ** -3):
        """Set the number densities for each state."""
        if np.any(value.value < 0):
            raise AtomicError("Number densities cannot be negative.")
        if len(value) != self.atomic_number + 1:
            raise AtomicError(f"Incorrect number of charge states for {self.particle}")
        try:
            value = value.to(u.m ** -3)
        except (AttributeError, u.UnitsError):
            raise AtomicError

        self._n_elem = value.sum()
        self._ionic_fractions = value / self._n_elem

    @property
    @u.quantity_input(equivalencies=u.temperature_energy())
    def T_e(self) -> u.K:
        """Return the electron temperature."""
        if self._T_e is None:
            raise AtomicError("No electron temperature has been specified.")
        return self._T_e.to(u.K, equivalencies=u.temperature_energy())

    @T_e.setter
    def T_e(self, value: u.K):
        """Set the electron temperature."""
        try:
            value = value.to(u.K, equivalencies=u.temperature_energy())
        except (AttributeError, u.UnitsError, u.UnitConversionError):
            raise AtomicError("Invalid temperature.") from None
        else:
            if value < 0 * u.K:
                raise AtomicError("T_e cannot be negative.")
        self._T_e = value

    @property
    def equil_ionic_fractions(self, T_e: u.K = None):
        """
        Return the equilibrium ionic fractions for temperature `T_e` or
        the temperature set in the IonizationState instance.  Not
        implemented.
        """
        raise NotImplementedError

    def equilibrate(self, T_e: u.K = None):
        """
        Set the ionic fractions to collisional ionization equilibrium
        for temperature `T_e`.  Not implemented.
        """
        # self.ionic_fractions = self.equil_ionic_fractions
        raise NotImplementedError

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
        # Returning None when not an isotope makes it easier to do
        # things like: symbol = x.isotope if x.isotope else x.element
        return self._particle.isotope

    @property
    def particle(self) -> str:
        return self.isotope if self.isotope else self.element

    @property
    def particles(self) -> typing.List[Particle]:
        """
        Return a list of the `~plasmapy.atomic.Particle` class
        instances corresponding to each ion.
        """
        return [Particle(self._particle.particle, Z=i) for i in range(self.atomic_number + 1)]

    @property
    def ionic_symbols(self) -> typing.List[str]:
        """Return the ionic symbols for all charge states."""
        return [particle.ionic_symbol for particle in self.particles]

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
        """Return the root mean square integer charge."""
        return np.sqrt(np.sum(self.ionic_fractions * np.arange(self.atomic_number + 1) ** 2))

    @property
    def Z_mode(self) -> np.float64:
        return np.argmax

    def _is_normalized(self, tol: numbers.Integral = None) -> bool:
        """
        Return `True` if the sum of the ionization fractions is equal to
        one within the allowed tolerance, and `False` otherwise.
        """
        tol = tol if tol is not None else self.tol
        if not isinstance(tol, numbers.Real):
            raise TypeError("tol must be an int or float.")
        if not 0 <= tol < 1:
            raise ValueError("Need 0 <= tol < 1.")
        total = np.sum(self._ionic_fractions)
        return np.isclose(total, 1, atol=tol, rtol=0)

    def normalize(self) -> None:
        """
        Normalize the ionization state distribution so that the sum
        becomes equal to one.
        """
        self._ionic_fractions = self._ionic_fractions / np.sum(self._ionic_fractions)

    @property
    def tol(self) -> numbers.Real:
        """Return the absolute tolerance for comparisons."""
        return self._tol

    @tol.setter
    def tol(self, atol: numbers.Real):
        """Set the absolute tolerance for comparisons."""
        if not isinstance(atol, numbers.Real):
            raise TypeError("The attribute tol must be a real number.")
        if 0 <= atol < 1:
            self._tol = atol
        else:
            raise ValueError("Need 0 <= tol < 1.")
