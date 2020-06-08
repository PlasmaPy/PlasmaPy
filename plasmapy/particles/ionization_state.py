"""
Objects for storing ionization state data for a single element or for
a single ionization level.
"""
__all__ = ["IonizationState", "State"]

import collections
import warnings
from numbers import Integral, Real
from typing import List, Optional, Union

import numpy as np
from astropy import units as u

from plasmapy.particles.exceptions import AtomicError, ChargeError, InvalidParticleError
from plasmapy.particles.particle_class import Particle
from plasmapy.particles.particle_input import particle_input
from plasmapy.utils.decorators import validate_quantities

_number_density_errmsg = (
    "Number densities must be Quantity objects with units of inverse " "volume."
)

# TODO: Change `State` into a class with validations for all of the
# TODO: attributes.

State = collections.namedtuple(
    "State", ["integer_charge", "ionic_fraction", "ionic_symbol", "number_density"]
)


class IonizationState:
    """
    Representation of the ionization state distribution of a single
    element or isotope.

    Parameters
    ----------
    particle: str, integer, or ~plasmapy.particles.Particle
        A `str` or `~plasmapy.particles.Particle` instance representing
        an element or isotope, or an integer representing the atomic
        number of an element.

    ionic_fractions: ~numpy.ndarray, list, tuple, or ~astropy.units.Quantity; optional
        The ionization fractions of an element, where the indices
        correspond to integer charge.  This argument should contain the
        atomic number plus one items, and must sum to one within an
        absolute tolerance of ``tol`` if dimensionless.  Alternatively,
        this argument may be a `~astropy.units.Quantity` that represents
        the number densities of each neutral/ion.

    T_e: ~astropy.units.Quantity, keyword-only, optional
        The electron temperature or thermal energy per particle.

    n_elem: ~astropy.units.Quantity, keyword-only, optional
        The number density of the element, including neutrals and all
        ions.

    tol: float or integer, keyword-only, optional
        The absolute tolerance used by `~numpy.isclose` when testing
        normalizations and making comparisons.  Defaults to ``1e-15``.

    Raises
    ------
    ~plasmapy.utils.AtomicError
        If the ionic fractions are not normalized or contain invalid
        values, or if number density information is provided through
        both ``ionic_fractions`` and ``n_elem``.

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
    # TODO: ions.  There are instances where singly negatively charged
    # TODO: ions are important in astrophysical plasmas, such as H- in
    # TODO: the atmospheres of relatively cool stars.  There may be some
    # TODO: rare situations where doubly negatively charged ions show up
    # TODO: too, but triply negatively charged ions are very unlikely.

    # TODO: Add in functionality to find equilibrium ionization states.

    @validate_quantities(T_e={"equivalencies": u.temperature_energy()})
    @particle_input(require="element", exclude="ion")
    def __init__(
        self,
        particle: Particle,
        ionic_fractions=None,
        *,
        T_e: u.K = np.nan * u.K,
        kappa: Real = np.inf,
        n_elem: u.m ** -3 = np.nan * u.m ** -3,
        tol: Union[float, int] = 1e-15,
    ):
        """Initialize an `~plasmapy.particles.IonizationState` instance."""

        self._particle_instance = particle

        try:
            self.tol = tol
            self.T_e = T_e
            self.kappa = kappa

            if (
                not np.isnan(n_elem)
                and isinstance(ionic_fractions, u.Quantity)
                and ionic_fractions.si.unit == u.m ** -3
            ):
                raise AtomicError(
                    "Cannot simultaneously provide number density "
                    "through both n_elem and ionic_fractions."
                )

            self.n_elem = n_elem
            self.ionic_fractions = ionic_fractions

            if ionic_fractions is None and not np.isnan(self.T_e):
                warnings.warn(
                    "Collisional ionization equilibration has not yet "
                    "been implemented in IonizationState; cannot set "
                    "ionic fractions."
                )

        except Exception as exc:
            raise AtomicError(
                f"Unable to create IonizationState instance for "
                f"{particle.particle}."
            ) from exc

    def __str__(self) -> str:
        return f"<IonizationState instance for {self.base_particle}>"

    def __repr__(self) -> str:
        return self.__str__()

    def __getitem__(self, value) -> State:
        """Return information for a single ionization level."""
        if isinstance(value, slice):
            raise TypeError("IonizationState instances cannot be sliced.")

        if isinstance(value, Integral) and 0 <= value <= self.atomic_number:
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
                        f"{value} is not a valid integer charge or " f"particle."
                    ) from exc

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

    def __setitem__(self, key, value):
        raise NotImplementedError(
            "Item assignment of an IonizationState instance is not "
            "allowed because the ionic fractions for different "
            "ionization levels must be set simultaneously due to the "
            "normalization constraint."
        )

    def __iter__(self):
        """Initialize an instance prior to iteration."""
        self._charge_index = 0
        return self

    def __next__(self):
        """
        Return a `~plasmapy.particles.State` instance that contains
        information about a particular ionization level.
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
        Return `True` if the ionic fractions, number density scaling
        factor (if set), and electron temperature (if set) are all
        equal, and `False` otherwise.

        Raises
        ------
        TypeError
            If ``other`` is not an `~plasmapy.particles.IonizationState`
            instance.

        AtomicError
            If ``other`` corresponds to a different element or isotope.

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
                "compared with another IonizationState instance."
            )

        same_element = self.element == other.element
        same_isotope = self.isotope == other.isotope

        if not same_element or not same_isotope:
            raise AtomicError(
                "An instance of the IonizationState class may only be "
                "compared with another IonizationState instance if "
                "both correspond to the same element and/or isotope."
            )

        # Use the tighter of the two tolerances. For thermodynamic
        # quantities, use it as a relative tolerance because the values
        # may substantially depart from order unity.

        min_tol = np.min([self.tol, other.tol])

        same_T_e = (
            np.isnan(self.T_e)
            and np.isnan(other.T_e)
            or u.allclose(self.T_e, other.T_e, rtol=min_tol * u.K, atol=0 * u.K)
        )

        same_n_elem = (
            np.isnan(self.n_elem)
            and np.isnan(other.n_elem)
            or u.allclose(
                self.n_elem, other.n_elem, rtol=min_tol * u.m ** -3, atol=0 * u.m ** -3
            )
        )

        # For the next line, recall that np.nan == np.nan is False (sigh)

        same_fractions = np.any(
            [
                np.allclose(
                    self.ionic_fractions, other.ionic_fractions, rtol=0, atol=min_tol
                ),
                np.all(np.isnan(self.ionic_fractions))
                and np.all(np.isnan(other.ionic_fractions)),
            ]
        )

        return np.all(
            [same_element, same_isotope, same_T_e, same_n_elem, same_fractions]
        )

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
        if fractions is None or np.all(np.isnan(fractions)):
            self._ionic_fractions = np.full(
                self.atomic_number + 1, np.nan, dtype=np.float64
            )
            return

        try:
            if np.min(fractions) < 0:
                raise AtomicError("Cannot have negative ionic fractions.")

            if len(fractions) != self.atomic_number + 1:
                raise AtomicError(
                    "The length of ionic_fractions must be "
                    f"{self.atomic_number + 1}."
                )

            if isinstance(fractions, u.Quantity):
                fractions = fractions.to(u.m ** -3)
                self.n_elem = np.sum(fractions)
                self._ionic_fractions = np.array(fractions / self.n_elem)
            else:
                fractions = np.array(fractions, dtype=np.float64)
                sum_of_fractions = np.sum(fractions)
                all_nans = np.all(np.isnan(fractions))

                if not all_nans:
                    if np.any(fractions < 0) or np.any(fractions > 1):
                        raise AtomicError("Ionic fractions must be between 0 and 1.")

                    if not np.isclose(sum_of_fractions, 1, rtol=0, atol=self.tol):
                        raise AtomicError("Ionic fractions must sum to one.")

                self._ionic_fractions = fractions

        except Exception as exc:
            raise AtomicError(
                f"Unable to set ionic fractions of {self.element} " f"to {fractions}."
            ) from exc

    def _is_normalized(self, tol: Optional[Real] = None) -> bool:
        """
        Return `True` if the sum of the ionization fractions is equal to
        one within the allowed tolerance, and `False` otherwise.
        """
        tol = tol if tol is not None else self.tol
        if not isinstance(tol, Real):
            raise TypeError("tol must be an int or float.")
        if not 0 <= tol < 1:
            raise ValueError("Need 0 <= tol < 1.")
        total = np.sum(self._ionic_fractions)
        return np.isclose(total, 1, atol=tol, rtol=0)

    def normalize(self) -> None:
        """
        Normalize the ionization state distribution (if set) so that the
        sum becomes equal to one.
        """
        self._ionic_fractions = self._ionic_fractions / np.sum(self._ionic_fractions)

    @property
    def equil_ionic_fractions(self, T_e: u.K = None):
        """
        Return the equilibrium ionic fractions for temperature ``T_e``
        or the temperature set in the IonizationState instance.  Not
        implemented.
        """
        raise NotImplementedError

    @u.quantity_input(equivalencies=u.temperature_energy())
    def equilibrate(self, T_e: u.K = np.nan * u.K):
        """
        Set the ionic fractions to collisional ionization equilibrium
        for temperature ``T_e``.  Not implemented.
        """
        # self.ionic_fractions = self.equil_ionic_fractions
        raise NotImplementedError

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
        """Return the total number density of neutrals and all ions."""
        return self._n_elem.to(u.m ** -3)

    @n_elem.setter
    @u.quantity_input
    def n_elem(self, value: u.m ** -3):
        """Set the number density of neutrals and all ions."""
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
            raise AtomicError(
                f"Incorrect number of charge states for " f"{self.base_particle}"
            )
        value = value.to(u.m ** -3)

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
    @u.quantity_input(equivalencies=u.temperature_energy())
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
    def kappa(self) -> np.real:
        """
        Return the kappa parameter for a kappa distribution function
        for electrons.

        The value of ``kappa`` must be greater than ``1.5`` in order to
        have a valid distribution function.  If ``kappa`` equals
        `~numpy.inf`, then the distribution function reduces to a
        Maxwellian.

        """
        return self._kappa

    @kappa.setter
    def kappa(self, value: Real):
        """
        Set the kappa parameter for a kappa distribution function for
        electrons.  The value must be between ``1.5`` and `~numpy.inf`.
        """
        kappa_errmsg = "kappa must be a real number greater than 1.5"
        if not isinstance(value, Real):
            raise TypeError(kappa_errmsg)
        if value <= 1.5:
            raise ValueError(kappa_errmsg)
        self._kappa = np.real(value)

    @property
    def element(self) -> str:
        """Return the atomic symbol of the element."""
        return self._particle_instance.element

    @property
    def isotope(self) -> Optional[str]:
        """
        Return the isotope symbol for an isotope, or `None` if the
        particle is not an isotope.
        """
        return self._particle_instance.isotope

    @property
    def base_particle(self) -> str:
        """Return the symbol of the element or isotope."""
        return self.isotope if self.isotope else self.element

    @property
    def atomic_number(self) -> int:
        """Return the atomic number of the element."""
        return self._particle_instance.atomic_number

    @property
    def _particle_instances(self) -> List[Particle]:
        """
        Return a list of the `~plasmapy.particles.Particle` class
        instances corresponding to each ion.
        """
        return [
            Particle(self._particle_instance.particle, Z=i)
            for i in range(self.atomic_number + 1)
        ]

    @property
    def ionic_symbols(self) -> List[str]:
        """Return the ionic symbols for all charge states."""
        return [particle.ionic_symbol for particle in self._particle_instances]

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
                f"information is available for {self.base_particle}."
            )
        return np.sum(self.ionic_fractions * np.arange(self.atomic_number + 1))

    @property
    def Z_rms(self) -> np.float64:
        """Return the root mean square integer charge."""
        return np.sqrt(
            np.sum(self.ionic_fractions * np.arange(self.atomic_number + 1) ** 2)
        )

    @property
    def Z_most_abundant(self) -> List[Integral]:
        """
        Return a `list` of the integer charges with the highest ionic
        fractions.

        Examples
        --------
        >>> He = IonizationState('He', [0.2, 0.5, 0.3])
        >>> He.Z_most_abundant
        [1]
        >>> Li = IonizationState('Li', [0.4, 0.4, 0.2, 0.0])
        >>> Li.Z_most_abundant
        [0, 1]

        """
        if np.any(np.isnan(self.ionic_fractions)):
            raise AtomicError(
                f"Cannot find most abundant ion of {self.base_particle} "
                f"because the ionic fractions have not been defined."
            )

        return np.flatnonzero(
            self.ionic_fractions == self.ionic_fractions.max()
        ).tolist()

    @property
    def tol(self) -> Real:
        """Return the absolute tolerance for comparisons."""
        return self._tol

    @tol.setter
    def tol(self, atol: Real):
        """Set the absolute tolerance for comparisons."""
        if not isinstance(atol, Real):
            raise TypeError("The attribute tol must be a real number.")
        if 0 <= atol < 1:
            self._tol = atol
        else:
            raise ValueError("Need 0 <= tol < 1.")

    def _get_states_info(self, minimum_ionic_fraction=0.01) -> List[str]:
        """
        Return a `list` containing the ion symbol, ionic fraction, and
        (if available) the number density for that ion.
        """

        states_info = []

        for state in self:
            if state.ionic_fraction > minimum_ionic_fraction:
                state_info = ""
                symbol = state.ionic_symbol
                if state.integer_charge < 10:
                    symbol = symbol[:-2] + " " + symbol[-2:]
                fraction = "{:.3f}".format(state.ionic_fraction)

                state_info += f"{symbol}: {fraction}"

                if np.isfinite(self.n_elem):
                    value = "{:.2e}".format(state.number_density.si.value)
                    state_info += f"    n_i = {value} m**-3"

                states_info.append(state_info)

        return states_info

    def info(self, minimum_ionic_fraction: Real = 0.01) -> None:
        """
        Print quicklook information for an
        `~plasmapy.particles.IonizationState` instance.

        Parameters
        ----------
        minimum_ionic_fraction: Real
            If the ionic fraction for a particular ionization state is
            below this level, then information for it will not be
            printed.  Defaults to 0.01.

        Example
        -------
        >>> He_states = IonizationState(
        ...     'He',
        ...     [0.941, 0.058, 0.001],
        ...     T_e = 5.34 * u.K,
        ...     kappa = 4.05,
        ...     n_elem = 5.51e19 * u.m ** -3,
        ... )
        >>> He_states.info()
        IonizationState instance for He with Z_mean = 0.06
        ----------------------------------------------------------------
        He  0+: 0.941    n_i = 5.18e+19 m**-3
        He  1+: 0.058    n_i = 3.20e+18 m**-3
        ----------------------------------------------------------------
        n_elem = 5.51e+19 m**-3
        n_e = 3.31e+18 m**-3
        T_e = 5.34e+00 K
        kappa = 4.05
        ----------------------------------------------------------------

        """
        separator_line = [64 * "-"]

        scientific = "{:.2e}"
        floaty = "{:.2f}"

        n_elem = scientific.format(self.n_elem.value)
        n_e = scientific.format(self.n_e.value)
        T_e = scientific.format(self.T_e.value)
        kappa = floaty.format(self.kappa)
        Z_mean = floaty.format(self.Z_mean)

        output = [
            f"IonizationState instance for {self.base_particle} with Z_mean = {Z_mean}"
        ]
        attributes = []

        if not np.all(np.isnan(self.ionic_fractions)):
            output += separator_line
            output += self._get_states_info(minimum_ionic_fraction)
            output += separator_line

        if not np.isnan(self.n_elem):
            attributes.append(f"n_elem = {n_elem} m**-3")
            attributes.append(f"n_e = {n_e} m**-3")
        if not np.isnan(self.T_e):
            attributes.append(f"T_e = {T_e} K")
        if np.isfinite(self.kappa):
            attributes.append(f"kappa = {kappa}")

        if attributes:
            attributes += separator_line
            output += attributes

        for line in output:
            print(line)
