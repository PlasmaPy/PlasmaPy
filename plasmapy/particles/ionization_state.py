"""
Objects for storing ionization state data for a single element or for
a single ionization level.
"""

__all__ = ["IonicLevel", "IonizationState"]

import numpy as np
import warnings

from astropy import units as u
from numbers import Integral, Real
from typing import List, NoReturn, Optional, Union

from plasmapy.particles.atomic import ionic_levels
from plasmapy.particles.decorators import particle_input
from plasmapy.particles.exceptions import (
    ChargeError,
    InvalidParticleError,
    ParticleError,
)
from plasmapy.particles.particle_class import CustomParticle, Particle
from plasmapy.particles.particle_collections import ParticleList
from plasmapy.utils.decorators import validate_quantities

_number_density_errmsg = (
    "Number densities must be Quantity objects with units of inverse volume."
)


class IonicLevel:
    """
    Representation of the ionic fraction for a single ion.

    Parameters
    ----------
    ion : |atom-like|
        The ion for the corresponding ionic fraction.

    ionic_fraction : real number, optional
        The fraction of an element or isotope that is at this ionization
        level. Must be between 0 and 1, inclusive.

    number_density : `~astropy.units.Quantity`, optional
        The number density of this ion.

    See Also
    --------
    IonizationState
    ~plasmapy.particles.ionization_state_collection.IonizationStateCollection

    Examples
    --------
    >>> alpha_fraction = IonicLevel("alpha", ionic_fraction=0.31)
    >>> alpha_fraction.ionic_symbol
    'He-4 2+'
    >>> alpha_fraction.charge_number
    2
    >>> alpha_fraction.ionic_fraction
    0.31
    """

    def __eq__(self, other):

        if not isinstance(other, IonicLevel):
            return False

        if self.ionic_symbol != other.ionic_symbol:
            return False

        ionic_fraction_within_rtol = u.isclose(
            self.ionic_fraction,
            other.ionic_fraction,
            rtol=1e-15,
            equal_nan=True,
        )

        number_density_within_rtol = u.isclose(
            self.number_density,
            other.number_density,
            rtol=1e-15,
            equal_nan=True,
        )

        return ionic_fraction_within_rtol and number_density_within_rtol

    @particle_input
    def __init__(
        self, ion: Particle, ionic_fraction=None, number_density=None, T_i=None
    ):
        try:
            self.ion = ion
            self.ionic_fraction = ionic_fraction
            self.number_density = number_density
            self.T_i = T_i
        except (ValueError, TypeError) as exc:
            raise ParticleError("Unable to create IonicLevel object") from exc

    def __repr__(self):
        return (
            f"IonicLevel({repr(self.ionic_symbol)}, "
            f"ionic_fraction={self.ionic_fraction})"
        )

    @property
    def ionic_symbol(self) -> str:
        """The symbol of the ion."""
        return self.ion.ionic_symbol

    @property
    def charge_number(self) -> Integral:
        """The charge number of the ion."""
        return self.ion.charge_number

    @property
    def ionic_fraction(self) -> Real:
        r"""
        The fraction of particles of an element that are at this
        ionization level.

        Notes
        -----
        An ionic fraction must be in the interval :math:`[0, 1]`.

        If no ionic fraction is specified, then this attribute will be
        assigned the value of |nan|.
        """
        return self._ionic_fraction

    @ionic_fraction.setter
    def ionic_fraction(self, ionfrac: Optional[Real]):
        if ionfrac is None or np.isnan(ionfrac):
            self._ionic_fraction = np.nan
        else:
            try:
                out_of_range = ionfrac < 0 or ionfrac > 1
            except TypeError as ex:
                raise TypeError(f"Invalid ionic fraction: {ionfrac}") from ex
            else:
                if out_of_range:
                    raise ValueError("The ionic fraction must be between 0 and 1.")
                else:
                    self._ionic_fraction = ionfrac

    @property
    def number_density(self) -> u.m**-3:
        """The number density of the ion."""
        return self._number_density

    @number_density.setter
    @validate_quantities(
        n={"can_be_negative": False, "can_be_inf": False, "none_shall_pass": True},
    )
    def number_density(self, n: u.m**-3):
        self._number_density = np.nan * u.m**-3 if n is None else n

    @property
    def T_i(self) -> u.K:
        """The ion temperature of this particular charge state."""
        return self._T_i

    @T_i.setter
    @validate_quantities(
        T={"can_be_negative": False, "can_be_inf": False, "none_shall_pass": True},
    )
    def T_i(self, T: u.K):
        self._T_i = np.nan * u.K if T is None else T


class IonizationState:
    """
    Representation of the ionization state distribution of a single
    element or isotope.

    Parameters
    ----------
    particle : |particle-like|
        A `str` or `~plasmapy.particles.particle_class.Particle` instance
        representing an element, isotope, or ion; or an integer representing
        the atomic number of an element.

    ionic_fractions : `~numpy.ndarray`, `list`, `tuple`, or `~astropy.units.Quantity`; optional
        The ionization fractions of an element, where the indices
        correspond to the charge number.  This argument should contain the
        atomic number plus one items, and must sum to one within an
        absolute tolerance of ``tol`` if dimensionless.  Alternatively,
        this argument may be a `~astropy.units.Quantity` that represents
        the number densities of each neutral/ion.  This argument cannot
        be specified when ``particle`` is an ion.

    T_e : `~astropy.units.Quantity`, |keyword-only|, optional
        The electron temperature or thermal energy per electron.

    n_elem : `~astropy.units.Quantity`, |keyword-only|, optional
        The number density of the element, including neutrals and all
        ions.

    tol : `float` or integer, |keyword-only|, optional, default: ``1e-15``
        The absolute tolerance used by `~numpy.isclose` and similar
        functions when testing normalizations and making comparisons.

    Raises
    ------
    `~plasmapy.particles.exceptions.ParticleError`
        If the ionic fractions are not normalized or contain invalid
        values, or if number density information is provided through
        both ``ionic_fractions`` and ``n_elem``.

    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the particle is invalid.

    See Also
    --------
    IonicLevel
    plasmapy.particles.ionization_state_collection.IonizationStateCollection

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

    If the input particle is an ion, then the ionization state for the
    corresponding element or isotope will be set to ``1.0`` for that
    ion.  For example, when the input particle is an alpha particle, the
    base particle will be He-4, and all He-4 particles will be set as
    doubly charged.

    >>> states = IonizationState('alpha')
    >>> states.base_particle
    'He-4'
    >>> states.ionic_fractions
    array([0., 0., 1.])
    """

    # TODO: Allow this class to handle negatively charged

    # TODO: Add in functionality to find equilibrium ionization states.

    @particle_input(require="element")
    @validate_quantities(
        T_e={"unit": u.K, "equivalencies": u.temperature_energy()},
        T_i={
            "unit": u.K,
            "equivalencies": u.temperature_energy(),
            "none_shall_pass": True,
        },
    )
    def __init__(
        self,
        particle: Particle,
        ionic_fractions=None,
        *,
        T_e: u.K = np.nan * u.K,
        T_i: u.K = None,
        kappa: Real = np.inf,
        n_elem: u.m**-3 = np.nan * u.m**-3,
        tol: Union[float, int] = 1e-15,
    ):
        self._number_of_particles = particle.atomic_number + 1

        if particle.is_ion or particle.is_category(require=("uncharged", "element")):
            if ionic_fractions is not None:
                raise ParticleError(
                    "The ionic fractions must not be specified when "
                    "the input particle to IonizationState is an ion."
                )

            ionic_fractions = np.zeros(self._number_of_particles)
            ionic_fractions[particle.charge_number] = 1.0
            particle = Particle(particle.isotope or particle.element)
        self._particle = particle

        try:
            self.tol = tol
            self.T_e = T_e
            self.T_i = T_i
            self.kappa = kappa
            self.n_elem = n_elem
            self.ionic_fractions = ionic_fractions

        except TypeError as exc:
            raise ParticleError(
                f"Unable to create IonizationState object for {particle.symbol}."
            ) from exc
        else:
            if (
                not np.isnan(n_elem)
                and isinstance(ionic_fractions, u.Quantity)
                and ionic_fractions.si.unit == u.m**-3
            ):
                raise ParticleError(
                    "Cannot simultaneously provide number density "
                    "through both n_elem and ionic_fractions."
                )

            if ionic_fractions is None and not np.isnan(self.T_e):
                warnings.warn(
                    "Collisional ionization equilibration has not yet "
                    "been implemented in IonizationState; cannot set "
                    "ionic fractions."
                )

    def __str__(self) -> str:
        return f"<IonizationState instance for {self.base_particle}>"

    def __repr__(self) -> str:
        return self.__str__()

    def __getitem__(self, value) -> List[IonicLevel]:
        """Return information for a single ionization level."""
        if isinstance(value, slice):
            return [
                IonicLevel(
                    ion=Particle(self.base_particle, Z=val),
                    ionic_fraction=self.ionic_fractions[val],
                    number_density=self.number_densities[val],
                    T_i=self.T_i[val],
                )
                for val in range(self._number_of_particles)[value]
            ]

        if isinstance(value, Integral) and 0 <= value <= self.atomic_number:
            result = IonicLevel(
                ion=Particle(self.base_particle, Z=value),
                ionic_fraction=self.ionic_fractions[value],
                number_density=self.number_densities[value],
                T_i=self.T_i[value],
            )
        else:
            if not isinstance(value, Particle):
                try:
                    value = Particle(value)
                except InvalidParticleError as exc:
                    raise InvalidParticleError(
                        f"{value} is not a valid charge number or particle."
                    ) from exc

            same_element = value.element == self.element
            same_isotope = value.isotope == self.isotope
            has_charge_info = value.is_category(any_of=["charged", "uncharged"])

            if same_element and same_isotope and has_charge_info:
                Z = value.charge_number
                result = IonicLevel(
                    ion=Particle(self.base_particle, Z=Z),
                    ionic_fraction=self.ionic_fractions[Z],
                    number_density=self.number_densities[Z],
                    T_i=self.T_i[Z],
                )
            elif not (same_element and same_isotope):
                raise ParticleError("Inconsistent element or isotope.")
            else:
                raise ChargeError("No charge number provided.")

        return result

    def __setitem__(self, key, value):
        raise NotImplementedError(
            "Item assignment of an IonizationState instance is not "
            "allowed because the ionic fractions for different "
            "ionization levels must be set simultaneously due to the "
            "normalization constraint."
        )

    def __iter__(self):
        yield from [self[i] for i in range(self.atomic_number + 1)]

    def __eq__(self, other):
        """
        Return `True` if the ionic fractions, number density scaling
        factor (if set), and electron temperature (if set) are all
        equal, and `False` otherwise.

        Raises
        ------
        `TypeError`
            If ``other`` is not an `~plasmapy.particles.ionization_state.IonizationState`
            instance.

        `ParticleError`
            If ``other`` corresponds to a different element or isotope.

        Examples
        --------
        >>> IonizationState('H', [1, 0], tol=1e-6) == IonizationState('H', [1, 1e-6], tol=1e-6)
        True
        >>> IonizationState('H', [1, 0], tol=1e-8) == IonizationState('H', [1, 1e-6], tol=1e-5)
        False

        """
        if not isinstance(other, IonizationState):
            return False

        same_element = self.element == other.element
        same_isotope = self.isotope == other.isotope

        if not same_element or not same_isotope:
            return False

        # Use the tighter of the two tolerances. For thermodynamic
        # quantities, use it as a relative tolerance because the values
        # may substantially depart from order unity.

        min_tol = np.min([self.tol, other.tol])

        same_T_e = (
            np.isnan(self.T_e)
            and np.isnan(other.T_e)
            or u.allclose(self.T_e, other.T_e, rtol=min_tol, atol=0 * u.K)
        )

        same_n_elem = (
            np.isnan(self.n_elem)
            and np.isnan(other.n_elem)
            or u.allclose(self.n_elem, other.n_elem, rtol=min_tol, atol=0 * u.m**-3)
        )

        # For the next line, recall that np.nan == np.nan is False

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
        The ionic fractions, where the index corresponds to the charge
        number.

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
                raise ParticleError(  # noqa: TC301
                    "Cannot have negative ionic fractions."
                )

            if len(fractions) != self.atomic_number + 1:
                raise ParticleError(  # noqa: TC301
                    "The length of ionic_fractions must be "
                    f"{self.atomic_number + 1}."
                )

            if isinstance(fractions, u.Quantity):
                fractions = fractions.to(u.m**-3)
                self.n_elem = np.sum(fractions)
                self._ionic_fractions = np.array(fractions / self.n_elem)
            else:
                fractions = np.array(fractions, dtype=np.float64)
                sum_of_fractions = np.sum(fractions)
                all_nans = np.all(np.isnan(fractions))

                if not all_nans:
                    if np.any(fractions < 0) or np.any(fractions > 1):
                        raise ParticleError("Ionic fractions must be between 0 and 1.")

                    if not np.isclose(sum_of_fractions, 1, rtol=0, atol=self.tol):
                        raise ParticleError("Ionic fractions must sum to one.")

                self._ionic_fractions = fractions

        except ParticleError as exc:
            raise ParticleError(
                f"Unable to set ionic fractions of {self.element} to {fractions}."
            ) from exc

    def _is_normalized(self, tol: Optional[Real] = None) -> bool:
        """
        `True` if the sum of the ionization fractions is equal to
        ``1`` within the allowed tolerance, and `False` otherwise.
        """
        tol = tol if tol is not None else self.tol
        if not isinstance(tol, Real):
            raise TypeError("tol must be an int or float.")
        if not 0 <= tol < 1:
            raise ValueError("Need 0 <= tol < 1.")
        total = np.sum(self._ionic_fractions)
        return np.isclose(total, 1, atol=tol, rtol=0)

    def normalize(self) -> NoReturn:
        """
        Normalize the ionization state distribution (if set) so that the
        sum of the ionic fractions becomes equal to one.

        This method may be used, for example, to correct for rounding
        errors.
        """
        self._ionic_fractions = self._ionic_fractions / np.sum(self._ionic_fractions)

    @property
    @validate_quantities
    def n_e(self) -> u.m**-3:
        """
        The electron number density assuming a single species plasma.
        """
        return np.sum(self._n_elem * self.ionic_fractions * self.charge_numbers)

    @property
    @validate_quantities
    def n_elem(self) -> u.m**-3:
        """The total number density of neutrals and all ions."""
        return self._n_elem.to(u.m**-3)

    @n_elem.setter
    @validate_quantities
    def n_elem(self, value: u.m**-3):
        """Set the number density of neutrals and all ions."""
        if value < 0 * u.m**-3:
            raise ParticleError
        if 0 * u.m**-3 < value <= np.inf * u.m**-3:
            self._n_elem = value.to(u.m**-3)
        elif np.isnan(value):
            self._n_elem = np.nan * u.m**-3

    @property
    @validate_quantities
    def number_densities(self) -> u.m**-3:
        """The number densities for each state."""
        try:
            return (self.n_elem * self.ionic_fractions).to(u.m**-3)
        except Exception:
            return np.full(self.atomic_number + 1, np.nan) * u.m**-3

    @number_densities.setter
    @validate_quantities
    def number_densities(self, value: u.m**-3):
        """Set the number densities for each state."""
        if np.any(value.value < 0):
            raise ParticleError("Number densities cannot be negative.")
        if len(value) != self.atomic_number + 1:
            raise ParticleError(
                f"Incorrect number of charge states for {self.base_particle}"
            )
        value = value.to(u.m**-3)

        self._n_elem = value.sum()
        self._ionic_fractions = value / self._n_elem

    @property
    def T_e(self) -> u.K:
        """The electron temperature."""
        if self._T_e is None:
            raise ParticleError("No electron temperature has been specified.")
        return self._T_e.to(u.K, equivalencies=u.temperature_energy())

    @T_e.setter
    @validate_quantities(value=dict(equivalencies=u.temperature_energy()))
    def T_e(self, value: u.K):
        """Set the electron temperature."""
        try:
            value = value.to(u.K, equivalencies=u.temperature_energy())
        except (AttributeError, u.UnitsError, u.UnitConversionError):
            raise ParticleError("Invalid temperature.") from None
        else:
            if value < 0 * u.K:
                raise ParticleError("T_e cannot be negative.")
        self._T_e = value

    @property
    @validate_quantities(
        validations_on_return=dict(
            equivalencies=u.temperature_energy(),
        )
    )
    def T_i(self) -> u.K:
        """
        The ion temperature. If the ion temperature has not been provided,
        then this attribute will provide the electron temperature.
        """
        return self._T_i

    @T_i.setter
    @validate_quantities(
        value=dict(
            equivalencies=u.temperature_energy(),
            none_shall_pass=True,
            can_be_negative=False,
        )
    )
    def T_i(self, value: u.K):
        """Set the ion temperature."""
        if value is None:
            self._T_i = np.repeat(self._T_e, self._number_of_particles)
            return

        if value.size == 1:
            self._T_i = np.repeat(value, self._number_of_particles)
        elif value.size == self._number_of_particles:
            self._T_i = value
        else:
            error_str = (
                "T_i must be set with either one common temperature"
                f" for all ions, or a set of {self._number_of_particles} of them. "
            )

            if value.size == 5 and self._number_of_particles != 5:
                error_str += f" For {self.base_particle}, five is right out."
            raise ParticleError(error_str)

    @property
    def kappa(self) -> Real:
        """
        The κ parameter for a kappa distribution function for electrons.

        The value of ``kappa`` must be greater than ``1.5`` in order to
        have a valid distribution function.  If ``kappa`` is
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
        """The atomic symbol of the element."""
        return self._particle.element

    @property
    def isotope(self) -> Optional[str]:
        """
        The isotope symbol for an isotope, or `None` if the particle is
        not an isotope.
        """
        return self._particle.isotope

    @property
    def base_particle(self) -> str:
        """The symbol of the element or isotope."""
        return self.isotope or self.element

    def to_list(self) -> ParticleList:
        """
        Return a `~plasmapy.particles.particle_collections.ParticleList`
        of the ionic levels.
        """
        return ionic_levels(self.base_particle)

    @property
    def atomic_number(self) -> int:
        """The atomic number of the element."""
        return self._particle.atomic_number

    def __len__(self):
        return self._number_of_particles

    @property
    def ionic_symbols(self) -> List[str]:
        """The ionic symbols for all charge states."""
        return self.to_list().symbols

    @property
    def charge_numbers(self) -> np.ndarray:
        """An array of the charge numbers."""
        return self.to_list().charge_number

    @property
    def Z_mean(self) -> np.float64:
        """Return the mean charge number."""
        if np.nan in self.ionic_fractions:
            raise ChargeError(
                "Z_mean cannot be found because no ionic fraction "
                f"information is available for {self.base_particle}."
            )
        return np.sum(self.ionic_fractions * self.charge_numbers)

    @property
    def Z_rms(self) -> np.float64:
        """The root-mean-square charge number."""
        return np.sqrt(np.sum(self.ionic_fractions * self.charge_numbers**2))

    @property
    def Z_most_abundant(self) -> List[Integral]:
        """
        A `list` of the charge numbers with the highest ionic fractions.

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
            raise ParticleError(
                f"Cannot find most abundant ion of {self.base_particle} "
                f"because the ionic fractions have not been defined."
            )

        return np.flatnonzero(
            self.ionic_fractions == self.ionic_fractions.max()
        ).tolist()

    @property
    def tol(self) -> Real:
        """
        The absolute tolerance for comparisons.

        This attribute is used as the ``atol`` parameter in
        `numpy.isclose`, `numpy.allclose`,
        `astropy.units.isclose`, and `astropy.units.allclose`
        when testing normalizations and making comparisons.
        """
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
        (if available) the number density and temperature for that ion.

        Parameters
        ----------
        minimum_ionic_fraction
            The minimum ionic fraction to return state information for.
        """

        states_info = []

        for state in self:
            if state.ionic_fraction >= minimum_ionic_fraction:
                state_info = ""
                symbol = state.ionic_symbol
                if state.charge_number < 10:
                    symbol = f"{symbol[:-2]} {symbol[-2:]}"
                fraction = f"{state.ionic_fraction:.3f}"

                state_info += f"{symbol}: {fraction}"

                if np.isfinite(self.n_elem):
                    value = f"{state.number_density.si.value:.2e}"
                    state_info += f"    n_i = {value} m**-3"

                if np.isfinite(state.T_i):
                    value = f"{state.T_i.si.value:.2e}"
                    state_info += f"    T_i = {value} K"

                states_info.append(state_info)

        return states_info

    def average_ion(
        self,
        *,
        include_neutrals: bool = True,
        use_rms_charge: bool = False,
        use_rms_mass: bool = False,
    ) -> CustomParticle:
        """
        Return a |CustomParticle| instance representing the average
        particle in this ionization state.

        By default, the weighted mean will be used as the average, with
        the ionic fractions as the weights. If ``use_rms_charge`` or
        ``use_rms_mass`` is `True`, then this method will return the root
        mean square of the charge or mass, respectively.

        Parameters
        ----------
        include_neutrals : `bool`, optional, |keyword-only|, default: `True`
            If `True`, include neutrals when calculating the mean values
            of the different particles.  If `False`, exclude neutrals.

        use_rms_charge : `bool`, optional, |keyword-only|, default: `False`
            If `True`, use the root-mean-square charge instead of the
            mean charge.

        use_rms_mass : `bool`, optional, |keyword-only|, default: `False`
            If `True`, use the root-mean-square mass instead of the mean
            mass.

        Returns
        -------
        ~plasmapy.particles.particle_class.CustomParticle

        Examples
        --------
        >>> state = IonizationState("He", [0.1, 0.9, 0.0])
        >>> state.average_ion()
        CustomParticle(mass=6.645657...e-27 kg, charge=1.44...e-19 C)
        >>> state.average_ion(include_neutrals=False)
        CustomParticle(mass=6.6455660...e-27 kg, charge=1.602...e-19 C)
        >>> state.average_ion(use_rms_charge=True, use_rms_mass=True)
        CustomParticle(mass=6.645657...e-27 kg, charge=1.519958...e-19 C)
        """
        min_charge = 0 if include_neutrals else 1

        particle_list = self.to_list()[min_charge:]
        abundances = self.ionic_fractions[min_charge:]

        return particle_list.average_particle(
            abundances=abundances,
            use_rms_charge=use_rms_charge,
            use_rms_mass=use_rms_mass,
        )

    def summarize(self, minimum_ionic_fraction: Real = 0.01) -> None:
        """
        Print quicklook information.

        Parameters
        ----------
        minimum_ionic_fraction : real number, default: ``0.01``
            If the ionic fraction for a particular ionization state is
            below this level, then information for it will not be
            printed.

        Examples
        --------
        >>> He_states = IonizationState(
        ...     'He',
        ...     [0.941, 0.058, 0.001],
        ...     T_e = 5.34 * u.K,
        ...     kappa = 4.05,
        ...     n_elem = 5.51e19 * u.m ** -3,
        ... )
        >>> He_states.summarize()
        IonizationState instance for He with Z_mean = 0.06
        ----------------------------------------------------------------
        He  0+: 0.941    n_i = 5.18e+19 m**-3    T_i = 5.34e+00 K
        He  1+: 0.058    n_i = 3.20e+18 m**-3    T_i = 5.34e+00 K
        ----------------------------------------------------------------
        n_elem = 5.51e+19 m**-3
        n_e = 3.31e+18 m**-3
        T_e = 5.34e+00 K
        kappa = 4.05
        ----------------------------------------------------------------

        """
        separator_line = [64 * "-"]

        output = [
            f"IonizationState instance for {self.base_particle} with Z_mean = {self.Z_mean:.2f}"
        ]
        attributes = []

        if not np.all(np.isnan(self.ionic_fractions)):
            output += separator_line
            output += self._get_states_info(minimum_ionic_fraction)
            output += separator_line
            # TODO add T_i somewhere around here, probably

        if not np.isnan(self.n_elem):
            attributes.extend(
                [
                    f"n_elem = {self.n_elem.value:.2e} m**-3",
                    f"n_e = {self.n_e.value:.2e} m**-3",
                ]
            )

        if not np.isnan(self.T_e):
            attributes.append(f"T_e = {self.T_e.value:.2f} K")
        if np.isfinite(self.kappa):
            attributes.append(f"kappa = {self.kappa:.2f}")

        if attributes:
            attributes += separator_line
            output += attributes

        for line in output:
            print(line)
