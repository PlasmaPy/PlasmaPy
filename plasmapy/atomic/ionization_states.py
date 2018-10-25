"""
A class for storing ionization state data for multiple elements or
isotopes.
"""

from numbers import Real, Integral
from typing import Dict, List, Optional, Tuple, Union
import collections

import numpy as np
import astropy.units as u

from plasmapy.atomic import atomic_number, Particle, particle_symbol, IonizationState, State
from plasmapy.utils import AtomicError, ChargeError, InvalidParticleError, check_quantity

__all__ = ["IonizationStates"]


class IonizationStates:
    """
    Describe the ionization state distributions of multiple elements
    or isotopes.

    Parameters
    ----------
    inputs: `list`, `tuple`, or `dict`
        A `list` or `tuple` of elements or isotopes (if `T_e` is
        provided); a `list` of `~plasmapy.atomic.IonizationState`
        instances; a `dict` with elements or isotopes as keys and
        a `~numpy.ndarray` of ionic fractions as the values; or a `dict`
        with elements or isotopes as keys and `~astropy.units.Quantity`
        instances with units of number density.

    abundances: `dict` or `str`, optional, keyword-only
        The relative abundances of each element in the plasma.

    log_abundances: `dict`, optional, keyword-only
        The base 10 logarithm of the relative abundances of each element
        in the plasma.

    n: ~astropy.units.Quantity, optional, keyword-only
        The number density scaling factor.  The number density of an
        element will be the product of its abundance and ``n``.

    T_e: `~astropy.units.Quantity`, optional, keyword-only
        The electron temperature in units of temperature or thermal
        energy per base_particle.

    kappa: float, optional, keyword-only
        The value of kappa for a kappa distribution function.

    tol: float or integer, keyword-only, optional
        The absolute tolerance used by `~numpy.isclose` when testing
        normalizations and making comparisons.  Defaults to ``1e-15``.

    equilibrate: `bool`, optional, keyword-only
        Set the ionic fractions to the estimated collisional ionization
        equilibrium.  Not implemented.

    Raises
    ------
    AtomicError
        If `~plasmapy.atomic.IonizationStates` cannot be instantiated.

    Examples
    --------
    >>> from plasmapy.atomic import IonizationStates
    >>> solar_corona = IonizationStates(['H', 'He', 'Fe'])

    Notes
    -----
    No more than one of ``abundances``, ``log_abundances``, and
    ``number_densities`` may be specified.

    Collisional ionization equilibrium is based on atomic data that
    has relative errors of order 20%.

    """

    @check_quantity(
        T_e={"units": u.K},
        n={"units": u.m ** -3},
    )
    def __init__(
            self,
            inputs: Union[Dict[str, np.ndarray], List, Tuple],
            *,
            T_e: u.K = np.nan * u.K,
            equilibrate: Optional[bool] = None,
            abundances: Optional[Dict[str, Real]] = None,
            log_abundances: Optional[Dict[str, Real]] = None,
            n: u.m ** -3 = np.nan * u.m ** -3,
            tol: Real = 1e-15,
            kappa: Real = np.inf):
        """
        Initialize an `~plasmapy.atomic.IonizationStates`.
        """

        try:
            self._pars = collections.defaultdict(lambda: None)
            self.T_e = T_e
            self.n = n
            self.tol = tol
            self.ionic_fractions = inputs
            self.abundances = abundances
            self.log_abundances = log_abundances
            self.kappa = kappa
        except Exception as exc:
            raise AtomicError("Unable to create IonizationStates instance.") from exc

    def __str__(self) -> str:
        return f"<IonizationStates for: {', '.join(self.base_particles)}>"

    def __repr__(self) -> str:
        """Show diagnostic information."""
        output = []

        output.append(f"IonizationStates instance for: {', '.join(self.base_particles)}")

        # Get the ionic symbol with the corresponding ionic fraction and
        # number density (if available), but only for the most abundant
        # ionization levels for each element.

        for ionization_state in self:
            states_info = ionization_state._get_states_info(minimum_ionic_fraction=0.01)
            if len(states_info) > 0:
                output += states_info
                output[-1] += "\n"

        attributes = []
        if np.isfinite(self.T_e):
            attributes.append("   T_e = " + "{:.2e}".format(self.T_e.value) + " K")
        if np.isfinite(self.kappa):
            attributes.append(" kappa = " + "{:.3f}".format(self.kappa))
        if np.isfinite(self.n_e):
            attributes.append("   n_e = " + "{:.2e}".format(self.T_e.value) + " m ** -3")
        output += ["\n".join(attributes)]

        if len(output) > 1:
            output[0] += "\n"
            output_string = "\n".join(output)
        else:
            output_string = output[0]

        return output_string

    def __getitem__(self, *values) -> IonizationState:

        errmsg = f"Invalid indexing for IonizationStates instance: {values[0]}"

        one_input = not isinstance(values[0], tuple)
        two_inputs = len(values[0]) == 2

        if not one_input and not two_inputs:
            raise TypeError(errmsg)

        try:
            arg1 = values[0] if one_input else values[0][0]
            int_charge = None if one_input else values[0][1]
            particle = arg1 if arg1 in self.base_particles else particle_symbol(arg1)

            if int_charge is None:
                return IonizationState(
                    particle=particle,
                    ionic_fractions=self.ionic_fractions[particle],
                    T_e=self._pars["T_e"],
                    n_elem=np.sum(self.number_densities[particle]),
                    tol=self.tol,
                )
            else:
                if not isinstance(int_charge, Integral):
                    raise TypeError(f"{int_charge} is not a valid charge for {base_particle}.")
                elif not 0 <= int_charge <= atomic_number(particle):
                    raise ChargeError(f"{int_charge} is not a valid charge for {base_particle}.")
                return State(
                    integer_charge=int_charge,
                    ionic_fraction=self.ionic_fractions[particle][int_charge],
                    ionic_symbol=particle,
                    number_density=self.number_densities[particle][int_charge]
                )
        except Exception as exc:
            raise AtomicError(errmsg) from exc

    def __setitem__(self, key, value):
        if isinstance(value, dict):
            raise NotImplementedError("Dictionary assignment not implemented.")
        else:
            try:
                particle = particle_symbol(key)
                if particle not in self.base_particles:
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
                    f"{repr(value)}") from exc

    def __iter__(self):
        self._element_index = 0
        return self

    def __next__(self):
        if self._element_index < len(self.base_particles):
            particle = self.base_particles[self._element_index]
            result = IonizationState(
                particle,
                self.ionic_fractions[particle],
                T_e=self.T_e,
                n_elem=np.sum(self.number_densities[particle]),
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

        if self.base_particles != other.base_particles:
            raise AtomicError

        tol = np.min([self.tol, other.tol])

        for element in self.base_particles:

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
    @u.quantity_input
    def n_e(self) -> u.m ** -3:
        """
        Return the electron number density under the assumption of
        quasineutrality.
        """
        number_densities = self.number_densities
        n_e = 0.0 * u.m ** -3
        for elem in self.base_particles:
            atomic_numb = atomic_number(elem)
            number_of_ionization_states = atomic_numb + 1
            integer_charges = np.linspace(0, atomic_numb, number_of_ionization_states)
            n_e += np.sum(number_densities[elem] * integer_charges)
        return n_e

    @property
    @u.quantity_input
    def n(self) -> u.m ** -3:
        """Return the number density scaling factor."""
        if 'H' not in self.base_particles or self._pars['n'] is None:
            raise AtomicError("The number density of hydrogen is not ")
        return self._pars['n']

    @n.setter
    def n(self, n: u.m ** -3):
        """Set the number density scaling factor."""
        try:
            n = n.to(u.m ** -3)
        except u.UnitConversionError as exc:
            raise AtomicError("Units cannot be converted to u.m ** -3.") from exc
        except Exception as exc:
            raise AtomicError(f"{n} is not a valid number density") from exc

        if n < 0 * u.m ** -3:
            raise AtomicError("Number density cannot be negative.")

        self._pars['n'] = n.to(u.m ** -3)

    @property
    def number_densities(self) -> Dict:
        """
        Return a `dict` containing the number densities for element or
        isotope.
        """
        return {
            elem: self.n * self.abundances[elem] * self.ionic_fractions[elem]
            for elem in self.base_particles
        }

    @property
    def abundances(self) -> Optional[Dict]:
        """Return the elemental abundances."""
        return self._pars['abundances']

    @abundances.setter
    def abundances(self, abundances_dict: Optional[Dict]):
        """
        Set the elemental (or isotopic) abundances.  The elements and
        isotopes must be the same as or a superset of the elements whose
        ionization states are being tracked.
        """
        if abundances_dict is None:
            self._pars['abundances'] = {
                elem: np.full(Particle(elem).atomic_number + 1, np.nan)
                for elem in self.base_particles
            }
        elif not isinstance(abundances_dict, dict):
            raise TypeError(
                f"The abundances attribute must be a dict with "
                f"elements or isotopes as keys and real numbers "
                f"representing relative abundances as values.")
        else:
            old_keys = abundances_dict.keys()
            try:
                new_keys_dict = {particle_symbol(old_key): old_key for old_key in old_keys}
            except Exception:
                raise AtomicError(
                    "The key {repr(old_key)} in the abundances "
                    "dictionary is not a valid element or isotope.")

            new_elements = new_keys_dict.keys()

            old_elements_set = set(self.base_particles)
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
    def log_abundances(self) -> Optional[Dict[str, Real]]:
        """
        Return a `dict` with atomic or isotope symbols as keys and the
        base 10 logarithms of the relative abundances as the
        corresponding values.
        """
        if self._pars['abundances'] is not None:
            log_abundances_dict = {}
            for key in self.abundances.keys():
                log_abundances_dict[key] = np.log10(self.abundances[key])
            return log_abundances_dict
        else:
            raise AtomicError("No abundances are available.")

    @log_abundances.setter
    def log_abundances(self, value: Optional[Dict[str, Real]]):
        """
        Set the base 10 logarithm of the relative abundances.
        """
        if value is not None:
            try:
                new_abundances_input = {}
                for key in value.keys():
                    new_abundances_input[key] = 10 ** value[key]
                self.abundances = new_abundances_input
            except Exception:
                raise AtomicError("Invalid log_abundances.") from None

    @property
    def T_e(self) -> u.K:
        """Return the electron temperature."""
        return self._pars['T_e']

    @T_e.setter
    def T_e(self, electron_temperature: u.K):
        """Set the electron temperature."""
        try:
            temperature = electron_temperature.to(u.K, equivalencies=u.temperature_energy())
        except (AttributeError, u.UnitsError):
            raise AtomicError(
                f"{electron_temperature} is not a valid temperature.") from None
        if temperature < 0 * u.K:
            raise AtomicError("The electron temperature cannot be negative.")
        self._pars['T_e'] = temperature

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
        return self._pars['kappa']

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
        self._pars['kappa'] = np.real(value)

    @property
    def ionic_fractions(self) -> Dict[str, np.array]:
        """
        Return a `dict` containing the ionic fractions for each element
        and isotope.

        The keys of this `dict` are the symbols for each element or
        isotope.  The values will be `~numpy.ndarray` objects containing
        the ionic fractions for each ionization level corresponding to
        each element or isotope.

        """
        return self._ionic_fractions

    @ionic_fractions.setter
    def ionic_fractions(self, inputs: Union[Dict, List, Tuple]):
        """Set the ionic fractions."""
        if isinstance(inputs, dict):
            original_keys = inputs.keys()

            ionfrac_types = {type(inputs[key]) for key in original_keys}
            if u.Quantity in ionfrac_types and len(ionfrac_types) != 1:
                raise TypeError(
                    "Ionic fraction information may only be inputted "
                    "as a Quantity object if all ionic fractions are "
                    "Quantity arrays with units of inverse volume.")

            # Create a dictionary of Particle instances

            # TODO: clean up since Particle(Particle(x)) == Particle(x)

            particles = dict()
            for key in original_keys:
                try:
                    particles[key] = key if isinstance(key, Particle) else Particle(key)
                except (InvalidParticleError, TypeError) as exc:
                    raise AtomicError(
                        f"Unable to create IonizationStates instance "
                        f"because {key} is not a valid base_particle") from exc

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

            _elements_and_isotopes = []
            _particle_instances = []
            new_ionic_fractions = {}

            for key in sorted_keys:
                new_key = particles[key].particle
                _particle_instances.append(particles[key])
                if new_key in _elements_and_isotopes:
                    raise AtomicError("Repeated particles in IonizationStates.")

                nstates_input = len(inputs[key])
                nstates = particles[key].atomic_number + 1
                if nstates != nstates_input:
                    raise AtomicError(
                        f"The ionic fractions array for {key} must "
                        f"have a length of {nstates}.")

                _elements_and_isotopes.append(new_key)
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
                    except ValueError as exc:
                        raise AtomicError(f"Inappropriate ionic fractions for {key}.") from exc

            for key in _elements_and_isotopes:
                fractions = new_ionic_fractions[key]
                if not np.all(np.isnan(fractions)):
                    if np.min(fractions) < 0 or np.max(fractions) > 1:
                        raise AtomicError(f"Ionic fractions for {key} are not between 0 and 1.")
                    if not np.isclose(np.sum(fractions), 1, atol=self.tol, rtol=0):
                        raise AtomicError(f"Ionic fractions for {key} are not normalized to 1.")

        elif isinstance(inputs, (list, tuple)):

            try:
                _particle_instances = [Particle(particle) for particle in inputs]
            except (InvalidParticleError, TypeError) as exc:
                raise AtomicError("Invalid inputs to IonizationStates") from exc

            _particle_instances.sort(
                key=lambda p: (p.atomic_number, p.mass_number if p.isotope else 0)
            )
            _elements_and_isotopes = [particle.particle for particle in _particle_instances]
            new_ionic_fractions = {
                particle.particle: np.full(
                    particle.atomic_number + 1,
                    fill_value=np.nan,
                    dtype=np.float64
                ) for particle in _particle_instances
            }
        else:
            raise TypeError("Incorrect inputs to set ionic_fractions.")

        for i in range(1, len(_particle_instances)):
            if _particle_instances[i - 1].element == _particle_instances[i].element:
                if not _particle_instances[i - 1].isotope and _particle_instances[i].isotope:
                    raise AtomicError("Cannot have an element and isotopes of that element.")
            if _particle_instances[i - 1].atomic_number > _particle_instances[i].atomic_number:
                raise AtomicError("_particles has not been sorted.")

        self._particle_instances = _particle_instances
        self._base_particles = _elements_and_isotopes
        self._ionic_fractions = new_ionic_fractions

    def equilibrate(
            self,
            T_e: u.K = np.nan * u.K,
            particles: str = 'all',
            kappa: Real = np.inf):
        """
        Set the ionic fractions to collisional ionization equilibrium.
        Not implemented.

        The electron temperature used to calculate the new equilibrium
        ionic fractions will be the argument ``T_e`` to this method if
        given, and otherwise the attribute ``T_e`` if no electon
        temperature is provided to this method.

        Parameters
        ----------
        T_e: ~astropy.units.Quantity, optional
            The electron temperature.

        particles: `list`, `tuple`, or `str`, optional
            The elements and isotopes to be equilibrated.  If
            ``particles`` is ``'all'`` (default), then all
            elements and isotopes will be equilibrated.

        kappa: Real
            The value of kappa for a kappa distribution for electrons.

        """
        raise NotImplementedError

    @property
    def base_particles(self) -> List[str]:
        """
        Return a list of the elements and isotopes whose ionization
        states are being kept track of.
        """
        return self._base_particles

    @base_particles.setter
    def base_particles(self, particles):
        if hasattr(self, "_base_particles"):
            raise AtomicError(
                "Cannot change base particles once they have been set.")
        else:
            self._base_particles = particles

    @property
    def tol(self) -> np.real:
        """Return the absolute tolerance for comparisons."""
        return self._tol

    @tol.setter
    def tol(self, atol: Real):
        """
        Set the absolute tolerance for comparisons.
        """
        if not isinstance(atol, Real):
            raise TypeError("The attribute tol must be a real number.")
        if 0 <= atol <= 1.0:
            self._tol = np.real(atol)
        else:
            raise ValueError("Need 0 <= tol <= 1.")

    def normalize(self) -> None:
        """
        Normalize the ionic fractions so that the sum for each element
        equals one.
        """
        for particle in self.base_particles:
            tot = np.sum(self.ionic_fractions[particle])
            self.ionic_fractions[particle] = self.ionic_fractions[particle] / tot
