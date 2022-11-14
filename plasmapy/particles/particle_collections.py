"""Collections of particle objects."""

__all__ = ["ParticleList", "ParticleListLike"]

import astropy.units as u
import collections
import contextlib
import numpy as np

from typing import Callable, Iterable, List, Optional, Sequence, Union

from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.particles.particle_class import (
    CustomParticle,
    DimensionlessParticle,
    Particle,
    ParticleLike,
)


class ParticleList(collections.UserList):
    """
    A `list` like collection of |Particle| and |CustomParticle| objects.

    Parameters
    ----------
    particles : iterable of |particle-like|, optional
        An iterable that provides a sequence of |particle-like| objects.
        Objects that are not a |Particle| or |CustomParticle| will be
        cast into a |Particle| or |CustomParticle|. If not provided, an
        empty |ParticleList| will be created.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidParticleError`
        If an object supplied to |ParticleList| is not |particle-like|.

    TypeError
        If a `~plasmapy.particles.particle_class.DimensionlessParticle`
        is provided.

    See Also
    --------
    ~plasmapy.particles.particle_class.CustomParticle
    ~plasmapy.particles.particle_class.Particle

    Examples
    --------
    A |ParticleList| can be created by calling it with an iterable like
    a `list` or `tuple` that provides |particle-like| objects.

    >>> from plasmapy.particles import ParticleList
    >>> particle_list = ParticleList(["e-", "e+"])
    >>> particle_list[0]
    Particle("e-")

    Attributes such as
    `~plasmapy.particles.particle_collections.ParticleList.mass`
    and `~plasmapy.particles.particle_collections.ParticleList.charge`
    will return a `~astropy.units.Quantity` array containing the values
    of the corresponding attribute for each particle in the
    |ParticleList|.

    >>> particle_list.mass
    <Quantity [9.1093...e-31, 9.1093...e-31] kg>
    >>> particle_list.charge
    <Quantity [-1.60217663e-19,  1.60217663e-19] C>
    >>> particle_list.symbols
    ['e-', 'e+']

    |ParticleList| instances can also be created through addition and
    multiplication with |Particle|, |CustomParticle|, and |ParticleList|
    instances.

    >>> from plasmapy.particles import Particle, CustomParticle
    >>> import astropy.units as u
    >>> proton = Particle("p+")
    >>> custom_particle = CustomParticle(mass=1e-26*u.kg, charge=6e-19*u.C)
    >>> 2 * proton + custom_particle
    ParticleList(['p+', 'p+', 'CustomParticle(mass=1e-26 kg, charge=6e-19 C)'])

    These operations may also be performed using |particle-like| objects.

    >>> particle_list + "deuteron"
    ParticleList(['e-', 'e+', 'D 1+'])

    Normal `list` methods may also be used on |ParticleList| objects.
    When a |particle-like| object is appended to a |ParticleList|, that
    object will be cast into a |Particle| or |CustomParticle|.

    >>> noble_gases = ParticleList(["He", "Ar", "Kr", "Xe", "Rn"])
    >>> noble_gases.append("Og")
    >>> noble_gases[-1]
    Particle("Og")

    The ``>`` operator may be used with |Particle| and |ParticleList|
    instances to access the nuclear reaction energy.

    >>> reactants = ParticleList(["deuterium", "tritium"])
    >>> products = ParticleList(["alpha", "neutron"])
    >>> energy = reactants > products
    >>> energy.to("MeV")
    <Quantity 17.58925... MeV>
    """

    @staticmethod
    def _list_of_particles_and_custom_particles(
        particles: Optional[Iterable[ParticleLike]],
    ) -> List[Union[Particle, CustomParticle]]:  # TODO #687
        """
        Convert an iterable that provides |particle-like| objects into a
        `list` containing |Particle| and |CustomParticle| instances.
        """
        new_particles = []
        if particles is None:
            return new_particles
        for obj in particles:
            if isinstance(obj, (Particle, CustomParticle)):
                new_particles.append(obj)
            elif isinstance(obj, DimensionlessParticle):
                raise TypeError(
                    "ParticleList instances cannot include dimensionless particles."
                )
            else:
                try:
                    new_particles.append(Particle(obj))
                except (TypeError, InvalidParticleError) as exc:
                    raise InvalidParticleError(
                        f"The object {obj} supplied to ParticleList is not a "
                        f"particle-like object."
                    ) from exc

        return new_particles

    def __init__(self, particles: Optional[Iterable] = None):
        self._data = self._list_of_particles_and_custom_particles(particles)

    @staticmethod
    def _cast_other_as_particle_list(other):
        if isinstance(other, ParticleList):
            return other

        with contextlib.suppress(TypeError, InvalidParticleError):
            return ParticleList(other)

        try:
            return ParticleList([other])
        except (InvalidParticleError, TypeError) as ex:
            raise InvalidParticleError(
                f"Cannot cast {other} into a ParticleList"
            ) from ex

    def __add__(self, other):
        try:
            other_as_particle_list = self._cast_other_as_particle_list(other)
        except (TypeError, InvalidParticleError) as exc:
            raise InvalidParticleError(
                f"Cannot add {repr(other)} to a ParticleList."
            ) from exc
        return ParticleList(self.data + other_as_particle_list.data)

    def __radd__(self, other):
        other_as_particle_list = self._cast_other_as_particle_list(other)
        return other_as_particle_list.__add__(self)

    def __repr__(self):
        return f"ParticleList({repr(self.symbols)})"

    def __gt__(self, other):
        from plasmapy.particles.nuclear import nuclear_reaction_energy

        other_as_particle_list = self._cast_other_as_particle_list(other)
        return nuclear_reaction_energy(
            reactants=self.symbols, products=other_as_particle_list.symbols
        )

    def __str__(self):
        return self.__repr__()

    def _get_particle_attribute(self, attr, unit=None, default=None):
        """
        Get the values of a particular attribute from all of the particles.

        If a ``unit`` is provided, then this function will return a
        `~astropy.units.Quantity` array with that unit.
        """
        values = [getattr(particle, attr, default) for particle in self.data]
        if unit:
            values = u.Quantity(values)
        return values

    def append(self, particle: ParticleLike):
        """Append a particle to the end of the |ParticleList|."""
        # TODO: use particle_input when it works with CustomParticle and ParticleLike
        if not isinstance(particle, (Particle, CustomParticle)):
            particle = Particle(particle)
        self.data.append(particle)

    @property
    def charge(self) -> u.C:
        """
        The electric charges of the particles.

        Returns
        -------
        `~astropy.units.Quantity`
        """
        return self._get_particle_attribute("charge", unit=u.C, default=np.nan * u.C)

    @property
    def data(self) -> List[Union[Particle, CustomParticle]]:
        """
        A `list` containing the particles contained in the
        |ParticleList| instance.

        .. important::

           This attribute should not be modified directly.

        Returns
        -------
        `list` of |Particle| or |CustomParticle|
        """
        return self._data

    def extend(self, iterable: Iterable[ParticleLike]):
        """
        Extend the sequence by appending |particle-like| elements from
        ``iterable``.

        See Also
        --------
        list.extend
        """
        if isinstance(iterable, ParticleList):
            self.data.extend(iterable)
        else:
            for obj in iterable:
                self.append(obj)

    @property
    def half_life(self) -> u.s:
        """
        The half-lives of the particles.

        Returns
        -------
        `~astropy.units.Quantity`
        """
        return self._get_particle_attribute("half_life", unit=u.s, default=np.nan * u.s)

    def insert(self, index, particle: ParticleLike):
        """Insert a particle before an index."""
        # TODO: use particle_input when it works with CustomParticle and ParticleLike
        if not isinstance(particle, (Particle, CustomParticle)):
            particle = Particle(particle)
        self.data.insert(index, particle)

    def is_category(
        self,
        *category_tuple,
        require: Union[str, Iterable[str]] = None,
        any_of: Union[str, Iterable[str]] = None,
        exclude: Union[str, Iterable[str]] = None,
    ) -> List[bool]:
        """
        Determine element-wise if the particles in the |ParticleList|
        meet categorization criteria.

        Return a `list` in which each element will be `True` if the
        corresponding particle is consistent with the categorization
        criteria, and `False` otherwise.

        Please refer to the documentation of
        `~plasmapy.particles.particle_class.Particle.is_category`
        for information on the parameters and categories, as well as
        more extensive examples.

        Returns
        -------
        `list` of `bool`

        Examples
        --------
        >>> particles = ParticleList(["proton", "electron", "tau neutrino"])
        >>> particles.is_category("lepton")
        [False, True, True]
        >>> particles.is_category(require="lepton", exclude="neutrino")
        [False, True, False]
        >>> particles.is_category(any_of=["lepton", "charged"])
        [True, True, True]
        """
        return [
            particle.is_category(
                *category_tuple,
                require=require,
                any_of=any_of,
                exclude=exclude,
            )
            for particle in self
        ]

    @property
    def charge_number(self) -> np.array:
        """
        The charges of the particles in units of the elementary
        charge.

        Returns
        -------
        |ndarray|
        """
        return np.array(self._get_particle_attribute("charge_number", default=np.nan))

    @property
    def mass(self) -> u.kg:
        """
        The masses of the particles.

        Returns
        -------
        `~astropy.units.Quantity`
        """
        return self._get_particle_attribute("mass", unit=u.kg, default=np.nan * u.J)

    @property
    def mass_energy(self) -> u.J:
        """
        The mass energies of the particles.

        If the particle is an isotope or nuclide, return the mass energy
        of the nucleus only.

        Returns
        -------
        `~astropy.units.Quantity`
        """
        return self._get_particle_attribute(
            "mass_energy",
            unit=u.J,
            default=np.nan * u.J,
        )

    def sort(self, key: Callable = None, reverse: bool = False):
        """
        Sort the |ParticleList| in-place.

        Parameters
        ----------
        key : callable
            A function that accepts one argument that is used to extract
            a comparison key for each item in the |ParticleList|.

        reverse : `bool`, default: `False`
            If `True`, the items in the |ParticleList| are sorted as if
            each comparison were reversed.

        Raises
        ------
        TypeError
            If ``key`` is not provided.

        See Also
        --------
        list.sort
        sorted

        Examples
        --------
        To sort a |ParticleList| by atomic number, set ``key`` to
        |atomic_number|.

        >>> from plasmapy.particles import ParticleList, atomic_number
        >>> elements = ParticleList(["Li", "H", "He"])
        >>> elements.sort(key=atomic_number)
        >>> print(elements)
        ParticleList(['H', 'He', 'Li'])

        We can also create a function to pass to ``key``. In this
        example, we sort first by atomic number and second by mass
        number using different attributes of |Particle|.

        >>> def sort_key(isotope):
        ...     return isotope.atomic_number, isotope.mass_number
        >>> isotopes = ParticleList(["He-3", "T", "H-1", "He-4", "D"])
        >>> isotopes.sort(key=sort_key)
        >>> print(isotopes)
        ParticleList(['H-1', 'D', 'T', 'He-3', 'He-4'])
        """
        if key is None:
            raise TypeError("Unable to sort a ParticleList without a key.")
        super().sort(key=key, reverse=reverse)

    @property
    def symbols(self) -> List[str]:
        """
        A `list` of the symbols of the particles.

        Returns
        -------
        `list` of `str`
        """
        return self._get_particle_attribute("symbol")

    def average_particle(
        self,
        abundances=None,
        *,
        use_rms_charge: bool = False,
        use_rms_mass: bool = False,
    ) -> Union[CustomParticle, Particle]:
        """
        Return a particle with the average mass and charge.

        By default, the mean will be used as the average. If the
        ``abundances`` are provided, then this method will return the
        weighted mean. If ``use_rms_charge`` or ``use_rms_mass`` is
        `True`, then this method will return the root-mean-square of the
        charge or mass, respectively. If all items in the |ParticleList|
        are equal to each other, then this method will return the first
        item in the |ParticleList|.

        Parameters
        ----------
        abundances : |array_like|, optional
            Real numbers representing relative abundances of the
            particles in the |ParticleList|. Must have the same number
            of elements as the |ParticleList|. This parameter gets
            passed to `numpy.average` via that function's ``weights``
            parameter. If not provided, the particles contained in the
            |ParticleList| are assumed to be equally abundant.

        use_rms_charge : `bool`, optional, |keyword-only|, default: `False`
            If `True`, use the root-mean-square charge instead of the
            mean charge.

        use_rms_mass : `bool`, optional, |keyword-only|, default: `False`
            If `True`, use the root-mean-square mass instead of the mean
            mass.

        Returns
        -------
        |Particle| or |CustomParticle|

        Examples
        --------
        >>> reactants = ParticleList(["electron", "positron"])
        >>> reactants.average_particle()
        CustomParticle(mass=9.109383...e-31 kg, charge=0.0 C)
        >>> reactants.average_particle(abundances=[1, 0.5])
        CustomParticle(mass=9.109383...e-31 kg, charge=-5.34058...e-20 C)
        >>> reactants.average_particle(use_rms_charge=True)
        CustomParticle(mass=9.109383...e-31 kg, charge=1.6021766...-19 C)
        >>> protons = ParticleList(["p+", "p+", "p+"])
        >>> protons.average_particle()
        Particle("p+")
        """
        # If all items in the ParticleList are the same, return that item.
        if len(set(self)) == 1:
            return self[0]

        def _average(array, weights, use_rms):
            if use_rms:
                return np.sqrt(np.average(array**2, weights=weights))
            else:
                return np.average(array, weights=weights)

        new_mass = _average(self.mass, weights=abundances, use_rms=use_rms_mass)
        new_charge = _average(self.charge, weights=abundances, use_rms=use_rms_charge)

        return CustomParticle(mass=new_mass, charge=new_charge)


# Override the docstrings for the parent class

ParticleList.clear.__doc__ = """Remove all items from the |ParticleList|."""

ParticleList.copy.__doc__ = """Return a shallow copy of the |ParticleList|."""

ParticleList.count.__doc__ = """
Return the number of occurrences of ``item``. Here, ``item`` must be
|particle-like|.
"""

ParticleList.extend.__doc__ = """
Extend |ParticleList| by casting |particle-like| items from
``iterable`` into particle objects.
"""

ParticleList.index.__doc__ = """
Return first index of a |particle-like| value.

Raises
------
`ValueError`
    If the value is not present.
"""

ParticleList.pop.__doc__ = """
Remove and return item at index (default last).

Raises
------
`IndexError`
    If the |ParticleList| is empty or the index is out of range.
"""

ParticleList.remove.__doc__ = """
Remove the first occurrence of a |particle-like| item.

Raises
------
`ValueError`
    If the value is not present.
"""

ParticleList.reverse.__doc__ = """Reverse the |ParticleList| in place."""

ParticleListLike = Union[ParticleList, Sequence[ParticleLike]]

ParticleListLike.__doc__ = r"""
An `object` is |particle-list-like| if it can be identified as a
|ParticleList| or cast into one.

When used as a type hint annotation, |ParticleListLike| indicates that
the corresponding argument should represent a sequence of physical
particles. Each item in a |ParticleListLike| object must be
|particle-like|.

Notes
-----
`~plasmapy.particles.particle_class.DimensionlessParticle` instances do
not uniquely represent a physical particle, and are thus not
|ParticleLike| and cannot be contained in a |ParticleListLike| object.

See Also
--------
~plasmapy.particles.particle_collections.ParticleList
~plasmapy.particles.particle_class.ParticleLike
~plasmapy.particles.decorators.particle_input

Examples
--------
Using |ParticleListLike| as a type hint annotation indicates that an
argument or variable should represent a sequence of |ParticleLike|
objects.

>>> from plasmapy.particles import ParticleList, ParticleListLike
>>> def contains_only_leptons(particles: ParticleListLike):
...     particle_list = ParticleList(particles)
...     return all(particle_list.is_category("lepton"))
>>> contains_only_leptons(["electron", "muon"])
True
"""
