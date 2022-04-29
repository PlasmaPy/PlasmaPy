r"""Functions for calculating relativistic quantities (:math:`v \to c`)."""
__all__ = ["Lorentz_factor", "relativistic_energy", "RelativisticBody"]

import astropy.units as u
import numpy as np

from astropy.constants import c
from numbers import Integral, Real
from typing import Union

from plasmapy import utils
from plasmapy.particles._factory import _physical_particle_factory
from plasmapy.particles.particle_class import CustomParticle, Particle, ParticleLike
from plasmapy.particles.particle_collections import ParticleList
from plasmapy.utils.decorators import validate_quantities


@validate_quantities(V={"can_be_negative": True})
def Lorentz_factor(V: u.m / u.s):
    r"""
    Return the Lorentz factor.

    Parameters
    ----------
    V : `~astropy.units.Quantity`
        The velocity in units convertible to meters per second.

    Returns
    -------
    gamma : `float` or `~numpy.ndarray`
        The Lorentz factor associated with the inputted velocities.

    Raises
    ------
    `TypeError`
        If ``V`` is not a `~astropy.units.Quantity` and cannot be
        converted into a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the ``V`` is not in appropriate units.

    `ValueError`
        If the magnitude of ``V`` is faster than the speed of light.

    Warns
    -----
    `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    The Lorentz factor is a dimensionless number given by

    .. math::
        γ = \frac{1}{\sqrt{1-\frac{V^2}{c^2}}}

    The Lorentz factor is approximately one for sub-relativistic
    velocities, and :math:`γ → ∞` as the velocity approaches the
    speed of light.

    Examples
    --------
    >>> from astropy import units as u
    >>> velocity = 1.4e8 * u.m / u.s
    >>> Lorentz_factor(velocity)
    1.130885603948959
    >>> Lorentz_factor(299792458 * u.m / u.s)
    inf
    """

    if not np.all(np.abs(V) <= c):
        raise utils.RelativityError(
            "The Lorentz factor cannot be calculated for "
            "speeds faster than the speed of light. "
        )

    if V.size > 1:

        γ = np.zeros_like(V.value)

        equals_c = np.abs(V) == c
        is_slow = ~equals_c

        γ[is_slow] = ((1 - (V[is_slow] / c) ** 2) ** -0.5).value
        γ[equals_c] = np.inf

    else:
        γ = np.inf if np.abs(V) == c else ((1 - (V / c) ** 2) ** -0.5).value
    return γ


@validate_quantities(
    m={"can_be_negative": False}, validations_on_return={"can_be_negative": False}
)
def relativistic_energy(m: u.kg, v: u.m / u.s) -> u.Joule:
    """
    Calculate the relativistic energy (in joules) of an object of mass
    ``m`` and velocity ``v``.

    .. math::

        E = γ m c^2

    where :math:`γ` is the `Lorentz_factor`.

    Parameters
    ----------
    m : `~astropy.units.Quantity`
        The mass in units convertible to kilograms.

    v : `~astropy.units.Quantity`
        The velocity in units convertible to meters per second.

    Returns
    -------
    `~astropy.units.Quantity`
        The relativistic energy (in joules) of an object of mass ``m``
        moving at velocity ``v``.

    Raises
    ------
    `TypeError`
        If input arguments are not instances `~astropy.units.Quantity` or
        convertible to a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the ``v`` is not in appropriate units.

    `ValueError`
        If the magnitude of ``m`` is negative or arguments are complex.

    :exc:`~plasmapy.utils.exceptions.RelativityError`
        If the velocity ``v`` is greater than the speed of light.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Examples
    --------
    >>> from astropy import units as u
    >>> velocity = 1.4e8 * u.m / u.s
    >>> mass = 1 * u.kg
    >>> relativistic_energy(mass, velocity)
    <Quantity 1.01638929e+17 J>
    >>> relativistic_energy(mass, 299792458*u.m / u.s)
    <Quantity inf J>
    >>> relativistic_energy(1 * u.mg, 1.4e8 * u.m / u.s)
    <Quantity 1.01638929e+11 J>
    >>> relativistic_energy(-mass, velocity)
    Traceback (most recent call last):
        ...
    ValueError: The argument 'm' to function relativistic_energy() can not contain negative numbers.
    """
    γ = Lorentz_factor(v)
    return γ * m * c ** 2


class RelativisticBody:
    """
    A physical object that is moving a speed.

    Parameters
    ----------
    particle : |ParticleLike|, optional

    quantity : |Quantity|
        A |Quantity| with units of speed, energy,

    Z : integer, optional

    mass_numb : integer, optional
    """

    @validate_quantities
    def __init__(
        self,
        particle=None,
        *,
        kinetic_energy: u.J = None,
        Z: Integral = None,
        mass_numb: Integral = None,
    ):

        self._data = {
            "particle": _physical_particle_factory(particle, Z=Z, mass_numb=mass_numb)
        }
        self.kinetic_energy = kinetic_energy

    @property
    def particle(self) -> Union[CustomParticle, Particle, ParticleList]:
        print(self._data)
        return self._data["particle"]

    @particle.setter
    def particle(self, particle: ParticleLike):
        self._data["particle"] = _physical_particle_factory(particle)

    @property
    def rest_mass(self) -> u.Quantity[u.kg]:
        return self.particle.mass

    @property
    def mass_energy(self) -> u.Quantity[u.J]:
        """The rest mass energy, :math:`m_0 c^2`."""
        return self.rest_mass * c ** 2

    @property
    def total_energy(self) -> u.Quantity[u.J]:
        """The sum of the rest mass energy and the kinetic energy."""
        return self._data["total_energy"]

    @property
    def kinetic_energy(self) -> u.Quantity[u.J]:
        return self._data["kinetic_energy"]

    @kinetic_energy.setter
    def kinetic_energy(self, E: u.J):
        self._data["kinetic_energy"] = E
        self._data["total_energy"] = E + self.mass_energy
        self._data["lorentz_factor"] = (
            (self.total_energy / self.mass_energy).to("").value
        )
        self._data["v_over_c"] = np.sqrt(1.0 - 1.0 / self.lorentz_factor ** 2)
        self._data["speed"] = self.v_over_c * c

    @property
    def v_over_c(self) -> float:
        return self._data["v_over_c"]

    @property
    def speed(self) -> u.Quantity[u.m / u.s]:
        return self._data["speed"]

    @property
    def lorentz_factor(self) -> float:
        return self._data["lorentz_factor"]
