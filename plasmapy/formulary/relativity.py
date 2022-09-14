r"""Functions for calculating relativistic quantities (:math:`v \to c`)."""
__all__ = ["Lorentz_factor", "relativistic_energy", "RelativisticBody"]

import astropy.units as u
import numpy as np

from astropy.constants import c
from numbers import Integral, Real
from numpy.typing import DTypeLike
from typing import Dict, Optional, Union

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

    if not np.all((np.abs(V) <= c) | (np.isnan(V))):
        raise utils.RelativityError(
            "The Lorentz factor cannot be calculated for "
            "speeds faster than the speed of light."
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
    return γ * m * c**2


class RelativisticBody:
    r"""
    A physical body that is moving at a velocity relative to the speed
    of light.

    Parameters
    ----------
    particle : |ParticleLike|, |CustomParticle|, |ParticleList|, or |Quantity|
        A representation of a particle from which to get the mass
        of the relativistic body. If it is a |Quantity|, then it must
        have units of mass and describe the body's rest mass.

    V : |Quantity|, optional
        The velocity of the relativistic body in units convertible to
        m/s. The absolute magnitude of ``V`` cannot be greater than
        :math:`c`\ .

    momentum : |Quantity|, optional
        The momentum of the relativistic body in units convertible to
        kg·m/s.

    total_energy : |Quantity|, optional, |keyword-only|
       The sum of the mass energy and the kinetic energy in units
       convertible to joules. Must be non-negative.

    kinetic_energy : |Quantity|, optional, |keyword-only|
       The kinetic energy of the relativistic body in units convertible
       to joules. Must be non-negative.

    v_over_c : real number or |Quantity|, optional, |keyword-only|
       The ratio of the velocity to the speed of light. Must have an
       absolute magnitude :math:`≤ 1`\ .

    lorentz_factor : real number or |Quantity|, optional, |keyword-only|
       The Lorentz factor of the relativistic body. Must be
       :math:`≥ 1`\ .

    Z : integer, optional, |keyword-only|
        The charge number associated with ``particle``.

    mass_numb : integer, optional, |keyword-only|
        The mass number associated with ``particle``.

    dtype : |DTypeLike|, optional, |keyword-only|, default: `numpy.longdouble`
        The `numpy` data type to use to store the inputs.

    Notes
    -----
    At most one of ``V``, ``momentum``, ``total_energy``,
    ``kinetic_energy``, ``v_over_c``, and ``lorentz_factor`` must be
    provided.

    .. caution::

       For ultra-high-energy cosmic rays (UHECRs), the velocity may be
       within roundoff error of :math:`c` and :math:`\frac{V}{c}` may be
       within roundoff error of 1.

    Examples
    --------
    >>> import astropy.units as u
    >>> relativistic_proton = RelativisticBody("p+", total_energy = 1 * u.GeV)
    >>> relativistic_proton.particle
    Particle("p+")
    >>> relativistic_proton.velocity
    <Quantity 1.03697...e+08 m / s>
    >>> relativistic_proton.v_over_c
    0.3458980898746...
    >>> relativistic_proton.lorentz_factor
    1.0657889247888...
    >>> relativistic_proton.mass_energy.to("GeV")
    <Quantity 0.93827... GeV>
    >>> relativistic_proton.total_energy.to("GeV")
    <Quantity 1. GeV>
    >>> relativistic_proton.mass
    <Quantity 1.67262...e-27 kg>

    |RelativisticBody| also works with multiple particles and/or
    velocities.

    >>> particles = ["p+", "e-"]
    >>> velocities = [2e5, 2e8] * u.m / u.s
    >>> relativistic_particles = RelativisticBody(particles, velocities)
    >>> relativistic_particles.momentum
    <Quantity [3.345244...e-22, 2.445659...e-22] kg m / s>
    """

    @staticmethod
    def _get_speed_like_input(
        velocity_like_arguments: Dict[str, Union[u.Quantity, Real]]
    ):

        not_none_arguments = {
            key: value
            for key, value in velocity_like_arguments.items()
            if value is not None
        }

        if len(not_none_arguments) != 1:
            raise ValueError(
                "RelativisticBody can accept no more than one of the following "
                "arguments: V, v_over_c, momentum, total_energy, kinetic_energy, "
                "and lorentz_factor."
            )

        return not_none_arguments or {"velocity": np.nan * u.m / u.s}

    def _store_velocity_like_argument(
        self, speed_like_input: Dict[str, Union[u.Quantity, Real]]
    ):
        """
        Take the velocity-like argument and store it via the setter for
        the corresponding attribute.
        """
        name = list(speed_like_input.keys())[0]
        value = speed_like_input[name]
        if self._dtype:
            value = u.Quantity(value, dtype=self._dtype)
        setattr(self, name, value)

    @validate_quantities(
        V={"can_be_inf": False, "none_shall_pass": True},
        momentum={"can_be_inf": False, "none_shall_pass": True},
        total_energy={"can_be_negative": False, "none_shall_pass": True},
        kinetic_energy={"can_be_negative": False, "none_shall_pass": True},
    )
    def __init__(
        self,
        particle: Union[ParticleLike, u.Quantity],
        V: u.m / u.s = None,
        momentum: u.kg * u.m / u.s = None,
        *,
        total_energy: u.J = None,
        kinetic_energy: u.J = None,
        v_over_c: Optional[Real] = None,
        lorentz_factor: Optional[Real] = None,
        Z: Optional[Integral] = None,
        mass_numb: Optional[Integral] = None,
        dtype: Optional[DTypeLike] = np.longdouble,
    ):

        self._particle = _physical_particle_factory(particle, Z=Z, mass_numb=mass_numb)
        self._dtype = dtype

        velocity_like_inputs = {
            "velocity": V,
            "momentum": momentum,
            "total_energy": total_energy,
            "kinetic_energy": kinetic_energy,
            "v_over_c": v_over_c,
            "lorentz_factor": lorentz_factor,
        }

        speed_like_input = self._get_speed_like_input(velocity_like_inputs)
        self._store_velocity_like_argument(speed_like_input)

    def __repr__(self):
        return f"RelativisticBody({self.particle}, {self.velocity})"

    @property
    def particle(self) -> Union[CustomParticle, Particle, ParticleList]:
        """
        Representation of the particle(s).

        Returns
        -------
        |Particle|, |CustomParticle|, or |ParticleList|
        """
        return self._particle

    @property
    @validate_quantities
    def mass(self) -> u.kg:
        r"""
        The rest mass of the body, :math:`m_0`\ .

        Returns
        -------
        ~astropy.units.Quantity
        """
        return u.Quantity(self.particle.mass, dtype=self._dtype)

    @property
    @validate_quantities
    def mass_energy(self) -> u.J:
        r"""
        The rest mass energy of the body, :math:`m_0 c^2`\ .

        Returns
        -------
        ~astropy.units.Quantity
        """
        return self.mass * c**2

    @property
    @validate_quantities
    def total_energy(self) -> u.J:
        r"""
        The sum of the rest mass energy and the kinetic energy of the
        body,

        .. math::

           E_\mathrm{tot} ≡ γ m_0 c^2.

        Returns
        -------
        ~astropy.units.Quantity
        """
        return np.sqrt(self.momentum**2 * c**2 + self.mass_energy**2)

    @property
    @validate_quantities
    def kinetic_energy(self) -> u.J:
        """
        The kinetic energy of the body,

        .. math::

           E_K ≡ m_0 c^2 (γ-1).

        Returns
        -------
        ~astropy.units.Quantity
        """
        return self.total_energy - self.mass_energy

    @property
    @validate_quantities
    def v_over_c(self) -> Real:
        r"""
        The velocity of the body divided by the velocity of light:
        :math:`\frac{V}{c}`\ .

        Returns
        -------
        float
        """
        return (self.velocity / c).to(u.dimensionless_unscaled).value

    @property
    @validate_quantities
    def velocity(self) -> u.m / u.s:
        r"""
        The velocity of the body, :math:`V`\ .

        Returns
        -------
        ~astropy.units.Quantity
        """
        velocity = self.momentum / np.sqrt(self.mass**2 + self.momentum**2 / c**2)
        return velocity.to(u.m / u.s)

    @property
    @validate_quantities
    def lorentz_factor(self) -> Real:
        r"""
        The Lorentz factor of the body,

        .. math::

           γ ≡ \frac{1}{\sqrt{1 - \frac{V^2}{c^2}}}.

        Returns
        -------
        float
        """
        return Lorentz_factor(self.velocity)

    @property
    @validate_quantities
    def momentum(self) -> u.kg * u.m / u.s:
        r"""
        The magnitude of the momentum of the body,

        .. math::

           p ≡ γ m_0 V.

        Returns
        -------
        ~astropy.units.Quantity
        """
        return getattr(self, "_momentum")

    @kinetic_energy.setter
    @validate_quantities(E_K={"can_be_negative": False})
    def kinetic_energy(self, E_K: u.J):
        self.total_energy = E_K + self.mass_energy

    @total_energy.setter
    @validate_quantities(E_tot={"can_be_negative": False})
    def total_energy(self, E_tot: u.J):
        self._momentum = np.sqrt(E_tot**2 - self.mass_energy**2) / c

    @v_over_c.setter
    def v_over_c(self, v_over_c_: Real):
        self.velocity = v_over_c_ * c

    @velocity.setter
    @validate_quantities
    def velocity(self, V: u.m / u.s):
        self._momentum = (Lorentz_factor(V) * self.mass * V).to(u.kg * u.m / u.s)

    @lorentz_factor.setter
    def lorentz_factor(self, γ: Union[Real, u.Quantity]):
        if not isinstance(γ, (Real, u.Quantity)):
            raise TypeError("Invalid type for Lorentz factor")

        if isinstance(γ, u.Quantity):
            try:
                γ = γ.to(u.dimensionless_unscaled).value
            except u.UnitConversionError as exc:
                raise u.UnitConversionError(
                    "The Lorentz factor must be dimensionless."
                ) from exc

        if γ < 1:
            raise ValueError("The Lorentz factor must be ≥ 1")

        self.velocity = c * np.sqrt(1 - γ**-2)

    @momentum.setter
    @validate_quantities
    def momentum(self, p: u.kg * u.m / u.s):
        self._momentum = p.to(u.kg * u.m / u.s)

    def __eq__(self, other) -> bool:

        _attributes_to_compare = (
            "particle",
            "kinetic_energy",
            "mass_energy",
            "total_energy",
            "v_over_c",
            "momentum",
        )

        for attr in _attributes_to_compare:
            if not hasattr(other, attr):
                return False
            self_value = getattr(self, attr)
            other_value = getattr(other, attr)
            if self_value != other_value:
                return False
        return True
