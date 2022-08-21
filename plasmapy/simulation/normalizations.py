import astropy.units as u
import numpy as np

from astropy.constants import k_B
from astropy.constants import mu0 as μ0
from numbers import Integral
from typing import List, Union

from plasmapy.formulary import Alfven_speed
from plasmapy.particles import Particle, particle_input, ParticleLike
from plasmapy.particles._factory import _physical_particle_factory
from plasmapy.particles.exceptions import InvalidIonError, ParticleError
from plasmapy.simulation.abstractions import AbstractNormalizations
from plasmapy.utils._units_helpers import _get_physical_type_dict

_length = u.get_physical_type(u.m)
_magnetic_field = u.get_physical_type(u.T)
_time = u.get_physical_type(u.s)
_number_density = u.get_physical_type(u.m**-3)
_velocity = u.get_physical_type(u.m / u.s)


class MHDNormalizations(AbstractNormalizations):
    r"""
    :term:`Normalizations` commonly used for the equations of
    magnetohydrodynamics.

    Parameters
    ----------
    *args : `tuple` of |Quantity| or |ParticleLike|

    Z : integer
        The charge number of the ion.

    mass_numb : integer
        The mass number of the isotope of the ion.

    Examples
    --------
    >>> import astropy.units as u
    >>> from plasmapy.simulation import MHDNormalizations
    >>> n0 = 1e19 * u.m**-3
    >>> L0 = 1 * u.km
    >>> B0 = 10 * u.G
    >>> normalizations = MHDNormalizations("p+", n0, L0, B0)
    >>> normalizations.ion
    Particle("p+")
    >>> normalizations.number_density
    <Quantity 1.e+19 1 / m3>
    >>> normalizations.magnetic_field
    <Quantity 0.001 T>
    >>> normalizations.length
    <Quantity 1000. m>
    """

    @particle_input
    def _store_ion(self, ion: Particle):
        self._ion = ion

    def _process_quantities(self, quantities: List[u.Quantity]):
        self._quantities = _get_physical_type_dict(quantities)
        self._get_length()
        self._get_magnetic_field()
        self._get_number_density()

    def _get_length(self):
        if _length in self._quantities:
            return

    def _get_magnetic_field(self):
        if _magnetic_field not in self._quantities:
            return

    def _get_number_density(self):
        if _number_density not in self._quantities:
            ...

    # Make sure that when there are duplicates of a physical type in
    # args, that there is a ValueError or something.  This is (or should
    # be) done in _get_physical_type_dict.

    def __init__(
        self,
        *args: Union[u.Quantity, ParticleLike],
        Z: Integral = None,
        mass_numb: Integral = None,
    ):

        wrong_type_error_message = (
            "MHDNormalizations can only accept Quantity positional "
            "arguments and one particle-like argument that represents "
            "an ion."
        )

        quantities = []
        particles = []
        for arg in args:
            if isinstance(arg, u.Quantity):
                quantities.append(arg.si)
            elif isinstance(arg, ParticleLike):
                particles.append(Particle(arg, Z=Z, mass_numb=mass_numb))
            else:
                raise TypeError(wrong_type_error_message)

        # Should we allow a ParticleList?
        # Should we allow the ion to not be defined?

        if len(particles) != 1:
            raise TypeError(wrong_type_error_message)

        self._store_ion(particles[0])
        self._process_quantities(quantities)

    @property
    def acceleration(self) -> u.m * u.s**-2:
        r"""
        The acceleration :term:`normalization`,

        .. math::

           a_⭑ ≡ \frac{L_⭑}{t_⭑^2}.

        Double check this fits in with MHD equations!

        """

    @property
    def current_density(self) -> u.A * u.m**-2:
        r"""
        The current density :term:`normalization`,

        .. math::

           J_⭑ ≡ \frac{B_⭑}{μ_0 L_⭑}.

        Returns
        -------
        |Quantity|
        """
        return self.magnetic_field / (μ0 * self.length)

    @property
    def diffusivity(self) -> u.m**2 / u.s:
        r"""
        The diffusivity :term:`normalization`,

        .. math:

           D_⭑

        Returns
        -------
        |Quantity|
        """

    @property
    def dynamic_viscosity(self) -> u.Pa * u.s:
        r"""
        The dynamic viscosity :term:`normalization`,

        .. math::

           _⭑ ≡

        Returns
        -------
        |Quantity|
        """

    @property
    def electric_field(self) -> u.V / u.m:
        """
        The electric field :term:`normalization`,

        .. math::

           E_⭑ ≡ V_⭑ B_⭑.

        Returns
        -------
        |Quantity|
        """
        return self.velocity * self.magnetic_field

    @property
    def frequency(self) -> u.s**-1:
        """
        The frequency :term:`normalization`,

        .. math::

           f_⭑ ≡ t_⭑^{-1}.

        Returns
        -------
        |Quantity|
        """

    @property
    def heat_flux(self) -> u.J * u.m**-2 * u.s**-1:
        r"""
        The heat flux :term:`normalization`,

        .. math::

           _⭑ ≡

        Returns
        -------
        |Quantity|
        """

    @property
    def ion(self) -> Particle:
        """
        The particle...

        Returns
        -------
        |Particle|
        """
        return self._ion

    @property
    def length(self) -> u.m:
        r"""
        The length :term:`normalization`, :math:`L_⭑`\ .

        Returns
        -------
        |Quantity|
        """
        return self._quantities[_length]

    @property
    def magnetic_field(self) -> u.T:
        r"""
        The magnetic field :term:`normalization`, :math:`B_⭑`\ .

        Returns
        -------
        |Quantity|
        """
        return self._quantities[_magnetic_field]

    @property
    def magnetic_flux(self) -> u.Wb:
        r"""
        The magnetic flux :term:`normalization`,

        .. math::

           _⭑ ≡

        Returns
        -------
        |Quantity|
        """

    @property
    def mass(self) -> u.kg:
        r"""
        The mass :term:`normalization`,

        .. math::

           m_⭑ ≡ ???

        Returns
        -------
        |Quantity|
        """
        raise NotImplementedError

    @property
    def mass_density(self) -> u.kg * u.m**-3:
        r"""
        The mass density :term:`normalization`, :math:`m_⭑=`\ ...

        Returns
        -------
        |Quantity|
        """
        raise NotImplementedError

    @property
    def number_density(self) -> u.m**-3:
        r"""
        The number density :term:`normalization`, :math:`n_⭑`\ .

        Returns
        -------
        |Quantity|
        """
        return self._quantities[_number_density]

    @property
    def pressure(self):
        r"""
        The pressure :term:`normalization`,

        .. math::

           p_⭑ ≡ \frac{B_⭑^2}{μ_0}.

        Returns
        -------
        |Quantity|
        """
        return self.magnetic_field**2 / μ0

    @property
    def resistivity(self) -> u.ohm * u.m:
        r"""
        The resistivity :term:`normalization`,

        .. math::

           η_⭑ ≡ μ_0 L_⭑ V_⭑.

        Returns
        -------
        |Quantity|
        """
        return μ0 * self.length

    @property
    def temperature(self) -> u.K:
        r"""
        The temperature :term:`normalization`,

        .. math::

           T_⭑ ≡ \frac{B_⭑^2}{k_B μ_0 n_⭑}

        Returns
        -------
        |Quantity|
        """
        return self.magnetic_field**2 / (k_B * μ0 * self.number_density)

    @property
    def thermal_conduction(self) -> u.W * u.K**-1 * u.m**-1:
        r"""
        The thermal conduction :term:`normalization`,

        .. math::

           _⭑ ≡

        Returns
        -------
        |Quantity|
        """

    @property
    def time(self) -> u.s:
        r"""
        The time :term:`normalization`,

        t_⭑ ≡ \frac{L_⭑}{V_⭑}.

        Returns
        -------
        |Quantity|
        """
        return self.length / self.velocity

    @property
    def velocity(self) -> u.m / u.s:
        r"""
        The velocity :term:`normalization`,

        .. math::

           V_⭑ ≡ \frac{B_⭑}{\sqrt{μ_0 ρ_⭑}}.

        Returns
        -------
        |Quantity|
        """
        return Alfven_speed(B=self.magnetic_field, density=self.mass_density)

    @property
    def volumetric_heating_rate(self) -> u.J * u.m**-3 * u.s**-1:
        r"""
        The volumetric heating rate :term:`normalization`,

        .. math::

           _⭑ ≡

        Returns
        -------
        |Quantity|

        Notes
        -----
        This :term:`normalization` is applicable to, for example, the number of
        collisions per unit volume per unit time.
        """

    @property
    def wavenumber(self) -> u.m**-1:
        r"""
        The wavenumber :term:`normalization`,

        .. math::

           k_⭑ ≡ \frac{1}{L_⭑}
        """
