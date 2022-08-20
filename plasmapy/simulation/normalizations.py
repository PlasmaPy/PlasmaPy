import astropy.units as u
import numpy as np

from astropy.constants import k_B
from astropy.constants import mu0 as μ0
from typing import List

from plasmapy.formulary import Alfven_speed
from plasmapy.particles import Particle, particle_input, ParticleLike
from plasmapy.particles._factory import _physical_particle_factory
from plasmapy.particles.exceptions import InvalidIonError, ParticleError
from plasmapy.simulation.abstractions import AbstractNormalizations
from plasmapy.utils._units_helpers import _get_physical_type_dict

_LENGTH = u.get_physical_type(u.m)
_MAGNETIC_FIELD = u.get_physical_type(u.T)
_TIME = u.get_physical_type(u.s)
_NUMBER_DENSITY = u.get_physical_type(u.m**-3)
_VELOCITY = u.get_physical_type(u.m / u.s)


class MHDNormalizations(AbstractNormalizations):
    r"""
    Normalizations commonly used for the equations of magnetohydrodynamics.

    Parameters
    ----------
    *args : `tuple` of |Quantity| or |ParticleLike|

    Z : integer
        The charge number of the ion.

    mass_numb : integer
        The mass number of the isotope of the ion.
    """

    def _store_ion(self, ion: Particle):
        pass

    def _store_quantites(self, quantities: List[u.Quantity]):
        pass

    def __init__(self, *args, Z=None, mass_numb=None):

        wrong_type_error_message = (
            "MHDNormalizations can only accept Quantity positinal "
            "arguments and one particle-like argument that represents "
            "an ion."
        )

        quantities = []
        particles = []
        for arg in args:
            if isinstance(arg, u.Quantity):
                quantities.append(arg)
            elif isinstance(arg, ParticleLike):
                particles.append(Particle(arg, Z=Z, mass_numb=mass_numb))
            else:
                raise TypeError(wrong_type_error_message)

        # Should we allow the ion to not be defined?
        if len(particles) != 1:
            raise TypeError(wrong_type_error_message)

        self._store_ion(particles[0])
        self._store_quantities(quantities)

    @property
    def current_density(self) -> u.A * u.m**-2:
        r"""
        The current density normalization,

        .. math::

           J_⭑ ≡ \frac{B_⭑}{μ_0 L_⭑}.

        Returns
        -------
        |Quantity|
        """
        return self.magnetic_field / (μ0 * self.length)

    @property
    def ion(self) -> Particle:
        """
        The particle...

        Returns
        -------
        |Particle|
        """
        pass

    @property
    def length(self) -> u.m:
        r"""
        The length normalization, :math:`L_⭑`\ .

        Returns
        -------
        |Quantity|
        """
        return self._length

    @property
    def magnetic_field(self) -> u.T:
        r"""
        The magnetic field normalization, :math:`B_⭑`\ .

        Returns
        -------
        |Quantity|
        """
        return self._magnetic_field

    @property
    def mass(self) -> u.kg:
        r"""
        The mass normalization,

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
        The mass density normalization, :math:`m_⭑=`\ ...

        Returns
        -------
        |Quantity|
        """
        raise NotImplementedError

    @property
    def number_density(self) -> u.m**-3:
        r"""
        The number density normalization, :math:`n_⭑`\ .

        Returns
        -------
        |Quantity|
        """
        return self._number_density

    @property
    def pressure(self):
        r"""
        The pressure normalization,

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
        The resistivity normalization,

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
        The temperature normalization,

        .. math::

           T_⭑ ≡ \frac{B_⭑^2}{k_B μ_0 n_⭑}

        Returns
        -------
        |Quantity|
        """
        return self.magnetic_field**2 / (k_B * μ0 * self.number_density)

    @property
    def time(self) -> u.s:
        r"""
        The time normalization,

        t_⭑ ≡ \frac{L_⭑}{V_⭑}.

        Returns
        -------
        |Quantity|
        """
        return self.length / self.velocity

    @property
    def velocity(self) -> u.m / u.s:
        r"""
        The velocity normalization,

        .. math::

           V_⭑ ≡ \frac{B_⭑}{\sqrt{μ_0 ρ_⭑}}.

        Returns
        -------
        |Quantity|
        """
        return Alfven_speed(B=self.magnetic_field, density=self.mass_density)
