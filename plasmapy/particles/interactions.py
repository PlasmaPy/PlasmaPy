"""Classes to represent particle interactions and reactions."""

__all__ = ["AbstractParticleInteraction"]

from abc import ABC, abstractmethod

import astropy.units as u
import numpy as np

from plasmapy.particles.particle_collections import ParticleList, ParticleListLike
from plasmapy.utils.decorators import validate_quantities


class AbstractParticleInteraction(ABC):
    """..."""

    def __init__(self, reactants: ParticleList, products: ParticleList):
        self.reactants = reactants
        self.products = products
        self.validate_interaction()

    @property
    def reactants(self) -> ParticleList:
        """
        The reactants to the particle interaction.

        Returns
        -------
        `~plasmapy.particles.particle_collection.ParticleList`
        """
        return self._reactants

    @reactants.setter
    def reactants(self, reactants_: ParticleListLike) -> None:
        self._reactants = ParticleList(reactants_)

    @property
    def products(self) -> ParticleList:
        """
        The products from the particle interaction.

        Returns
        -------
        `~plasmapy.particles.particle_collection.ParticleList`
        """
        return self._products

    @products.setter
    def products(self, products_: ParticleListLike) -> None:
        self._products = ParticleList(products_)

    @abstractmethod
    def validate_interaction(self) -> None:
        """
        If the interaction violates any required conservation laws,
        raise a `~plasmapy.particles.exceptions.InvalidReactionError`.
        """


class NuclearReaction(AbstractParticleInteraction):
    """..."""

    def validate_interaction(self) -> None:
        """..."""

    @validate_quantities
    def energy(self) -> u.Quantity[u.J]:
        """Return the mass energy of the nuclear reaction."""
        return np.sum(self.reactants.mass_energy) - np.sum(self.products.mass_energy)
