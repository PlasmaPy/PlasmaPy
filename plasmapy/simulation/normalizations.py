import astropy.units as u

from plasmapy.particles import Particle
from plasmapy.simulation.abstractions import AbstractNormalizations


class MHDNormalizations(AbstractNormalizations):
    """
    Normalizations commonly used for the equations of magnetohydrodynamics.

    Parameters
    ----------
    *args
        Multiple `~astropy.unit.Quantity` objects with different
        physical types.

    Notes
    -----

    """

    _all_normalizations = [
        "current_density",
    ]

    def __init__(self, *args):
        pass

    @property
    def current_density(self) -> u.A * u.m ** -2:
        pass

    @property
    def ion(self) -> Particle:
        pass

    @property
    def magnetic_field(self) -> u.T:
        pass

    @property
    def mass_density(self) -> u.kg * u.m ** -3:
        pass

    @property
    def number_density(self) -> u.m ** -3:
        pass
