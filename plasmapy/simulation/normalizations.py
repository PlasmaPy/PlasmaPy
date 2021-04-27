import astropy.units as u

from plasmapy.particles import Particle
from plasmapy.simulation.abstractions import AbstractNormalizations


class MHDNormalizations(AbstractNormalizations):
    """
    Normalizations commonly used for the equations of magnetohydrodynamics.

    Parameters
    ----------
    B : `~astropy.units.Quantity`, keyword-only, optional
        The magnetic field normalization in units convertible to tesla.

    n : `~astropy.units.Quantity`, keyword-only, optional
        The number density normalization in units convertible to
        m:sup:`-3`\ .

    L : `~astropy.units.Quantity`, keyword-only, optional
        The length normalization in units convertible to m.

    V : `~astropy.units.Quantity`, keyword-only, optional
        The velocity normalization in units convertible to m/s.

    time : `~astropy.units.Quantity` or callable, keyword-only, optional
        The time normalization, or a function that returns a

    Notes
    -----

    """

    _all_normalizations = [
        "current_density",
    ]

    def __init__(self, *, B=None, n=None, L=None, V=None, time=None):
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
