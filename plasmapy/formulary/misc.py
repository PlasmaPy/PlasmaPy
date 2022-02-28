"""Functions to for miscellaneous plasma parameter calculations."""

__all__ = []
__aliases__ = []

from plasmapy import particles
from plasmapy.particles import Particle


def _grab_charge(ion: Particle, z_mean=None):
    """
    Merge two possible inputs for particle charge.

    Parameters
    ----------
    ion : `~plasmapy.particles.particle_class.Particle`
        a string representing a charged particle, or a Particle object.

    z_mean : `float`
        An optional float describing the average ionization of a particle
        species.

    Returns
    -------
    float
        if ``z_mean`` was passed, ``z_mean``, otherwise, the charge number
        of ``ion``.

    """
    if z_mean is None:
        # warnings.warn("No z_mean given, defaulting to atomic charge",
        #               PhysicsWarning)
        Z = particles.charge_number(ion)
    else:
        # using average ionization provided by user
        Z = z_mean
    return Z

