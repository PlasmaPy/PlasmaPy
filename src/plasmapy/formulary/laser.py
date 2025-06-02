"""
Functions for calculating quantities associated with laser pulses.

.. attention::

   |expect-api-changes|
"""

__all__ = [
    "electric_field_amplitude",
]
__aliases__ = ["E0_"]

import astropy.units as u
import numpy as np
from astropy.constants.si import c, eps0

from plasmapy.utils.decorators import validate_quantities

__all__ += __aliases__

@validate_quantities(
    intensity={"can_be_negative": False},
)
def electric_field_amplitude(
    intensity: u.Quantity[u.watt / u.m**2],
) -> u.Quantity[u.V / u.m]:
    r"""
    Calculate the electric field amplitude :math:`E_0` from the intensity :math:`I`
    of a laser.

    The electric field amplitude of an electromagnetic plane wave in vacuum
    is calculated using:

    .. math::
        E_0=\sqrt{\frac{2I}{c ε_0}},

    where :math:`c` is the speed of light and
    :math:`ε_0` is the permittivity of free space.

    **Aliases:** `E0_`

    Parameters
    ----------
    intensity : `~astropy.units.Quantity`
        Intensity of the laser pulse (convertible to W / m\:sup:`2`).

    Returns
    -------
    E : `~astropy.units.Quantity`
        Maximum electric field amplitude for the intensity provided.

    Notes
    -----
    For details, see :cite:t:`ling:2016`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> electric_field_ampluitude(1e-3 * u.watt / u.m**2)  # Electric Field Strength
    <Quantity 0.8680211 V / m>
    """

    E = np.sqrt((2 * intensity) / (c * eps0))
    return E.to(u.V / u.m)


E0_ = electric_field_amplitude
"""Alias to `~plasmapy.formulary.laser.electric_field_amplitude`."""
