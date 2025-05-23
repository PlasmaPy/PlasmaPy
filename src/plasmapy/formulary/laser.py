"""
Functions for calculating quantities associated with laser pulses.

.. attention::

   |expect-api-changes|
"""

__all__ = [
    "E0",
]

import astropy.units as u
import numpy as np
from astropy.constants.si import c, eps0

from plasmapy.utils.decorators import validate_quantities


@validate_quantities(
    Intensity={"can_be_negative": False},
)
def E0(
    Intensity: u.Quantity[u.watt / u.m**2],
) -> u.Quantity[u.V / u.m]:
    r"""
    Calculate Electric Field :math:`E_0` from Intensity :math:`I`.

    The E-field is calculated using:

    .. math::
        E_0=\sqrt{\frac{2I}{c ε_0}},

    where :math:`c` is the speed of light and
    :math:`ε_0` is the permittivity of free space.

    Parameters
    ----------
    Intensity : `~astropy.units.Quantity`
        Intensity of the laser pulse (convertible to u.Watts/u.m**2).

    Returns
    -------
    E : `~astropy.units.Quantity`
        Maximum Electric Field strength for the intensity provided.

    Notes
    -----
    For details, see :cite:t:`ling:2016`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> import numpy as np
    >>> E0(1e-3 * u.watt / u.m**2)  # Electric Field Strength
    <Quantity 0.8680211 V / m>
    """

    E = np.sqrt((2 * Intensity) / (c * eps0))
    return E.to(u.V / u.m)
