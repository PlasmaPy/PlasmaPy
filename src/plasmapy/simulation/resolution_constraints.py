"""
Grid and time resolution constraints for numerical simulations.

.. attention::

   |expect-api-changes|
"""

__all__ = [
    "CFL_limit_electromagnetic_yee",
]


import astropy.units as u
import numpy as np
from astropy.constants.si import c

from plasmapy.utils.decorators import validate_quantities


@validate_quantities(
    dx={
        "can_be_negative": False,
        "can_be_inf": False,
        "can_be_nan": False,
        "can_be_zero": False,
    },
)
def CFL_limit_electromagnetic_yee(dx: u.Quantity[u.m]) -> u.Quantity[u.s]:
    r"""
    Calculates the limiting time-step for a finite difference time-domain
    electromagnetic Yee solver which uses a Cartesian grid.

    This limit is defined by the Courant-Friederichs-Lewy (CFL) Condition:

    .. math::
        \Delta t = \frac{1}{c\sqrt{{\sum_{i=1}^{n} \Delta x_i^{-2}}}}

    where :math:`\Delta x_i` corresponds to the grid resolution along the
    :math:`i`\ -th dimension and :math:`n` is the number of dimensions.
    For example, in 3D this corresponds to:

    .. math::
        \Delta t =
        \frac{1}{c\sqrt{\frac{1}{\Delta x^2} +\frac{1}{\Delta y^2} +\frac{1}{\Delta z^2}}}

    Parameters
    ----------
    dx: `~astropy.units.Quantity`
        Array with grid resolution along the different dimensions.

    Returns
    -------
    dt: `~astropy.units.Quantity`
        Computed CFL limiting time-step.

    Notes
    -----
    For details, see :cite:p:`courant:1928,courant:1967`

    Examples
    --------
    >>> import astropy.units as u
    >>> import numpy as np
    >>> CFL_limit_electromagnetic_yee(10 * u.nm)
    <Quantity 3.33564...e-17 s>
    >>> CFL_limit_electromagnetic_yee(np.array([5, 10, 15]) * u.nm)
    <Quantity 1.42956...e-17 s>
    """

    return 1 / (c * np.sqrt(np.sum(1 / (dx**2))))
