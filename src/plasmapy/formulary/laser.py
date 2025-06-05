"""
Functions for calculating quantities associated with laser pulses.

.. attention::

   |expect-api-changes|
"""

__all__ = [
    "electric_field_amplitude",
    "intensity",
]
__aliases__ = ["E0_", "I_"]

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
        Intensity of the laser pulse (convertible to W / m\ :sup:`2`).

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
    >>> electric_field_amplitude(1e-3 * u.watt / u.m**2)  # Electric Field Amplitude
    <Quantity 0.8680211 V / m>
    """

    E = np.sqrt((2 * intensity) / (c * eps0))
    return E.to(u.V / u.m)


E0_ = electric_field_amplitude
"""Alias to `~plasmapy.formulary.laser.electric_field_amplitude`."""


@validate_quantities(
    electric_field_amplitude={"can_be_negative": False},
)
def intensity(
    electric_field_amplitude: u.Quantity[u.V / u.m],
) -> u.Quantity[u.watt / u.m**2]:
    r"""
    Calculate the intensity :math:`I` of a laser from the
    electric field amplitude :math:`E_0`.

    The intensity of an electromagnetic plane wave in vacuum
    is calculated using:

    .. math::
        I=\frac{1}{2} c ε_0 E_0^2,

    where :math:`c` is the speed of light and
    :math:`ε_0` is the permittivity of free space.

    **Aliases:** `I_`

    Parameters
    ----------
    electric_field_amplitude: `~astropy.units.Quantity`
        Electric field amplitude of an electromagnetic plane wave
        (convertible to V / m).

    Returns
    -------
    Int : `~astropy.units.Quantity`
        Intensity for the electric field amplitude provided.

    Notes
    -----
    For details, see :cite:t:`ling:2016`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> intensity(0.8680211 * u.V / u.m)  # Intensity
    <Quantity 0.001 W / m2>
    """

    Int = (1 / 2) * c * eps0 * electric_field_amplitude**2
    return Int.to(u.Watt / u.m**2)


I_ = intensity
"""Alias to `~plasmapy.formulary.laser.intensity`."""
