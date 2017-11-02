"""Functions to deal with distribution : generate, fit, calculate"""

from astropy import units

from astropy.units import (UnitConversionError, UnitsError, quantity_input,
                           Quantity)

from ..constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, e)
from ..atomic import (ion_mass, charge_state)
from ..atomic.atomic import _is_electron as is_electron
import numpy as np

from ..utils import _check_quantity, check_relativistic, check_quantity


@units.quantity_input
def Maxwellian_1D(v: units.m/units.s,
                  T: units.K, particle="e",
                  V_drift=0*units.m/units.s):
    r"""Returns the probability at the velocity `v` in m/s
    to find a particle `particle` in a plasma of temperature `T`
    following the Maxwellian distribution function.

    Parameters
    ----------
    v: Quantity
        The velocity in units convertible to m/s

    T: Quantity
        The temperature in Kelvin

    particle: string, optional
        Representation of the particle species(e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for $He_4^{+1}$ : singly ionized helium-4),
        which defaults to electrons.

    Returns
    -------
    f : Quantity
        probability in Velocity^-1, normized so that: :math:`\int_{-\infty}^{+\infty} f(v) dv = 1`

    Raises
    ------
    TypeError
        The parameter arguments are not Quantities and
        cannot be converted into Quantities.

    UnitConversionError
        If the parameters is not in appropriate units.

    ValueError
        If the temperature is negative, or the particle mass or charge state
        cannot be found.

    Notes
    -----
    In one dimension, the Maxwellian distribution function for a particle of
    mass m, velocity v, a drift velocity V and with temperature T is:

    .. math::
        f = \sqrt{\frac{m}{2 \pi k_B T}} e^{-\frac{m}{2 k_B T} (v-V)^2}

    Examples
    --------
    >>> from plasmapy.physics import Maxwellian_1D
    >>> from astropy import units as u
    >>> v=1*u.m/u.s
    >>> Maxwellian_1D(v=v, T= 30000*u.K, particle='e',V_drift=0*u.m/u.s)
    <Quantity 5.916329687405701e-07 s / m>
    """

    # pass to Kelvin

    if T.unit.is_equivalent(units.K):
        T = T.to(units.K, equivalencies=units.temperature_energy())

    # Get mass

    if is_electron(particle):
        m_s = m_e

    else:
        try:
            m_s = ion_mass(particle)
        except Exception:
            raise ValueError("Unable to find ion mass in Maxwellian_1D")

    T = k_B*T

    f = np.sqrt((m_s/(2*pi*T)))*np.exp(-m_s/(2*T)*(v-V_drift)**2)

    return f.to(units.s/units.m)
