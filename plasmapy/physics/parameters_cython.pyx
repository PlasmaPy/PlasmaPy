"""
Optimized versions of functions in test_parameters.py.
The optimization is mainly due to static typing.
Be careful when using these function as they do have fixed precision, due
to using C types.
"""

# python modules
import numpy as np
from astropy import units

# plasmapy modules
from plasmapy.constants import k_B
import plasmapy.atomic as atomic
# from plasmapy.atomic import particle_mass, charge_state
# For future: change these into decorators.  _check_quantity does a
# bit more than @quantity_input as it can allow
import plasmapy.utils as utils
from plasmapy.utils.checks import _check_quantity
from plasmapy.physics.exceptions import PhysicsError  # , PhysicsWarning


# @utils.check_relativistic
# @utils.check_quantity({
#    'T': {'units': units.K, 'can_be_negative': False}
# })
def thermal_speed(double T, particle="e", method="most_probable"):
    r"""
    Returns the most probable speed for a particle within a Maxwellian
    distribution.

    Parameters
    ----------
    T : ~astropy.units.Quantity
        The particle temperature in either kelvin or energy per particle

    particle : string, optional
        Representation of the particle species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4),
        which defaults to electrons.  If no charge state information is
        provided, then the particles are assumed to be singly charged.

    method : string, optional
        Method to be used for calculating the thermal speed. Options are
        'most_probable' (default), 'rms', and 'mean_magnitude'.

    Returns
    -------
    V : ~astropy.units.Quantity
        particle thermal speed

    Raises
    ------
    TypeError
        The particle temperature is not a Quantity

    UnitConversionError
        If the particle temperature is not in units of temperature or
        energy per particle

    ValueError
        The particle temperature is invalid or particle cannot be used to
        identify an isotope or particle

    UserWarning
        If the particle thermal speed exceeds 10% of the speed of light, or
        if units are not provided and SI units are assumed.

    Notes
    -----
    The particle thermal speed is given by:

    .. math::
        V_{th,i} = \sqrt{\frac{2 k_B T_i}{m_i}}

    This function yields the most probable speed within a distribution
    function.  However, the definition of thermal velocity varies by
    the square root of two depending on whether or not this velocity
    absorbs that factor in the expression for a Maxwellian
    distribution.  In particular, the expression given in the NRL
    Plasma Formulary [1] is a square root of two smaller than the
    result from this function.

    Examples
    --------
    >>> from astropy import units as u
    >>> thermal_speed(5*u.eV, 'p')
    <Quantity 30949.690182856546 m / s>
    >>> thermal_speed(1e6*u.K, particle='p')
    <Quantity 128486.55193256242 m / s>
    >>> thermal_speed(5*u.eV)
    <Quantity 1326205.1212395933 m / s>
    >>> thermal_speed(1e6*u.K)
    <Quantity 5505693.988425379 m / s>
    >>> thermal_speed(1e6*u.K, method="rms")
    <Quantity 6743070.475775486 m / s>
    >>> thermal_speed(1e6*u.K, method="mean_magnitude")
    <Quantity 6212510.3969422 m / s>

    """
    # static typing of variables
    cdef double m, V
#    T = T.to(units.K, equivalencies=units.temperature_energy())

    try:
        m = atomic.particle_mass(particle).si.value
    except Exception:
        raise ValueError("Unable to find {particle} mass in thermal_speed")

    # different methods, as per https://en.wikipedia.org/wiki/Thermal_velocity
    if method == "most_probable":
        V = np.sqrt(2 * k_B.si.value * T / m)
    elif method == "rms":
        V = np.sqrt(3 * k_B.si.value * T / m)
    elif method == "mean_magnitude":
        V = np.sqrt(8 * k_B.si.value * T / (m * np.pi))
    else:
        raise ValueError("Method {method} not supported in thermal_speed")

#    return V.to(units.m / units.s)
    return V
