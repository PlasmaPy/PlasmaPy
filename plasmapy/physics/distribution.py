"""Functions to deal with distribution : generate, fit, calculate"""

from astropy import units as u

from astropy.units import (UnitConversionError, UnitsError, quantity_input,
                           Quantity)

from ..constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, e)
from ..atomic import (ion_mass, charge_state)
from ..atomic.atomic import _is_electron as is_electron
from .parameters import thermal_speed
import numpy as np

from ..utils import _check_quantity, check_relativistic, check_quantity


@u.quantity_input
def Maxwellian_1D(v: u.m/u.s,
                  T, 
                  particle="e",
                  V_drift=0*u.m/u.s):
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
        
    V_drift: Quantity, optional
        The drift velocity in units convertible to m/s

    Returns
    -------
    f : Quantity
        probability in Velocity^-1, normalized so that: :math:`\int_{-\infty}^{+\infty} f(v) dv = 1`

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
    	f = (\pi * v_Th^2)^{-1/2} e^{-(v - v_{drift})^2 / v_Th^2}
    	where v_Th = \sqrt(2 k_B T / m) is the thermal speed

    Examples
    --------
    >>> from plasmapy.physics import Maxwellian_1D
    >>> from astropy import units as u
    >>> v=1*u.m/u.s
    >>> Maxwellian_1D(v=v, T= 30000*u.K, particle='e',V_drift=0*u.m/u.s)
    <Quantity 5.916329687405703e-07 s / m>
    """
    # convert temperature to Kelvins
    T = T.to(u.K, equivalencies=u.temperature_energy())
    # get thermal velocity and thermal velocity squared
    vTh = thermal_speed(T, particle=particle, method="most_probable")
    vThSq = vTh ** 2
    # Get square of relative particle velocity
    vSq = (v - V_drift) ** 2
    # calculating distribution function
    try:
        coeff = (vThSq * np.pi) ** (-1 / 2)
        expTerm = np.exp(-vSq / vThSq)
        distFunc = coeff * expTerm
    except Exception:
        raise ValueError("Unable to get distribution function.")
    return distFunc.to(u.s/u.m)

#@u.quantity_input
#def Maxwellian_velocity_3D(vx: u.m/u.s,
#                           vy: u.m/u.s,
#                           vz: u.m/u.s,
#                           T, 
#                           particle="e",
#                           Vx_drift=0*u.m/u.s,
#                           Vy_drift=0*u.m/u.s,
#                           Vz_drift=0*u.m/u.s):
#    r"""Return the probability of finding a particle with velocity components
#    `v_x`, `v_y`, and `v_z`in m/s in an equilibrium plasma of temperature 
#    `T` which follows the 3D Maxwellian distribution function. This 
#    function assumes Cartesian coordinates.
#
#    Parameters
#    ----------
#    vx: Quantity
#        The velocity in x-direction units convertible to m/s
#        
#    vy: Quantity
#        The velocity in y-direction units convertible to m/s
#        
#    vz: Quantity
#        The velocity in z-direction units convertible to m/s
#
#    T: Quantity
#        The temperature, preferably in Kelvin
#
#    particle: string, optional
#        Representation of the particle species(e.g., 'p' for protons, 'D+'
#        for deuterium, or 'He-4 +1' for $He_4^{+1}$ : singly ionized helium-4),
#        which defaults to electrons.
#        
#    Vx_drift: Quantity
#        The drift velocity in x-direction units convertible to m/s
#        
#    Vy_drift: Quantity
#        The drift velocity in y-direction units convertible to m/s
#        
#    Vz_drift: Quantity
#        The drift velocity in z-direction units convertible to m/s
#
#    Returns
#    -------
#    f : Quantity
#        probability in Velocity^-1, normalized so that:
#        $\iiint_{0}^{\infty} f(\vec{v}) d\vec{v} = 1}
#
#    Raises
#    ------
#    TypeError
#        The parameter arguments are not Quantities and
#        cannot be converted into Quantities.
#
#    UnitConversionError
#        If the parameters is not in appropriate units.
#
#    ValueError
#        If the temperature is negative, or the particle mass or charge state
#        cannot be found.
#
#    Notes
#    -----
#    In one dimension, the Maxwellian speed distribution function describing
#    the distribution of particles with speed v in a plasma with temperature T
#    is given by:
#
#    .. math::
#    f = (\pi * v_Th^2)^{-3/2} \exp(-(\vec{v} - \vec{V_{drift}})^2 / v_Th^2)
#    where v_Th = \sqrt(2 k_B T / m) is the thermal speed
#
#    See also
#    --------
#    Maxwellian_1D
#
#    Example
#    -------
#    >>> from plasmapy.physics import Maxwellian_velocity_3D
#    >>> from astropy import units as u
#    >>> v=1*u.m/u.s
#    >>> Maxwellian_velocity_3D(vx=v,
#    ... vy=v,
#    ... vz=v,
#    ... T=30000*u.K,
#    ... particle='e',
#    ... Vx_drift=0*u.m/u.s,
#    ... Vy_drift=0*u.m/u.s,
#    ... Vz_drift=0*u.m/u.s)
#    <Quantity 3.985430307328086e-20 s3 / m3>
#    
#    
#    """
#    # convert temperature to Kelvins
#    T = T.to(u.K, equivalencies=u.temperature_energy())
#    # get thermal velocity and thermal velocity squared
#    vTh = thermal_speed(T, particle=particle, method="most_probable")
#    # accounting for thermal velocity in 3D
#    vTh3D = np.sqrt(3) * vTh
#    vThSq = vTh3D ** 2
#    # Get square of relative particle velocity
#    vSq = ((vx-Vx_drift) ** 2 + (vy-Vy_drift) ** 2 + (vz-Vz_drift) ** 2)
#    # calculating distribution function
#    try:
#        coeff = (vThSq * np.pi) ** (-3 / 2)
#        expTerm = np.exp(-vSq / vThSq)
#        distFunc = coeff * expTerm
#    except Exception:
#        raise ValueError("Unable to get distribution function.")
#    return distFunc.to((u.s/u.m) ** 3)

def Maxwellian_velocity_3D(vx,
                           vy,
                           vz,
                           T, 
                           particle="e",
                           Vx_drift=0,
                           Vy_drift=0,
                           Vz_drift=0):
    r"""Return the probability of finding a particle with velocity components
    `v_x`, `v_y`, and `v_z`in m/s in an equilibrium plasma of temperature 
    `T` which follows the 3D Maxwellian distribution function. This 
    function assumes Cartesian coordinates.

    Parameters
    ----------
    vx: Quantity
        The velocity in x-direction units convertible to m/s
        
    vy: Quantity
        The velocity in y-direction units convertible to m/s
        
    vz: Quantity
        The velocity in z-direction units convertible to m/s

    T: Quantity
        The temperature, preferably in Kelvin

    particle: string, optional
        Representation of the particle species(e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for $He_4^{+1}$ : singly ionized helium-4),
        which defaults to electrons.
        
    Vx_drift: Quantity
        The drift velocity in x-direction units convertible to m/s
        
    Vy_drift: Quantity
        The drift velocity in y-direction units convertible to m/s
        
    Vz_drift: Quantity
        The drift velocity in z-direction units convertible to m/s

    Returns
    -------
    f : Quantity
        probability in Velocity^-1, normalized so that:
        $\iiint_{0}^{\infty} f(\vec{v}) d\vec{v} = 1}

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
    In one dimension, the Maxwellian speed distribution function describing
    the distribution of particles with speed v in a plasma with temperature T
    is given by:

    .. math::
    f = (\pi * v_Th^2)^{-3/2} \exp(-(\vec{v} - \vec{V_{drift}})^2 / v_Th^2)
    where v_Th = \sqrt(2 k_B T / m) is the thermal speed

    See also
    --------
    Maxwellian_1D

    Example
    -------
    >>> from plasmapy.physics import Maxwellian_velocity_3D
    >>> from astropy import units as u
    >>> v=1*u.m/u.s
    >>> Maxwellian_velocity_3D(vx=v,
    ... vy=v,
    ... vz=v,
    ... T=30000*u.K,
    ... particle='e',
    ... Vx_drift=0*u.m/u.s,
    ... Vy_drift=0*u.m/u.s,
    ... Vz_drift=0*u.m/u.s)
    <Quantity 3.985430307328086e-20 s3 / m3>
    
    
    """
#    # convert temperature to Kelvins
#    T = T.to(u.K, equivalencies=u.temperature_energy())
#    # get thermal velocity and thermal velocity squared
#    vTh = (thermal_speed(T, particle=particle, method="most_probable")).value
    vTh = np.sqrt(1.6e-19 * T / (2 * 9.11e-31))
    # accounting for thermal velocity in 3D
    vThSq = 3 * vTh ** 2
#    vTh3D = np.sqrt(3) * vTh
#    vThSq = vTh3D ** 2
    # Get square of relative particle velocity
    vSq = ((vx-Vx_drift) ** 2 + (vy-Vy_drift) ** 2 + (vz-Vz_drift) ** 2)
    # calculating distribution function
    coeff = (vThSq * np.pi) ** (-3 / 2)
    expTerm = np.exp(-vSq / vThSq)
    distFunc = coeff * expTerm
    return distFunc

@u.quantity_input
def Maxwellian_speed_1D(v: u.m/u.s,
                        T, 
                        particle="e",
                        V_drift=0*u.m/u.s):
    r"""Return the probability of finding a particle with speed `v` in m/s
     in an equilibrium plasma of temperature `T` which follows the 
     Maxwellian distribution function.

    Parameters
    ----------
    v: Quantity
        The speed in units convertible to m/s

    T: Quantity
        The temperature, preferably in Kelvin

    particle: string, optional
        Representation of the particle species(e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for $He_4^{+1}$ : singly ionized helium-4),
        which defaults to electrons.
    
    V_drift: Quantity
        The drift speed in units convertible to m/s

    Returns
    -------
    f : Quantity
        probability in speed^-1, normalized so that:
        $\int_{0}^{\infty} f(v) dv = 1}

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
    In one dimension, the Maxwellian speed distribution function describing
    the distribution of particles with speed v in a plasma with temperature T
    is given by:

    .. math::
    f(v) = 4 \pi v^2 (\pi * v_Th^2)^{-3/2} \exp(-(v - V_{drift})^2 / v_Th^2)
    where v_Th = \sqrt(2 k_B T / m) is the thermal speed


    Example
    -------
    >>> from plasmapy.physics import Maxwellian_speed_1D
    >>> from astropy import units as u
    >>> v=1*u.m/u.s
    >>> Maxwellian_speed_1D(v=v, T= 30000*u.K, particle='e',V_drift=0*u.m/u.s)
    <Quantity 2.602357544747327e-18 s / m>
    
    """
    # convert temperature to Kelvins
    T = T.to(u.K, equivalencies=u.temperature_energy())
    # get thermal velocity and thermal velocity squared
    vTh = thermal_speed(T, particle=particle, method="most_probable")
    vThSq = vTh ** 2
    # get square of relative particle velocity
    vSq = (v - V_drift) ** 2
    # calculating distribution function
    try:
        coeff1 = (np.pi * vThSq) ** (-3 / 2)
        coeff2 = 4 * np.pi * vSq
        expTerm = np.exp(-vSq / vThSq)
        distFunc = coeff1 * coeff2 * expTerm
    except Exception:
        raise ValueError("Unable to get distribution function.")
    return distFunc.to(u.s/u.m)

@u.quantity_input
def Maxwellian_speed_3D(vx: u.m/u.s,
                        vy: u.m/u.s,
                        vz: u.m/u.s,
                        T, 
                        particle="e",
                        Vx_drift=0*u.m/u.s,
                        Vy_drift=0*u.m/u.s,
                        Vz_drift=0*u.m/u.s):
    r"""Return the probability of finding a particle with speed components
    `v_x`, `v_y`, and `v_z`in m/s in an equilibrium plasma of temperature 
    `T` which follows the 3D Maxwellian distribution function. This 
    function assumes Cartesian coordinates.

    Parameters
    ----------
    vx: Quantity
        The speed in x-direction units convertible to m/s
        
    vy: Quantity
        The speed in y-direction units convertible to m/s
        
    vz: Quantity
        The speed in z-direction units convertible to m/s

    T: Quantity
        The temperature, preferably in Kelvin

    particle: string, optional
        Representation of the particle species(e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for $He_4^{+1}$ : singly ionized helium-4),
        which defaults to electrons.
        
    Vx_drift: Quantity
        The drift speed in x-direction units convertible to m/s
        
    Vy_drift: Quantity
        The drift speed in y-direction units convertible to m/s
        
    Vz_drift: Quantity
        The drift speed in z-direction units convertible to m/s
        
    Returns
    -------
    f : Quantity
        probability in speed^-1, normalized so that:
        $\iiint_{0}^{\infty} f(\vec{v}) d\vec{v} = 1}

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
    In one dimension, the Maxwellian speed distribution function describing
    the distribution of particles with speed v in a plasma with temperature T
    is given by:

    .. math::
    f = 4 \pi \vec{v}^2 (\pi * v_Th^2)^{-3/2} \exp(-(\vec{v} - \vec{V_{drift}})^2 / v_Th^2)
    where v_Th = \sqrt(2 k_B T / m) is the thermal speed

    See also
    --------
    Maxwellian_speed_1D

    Example
    -------
    >>> from plasmapy.physics import Maxwellian_speed_3D
    >>> from astropy import units as u
    >>> v=1*u.m/u.s
    >>> Maxwellian_speed_3D(vx=v,
    ... vy=v,
    ... vz=v,
    ... T=30000*u.K,
    ... particle='e',
    ... Vx_drift=0*u.m/u.s,
    ... Vy_drift=0*u.m/u.s,
    ... Vz_drift=0*u.m/u.s)
    <Quantity 7.80707263422481e-18 s / m>
    
    """
    # convert temperature to Kelvins
    T = T.to(u.K, equivalencies=u.temperature_energy())
    # Get particle speed and drift speed
    v = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    V_drift = np.sqrt(Vx_drift ** 2 + Vy_drift ** 2 + Vz_drift ** 2)
    # calculating distribution function
    try:
        distFunc = Maxwellian_speed_1D(v=v,
                                       T=T, 
                                       particle=particle,
                                       V_drift=V_drift)
    except Exception:
        raise ValueError("Unable to get distribution function.")
    return distFunc.to(u.s/u.m)
