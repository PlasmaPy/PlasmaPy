"""Defines the basic magnetic pickup probe (Bdot probe) analysis."""

__all__ = ["compute_bfield"]

from warnings import warn
from scipy.integrate import cumulative_trapezoid
from typing import Tuple
import numpy as np


import astropy.units as u
from plasmapy.utils.decorators import validate_quantities
from typing import Tuple

@validate_quantities
def compute_bfield(bdot_voltage: u.Quantity[u.volt], time_array: u.Quantity[u.s], loop_area: u.Quantity[u.m**2], num_loop: int = 1, gain: float = 1.0) -> u.Quantity[u.T]:
    r"""
    Take the voltage output of a magnetic pickup probe (Bdot probe) and time base and 
    return the associated array of the magnetic field as a function of time.

    Parameters
    ----------
    bdot_voltage: `astropy.units.Quantity`
        A data array containing voltage fluctuations from a bdot probe
        in units convertible to volts. The array values are proportional
        to time changing magnetic field via Faradays Law.

    time_array: `astropy.units.Quantity`
        The time series to the data collection in units of seconds.

    loop_area: 'astropy.units.Quantity'
        The area through which the changing flux is measured in units
        convertible to square meters.

    num_loop: int
        The number of loops of the Bdot probe (assumed to be all the same area).
        
    gain: float
        Any dimensionless gains that need to be applied.
        
    Returns
    -------
    magnetic_field_array: `astropy.units.Quantity`
        The array containing the magnetic field in units convertible
        to tesla as a function of time in seconds.
        
    Notes
    -----
    The integral form of Faraday's Law states that the time change in flux through an area generates a voltage around a loop,
    .. math::
       V = -N\frac{d\Phi}{dt}
    where :math:`\Phi` is the magnetic flux
    .. math::
        \Phi = \int{BdA}
    and :math:`N` is the number of loops.


    A magnetic pickup probe works by providing a fixed area loop or loops using a length of conducting wire through 
    which magnetic fields penetrate. As these fields change in time, a voltage is generated that can be measured
    by some kind of data acquisition device (such as an oscilloscope). We can also assume that the loop is small 
    enough that the magnetic field through the loop is constant everywhere in the loop (and only changes in time).
    Given this, we can write
    .. math::
        V = -A\frac{dB}{dt}
    where the voltage becomes proportional only to the change in magnetic field and the area of the loop. Magnetic
    field as a function of time can be recovered by time integrating :math:`V(t)/NA`.
    
    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> time_array = np.linspace(0, 4*np.pi, 100) * u.s
    >>> bdot_voltage = np.sin(time_array / u.s) * u.V  # create a voltage oscillation
    >>> loop_area = 1 * u.m**2
    >>> num_loop = 1
    >>> gain = 1
    >>> bfield = compute_bfield(
    ...     bdot_voltage,
    ...     time_array,
    ...     loop_area,
    ...     num_loop=num_loop,
    ...     gain=gain,
    ... )
    >>> print(bfield)
    """
    if len(bdot_voltage) != len(time_array):
        raise ValueError("sizes of time_array and bdot_voltage are not equal.")
    
    bdot_voltage = bdot_voltage.value
    loop_area = loop_area.value
    time_array = time_array.value
    bdot_voltage /= (loop_area*num_loop*gain)
    # units are removed by cumulative_trapezoid
    magnetic_field_array = cumulative_trapezoid(bdot_voltage, time_array, initial=0)  # tesla

    return magnetic_field_array * u.T
