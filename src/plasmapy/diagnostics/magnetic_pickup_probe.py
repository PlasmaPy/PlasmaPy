"""Defines the basic magnetic pickup probe (Bdot probe) analysis
   
"""
__all__ = []

from warnings import warn
import scipy.integrate as sp
from typing import Tuple
import numpy as np


import astropy.units as u
from plasmapy.utils.decorators import validate_quantities
from typing import Tuple

@validate_quantities
def compute_bfield(bdot_data: np.ndarray, times: np.ndarray, loop_area: float) -> np.ndarray:
    r"""
    Takes the voltage output of a magnetic pickup probe (Bdot probe) and time base and 
    returns the associated array of the magnetic field as a function of time

    Parameters
    ----------
    bdot_data: `numpy.ndarray`
        A data array containing voltage fluctuations from a bdot probe. The
        array values are proportional to time changing magnetic field via Faradays Law
        The input values are assumed to be in volts.

    times: `numpy.ndarray`
        The time series to the data collection in units of seconds.

    loop_area: float
        The area through which the changing flux is measured in units of meters squared.

    Returns
    -------
    field_arr: `numpy.ndarray`
        The array containing the magnetic field fluctuations. 
        
    Notes
    -----
    The integral form of Faraday's Law states that the time change in flux through an area generates a voltage around a loop,
    .. math::
       V = -N\frac{d\Phi}{dt}
    where :math:`\Phi` is the magnetic flux
    .. math::
        \Phi = \int{BdA}
    and N is the number of loops.
    A magnetic pickup probe works by providing a fixed area loop or loops using a length of conducting wire through 
    which magnetic fields penetrate. As these fields change in time, a voltage is generated that can be measured
    by some kind of data aquisition device (such as an ossciloscope). We can also assume that the loop is small 
    enough that the magnetic field through the loop is constant everywhere in the loop (and only changes in time).
    Given this, we can write
    .. math::
        V = -A\frac{dB}{dt}
    where the voltage becomes preportional only to the change in magnetic field and the area of the loop. Magnetic
    field as a function of time can be recovered by time integrating :math:`V(t)/A`.
    """
    if len(bdot_data) != len(times):
        raise Exception("length of time and voltage arrays in not equal\n")

    bdot_data /= loop_area
    field_arr = sp.cumtrapz(bdot_data, times, initial=0)  # Tesla

    return field_arr

class Bdot:
    """Class representing a basic magnetic pickup probe (Bdot probe)

    """
    _settings = {}
    _signal = None
    _time = None
    _area = None

    @validate_quantities
    def __init__(
        self,
        signal: u.volt,
        time: u.s,
        area: u.m**2,
        *,
        nloops,
        gain,
    ):
        self._signal = signal.value
        self._time = time.value
        self._area = area.value
        self._settings = {
            "nloops": nloops,
            "gain": gain,
        }

    @validate_quantities
    def calc_bfield(self) -> u.T:
        """_summary_

        """
        bfield = compute_bfield(
            self._signal,
            self._time,
            self._area,
            #    ** self._settings,
        )
        return bfield * u.T