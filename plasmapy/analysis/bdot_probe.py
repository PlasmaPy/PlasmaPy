"""_summary_
   
"""
__all__ = []

from warnings import warn
import scipy.integrate as sp
import matplotlib.pylab as plt
from typing import Tuple
import numpy as np


def compute_bfield(bdot_data: np.ndarray, times: np.ndarray, tloop_area: float) -> Tuple[np.ndarray, np.ndarray]:
    """_summary_ : Function to compute magnetic field for bdots
         Args:
             Bdot (array): - volts  
             tloop_area (float): - meters squared 
             times_s (array): - seconds   
         Returns:
             field_arr : array of magnetic field values for the bdots
             new_time: corrsponding time array  

    """
    if len(bdot_data) != len(times):
        raise Exception("length of time and voltage arrays in not equal\n")

    bdot_data[:] = [x / tloop_area for x in bdot_data]
    field_arr = sp.cumtrapz(bdot_data, times)  # Tesla
    new_time = times[1:]

    return field_arr, new_time
