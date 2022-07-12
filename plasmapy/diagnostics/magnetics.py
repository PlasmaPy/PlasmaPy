"""
Defines the magnetics analysis module as part of `plasmapy.diagnostics`.
"""

__all__ = [
    
]

#import Packages
from warnings import warn
import scipy.io as spio
import matplotlib.pylab as plt
import numpy as np
import scipy.integrate as sp
import astropy.units as u
# class for magnetics
class Magnetics:
    """
    Group together values for magnetic analysis
    """

    def __init__(self):
        print("Class created")

    def bdot_field(self, bdot, tloop_area, time_s, unit_flag = "Gauss"):
        """
        Function to Compute Magnetic Field for Bdots
        @param Bdot: array - volts  
               tloop_area: float - meters squared 
               times_s: array - seconds 
               
        @return magnet field of bdot in gaus 
        """
        bdot[:] = [x / tloop_area for x in bdot]
        field_arr = sp.cumtrapz(bdot , time_s) #Tesla
        
        if unit_flag.lower() == "tesla":
            return np.array(field_arr) * u.tesla
        else:
            return np.array(field_arr*1e4) * u.gauss  # Gauss  