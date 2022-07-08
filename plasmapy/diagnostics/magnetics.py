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

# class for magnetics
class Magnetics:
    """
    Group together values for magnetic analysis
    """

    def __init__(self):
        print("Class created")

    def bdot_field(self, bdot, tloop_area, time_s):
        """
        Function to Compute Magnetic Field for Bdots
        @param Bdot: array, tloop_area: float, times_s: array
        @return
        """
        bdot[:] = [x / tloop_area for x in bdot]
        return sp.cumtrapz(bdot , time_s)*1e4  # Gauss 