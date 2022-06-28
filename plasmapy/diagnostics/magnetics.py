"""
Defines the magnetics analysis module as part of `plasmapy.diagnostics`.
"""

__all__ = [
    "Characterstic"
]

#import Packages
from warnings import warn

# class for magnetics
import scipy.io as spio
import matplotlib.pylab as plt
import numpy as np
import scipy.integrate as sp


class Magnetics:
    """
    Group together values for magnetic analysis
    """

    def __init__(self):
        print("Class created")

    # Function to Compute Magnetic Field for Bdots
    def Bdot_field(Bdot5t, tloop_area, time_s):
        return sp.cumtrapz(Bdot5t/tloop_area, time_s)*1e4  # Gauss
