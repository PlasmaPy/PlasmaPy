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
from scipy import signal
import astropy.units as u
# class for magnetics


class Magnetics:
    """
    Group together values for magnetic analysis
    """

    def __init__(self):
        print("Class created")

    def bdot_field(self, bdot, tloop_area, time_s, unit_flag="Gauss"):
        """_summary_ : Function to compute magnetic field for bdots

        Args:
            Bdot: array - volts  
            tloop_area: float - meters squared 
            times_s: array - seconds   

        Returns:
            _type_: _description_
        """
        bdot[:] = [x / tloop_area for x in bdot]
        field_arr = sp.cumtrapz(bdot, time_s)  # Tesla

        if unit_flag.lower() == "tesla":
            return np.array(field_arr) * u.tesla
        else:
            return np.array(field_arr*1e4) * u.gauss  # Gauss

    def mean_sub(self, data_array, start_idx, end_idx):
        """_summary_

        Args:     
            data_array (_type_): _description_
            start_idx (_type_): _description_
            end_idx (_type_): _description_

        Returns:
            _type_: _description_
        """
        data_array = np.array(data_array)
        meansub = np.mean(data_array[start_idx:end_idx])
        data_array[:] = [x - meansub for x in data_array]
        return data_array

    def band_pass_filter(self, data_array, cutoff, sampling_frequency, order):
        """_summary_

        Args:
            data_array (_type_): _description_
            start_freq (_type_): _description_
            end_freq (_type_): _description_
            order (_type_): _description_

        Returns:
            _type_: _description_

        def butter_highpass(cutoff, fs, order=5):
            nyq = 0.5 * fs
            normal_cutoff = cutoff/nyq
            b, a = signal.butter(order, normal_cutoff, btype='highpass', analog=False)
            return b, a

        def butter_highpass_filter(data, cutoff, fs, order=5):
            b, a = butter_highpass(cutoff, fs, order=order)
            y = signal.filtfilt(b, a, data)
            return y
        """
        nyq = 0.5 * sampling_frequency
        normal_cutoff = cutoff/nyq
        b, a = signal.butter(order, normal_cutoff,
                             btype='highpass', analog=False)
        y = signal.filtfilt(b, a, data_array)
        return y
