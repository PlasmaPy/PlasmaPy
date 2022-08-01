"""_summary_
"""

import numpy as np
import scipy.io as spio
import scipy.integrate as sp
from scipy import signal
from warnings import warn


def remove_offset(data_array: np.ndarray, start_idx: int, end_idx: int) -> np.ndarray:
    """_summary_ : Function which preforms mean subtraction

    Args:     
       data_array (array): array containing values from which to mean sub
       start_idx (int): idx of where to start obtaining mean from
       end_idx (int): idx of where to stop obtaining mean from
    Returns:
       data_array (array): array after mean_sub is preformed 
    """
    data_array = np.array(data_array)
    meansub = np.mean(data_array[start_idx:end_idx])
    data_array[:] = [x - meansub for x in data_array]
    return data_array


def high_pass_filter(data_array: np.ndarray, cutoff: float, sampling_frequency: float, order: int) -> np.ndarray:
    """_summary_ : Function which applies a band_pass_filter

    Args:
       data_array (array): _description_
       cuttoff (_type_): _description_
       sampling_frequency (_type_): _description_
       order (_type_): _description_
    Returns:
       filtered_arr : array with filter applied 
    """
    if cutoff > sampling_frequency:
        raise ValueError(
            "sampling frequency must be greater than cutoff value\n")

    nyq = 0.5 * sampling_frequency
    normal_cutoff = cutoff/nyq
    b, a = signal.butter(order, normal_cutoff, btype='highpass', analog=False)
    filt_arr = signal.filtfilt(b, a, data_array)
    return filt_arr
