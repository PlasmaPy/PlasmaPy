"""Functionality to calculate running moments of time series"""
import numpy as np


def running_mean(signal, radius):
    """
    Calculate the running mean of a sequence

    Parameters
    ----------
    signal: |array_like|
        signal to be averagedself.

    radius: int
        window size is 2*radius+1

    Returns
    -------
    |array_like|
        running mean of signal with length len(signal) - 2*radius

    Example
    -------
    >>> from plasmapy.analysis.time_series.running_moments import running_mean
    >>> running_mean([1,2,3,4], 1)
    array([2., 3.])
    """
    window = 2 * radius + 1
    run_mean = np.cumsum(signal, dtype=float)
    run_mean[window:] = run_mean[window:] - run_mean[:-window]
    return run_mean[window - 1 :] / window
