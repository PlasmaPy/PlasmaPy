"""Functionality to calculate running moments of time series"""
import numpy as np


def running_mean(signal, radius):
    """
    Calculate the running mean of a sequence

    Parameters
    ----------
    signal: |array_like|
        signal to be averaged.

    radius: int
        window size is 2*radius+1

    Returns
    -------
    |array_like|
        running mean of signal with length len(signal) - 2*radius

    Raises
    ------
    'ValueError'
        if len(signal <= 2 * radius)

    Example
    -------
    >>> from plasmapy.analysis.time_series.running_moments import running_mean
    >>> running_mean([1,2,3,4], 1)
    array([2., 3.])
    """
    if len(signal) <= 2 * radius:
        raise ValueError("len(signal) must be bigger than 2*radius")

    window = 2 * radius + 1
    run_mean = np.cumsum(signal, dtype=float)
    run_mean[window:] = run_mean[window:] - run_mean[:-window]
    return run_mean[window - 1 :] / window


def running_moment(signal, radius, moment=1, time=None):
    """
    Calculate either the running mean, standard deviation, skewness or excess kurtosis of a sequence

    Parameters
    ----------
    signal: |array_like|
       signal to be averaged

    radius: int
        window size is 2*radius +1 for running mean and 4*radius + 1 for higher moments

    moment: int
        which running moment to compute
                1: running mean.
                2: running standard deviation.
                3: running skewness
                4: running excess kurtosis

    time: |array_like|, optional
        time base of signal

    Returns
    -------
    run_moment: |array_like|
        running moment of signal
        length is (len(signal) - 2*radius) for running mean
        length is (len(signal) = 4*signal) for higher moments

    time: |array_like|
        time base corresponding to run_moment if time != None

    Raises
    ------
    'ValueError'
        if moment not in (1,2,3,4)

    'ValueError'
        if signal and time don't have same length

    'ValueError
        if len(singal) <= 4*radius if moment > 1


    Notes
    -----
    The running rms divides by window, not (window-1)

    Example
    -------
    >>> from plasmapy.analysis.time_series.running_moments import running_moment
    >>> running_moment([1, 2, 3, 2, 1], 1, 4, [1, 2, 3, 4, 5])
    (array([3.]), [3])
    """
    if moment not in range(1, 5):
        raise ValueError("Only first four moments implemented")

    if time and (len(signal) != len(time)):
        raise ValueError("signal and time must have same length")

    if moment == 1:
        if time:
            time = time[radius:-radius]
        return running_mean(signal, radius), time

    if len(signal) <= 4 * radius:
        raise ValueError("len(signal) must be bigger than 4*radius for chosen moment")

    tmp = signal[radius:-radius] - running_mean(signal, radius)

    if moment == 2:
        run_moment = np.sqrt(running_mean(tmp**2, radius))

    if moment == 3:
        run_moment = (
            running_mean(tmp**3, radius) / running_mean(tmp**2, radius) ** 1.5
        )

    if moment == 4:
        run_moment = (
            running_mean(tmp**4, radius) / running_mean(tmp**2, radius) ** 2
        )

    if time:
        time = time[2 * radius : -2 * radius]

    return run_moment, time
