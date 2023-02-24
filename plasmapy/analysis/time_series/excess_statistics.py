__all__ = ["excess_stat"]


import numbers
import numpy as np

from collections import namedtuple

_excess_stat_tuple = namedtuple(
    "Excess_Statistics",
    [
        "total_time_above_threshold",
        "number_of_crossings",
        "average_times",
        "rms_times",
        "hist",
        "bin_centers",
    ],
)


def excess_stat(signal, thresholds, time_step, pdf=False, bins=32):
    """
    For a given signal signal and a given threshold values a in thresholds, this function finds the excess statitics of signal.
    Input:
        signal: the signal, a 1d numpy array (1d np)
        thresholds: the threshold values, 1d np
        time_step: the time resolution, float
        pdf: if set to True, the function estimates the PDF of the time above threshold as well.
        bins: The number of bins in the estimation of the PDF above threshold.
    Output:
        Theta_array: The total time above threshold for each value in thresholds, 1d np.
        number_of_crossings: The total number of upwards crossings over the threshold for each value in thresholds, 1d np.
        avT_array: The average time above threshold for each value in thresholds, 1d np.
        rmsT_array: The rms-value of time above threshold for each value in thresholds, 1d np.
        If pdf is set to True, we additionaly output
        Tpdf: thresholds 2d numpy array, shape (thresholds.size,bins=32). For each value in thresholds, this is the estimated PDF of time above threshold.
        t: Time values for Tpdf, same shape.
    """
    if not isinstance(bins, numbers.Integral):
        raise TypeError("bins must be of type integer")

    if time_step <= 0:
        raise ValueError("time_step must be positive")

    total_time_above_threshold = []
    number_of_crossings = []
    average_times = []
    rms_times = []
    hist = None
    bin_centers = None

    if pdf:
        events_per_threshold = {}

    for threshold in thresholds:
        indices_above_threshold = np.where(np.array(signal) > threshold)[0]

        if len(indices_above_threshold) == 0:
            times_above_threshold = []
            total_time_above_threshold.append(0)
            number_of_crossings.append(0)
            average_times.append(0)
            rms_times.append(0)

        else:
            total_time_above_threshold.append(time_step * len(indices_above_threshold))

            distances_to_next_index = (
                indices_above_threshold[1:] - indices_above_threshold[:-1]
            )
            split_indices = np.where(distances_to_next_index != 1)[0]
            event_lengths = np.split(distances_to_next_index, split_indices)

            # set correct length for first event
            event_lengths[0] = np.append(event_lengths[0], 1)

            times_above_threshold = [
                time_step * len(event_lengths[i]) for i in range(0, len(event_lengths))
            ]

            number_of_crossings.append(len(event_lengths))
            if indices_above_threshold[0] == 0:
                # Don't count the first event if there is no crossing.
                number_of_crossings[-1] -= 1

            average_times.append(np.mean(times_above_threshold))
            rms_times.append(np.std(times_above_threshold))

        if pdf:
            events_per_threshold.update({threshold: times_above_threshold})

    if pdf:
        hist, bin_centers = calculate_pdfs(events_per_threshold, bins)

    return _excess_stat_tuple(
        total_time_above_threshold=total_time_above_threshold,
        number_of_crossings=number_of_crossings,
        average_times=average_times,
        rms_times=rms_times,
        hist=hist,
        bin_centers=bin_centers,
    )


def calculate_pdfs(events_per_threshold, bins):
    """
    Calculate the pdf P(time_step|thresholds) and avT from this pdf.
    time_step: dictionary. From above
    Returns the 2d time array t, and the 2d-array time_steppdf, containing the pdfs.
    t and time_steppdf are both 2d-arrays storing the values for each a along the axis. The pdf for thresholds[i] is time_steppdf[i,:], t[i,:].
    """
    hist = np.zeros((len(events_per_threshold), bins))
    bin_centers = np.zeros((len(events_per_threshold), bins))

    for i, threshold in enumerate(events_per_threshold.keys()):
        if len(events_per_threshold[threshold]) >= 1:
            hist[i, :], bin_edges = np.histogram(
                events_per_threshold[threshold], bins=bins, density=True
            )
            bin_centers[i, :] = (bin_edges[1:] + bin_edges[:-1]) / 2
    return hist, bin_centers
