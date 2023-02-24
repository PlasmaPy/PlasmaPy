"""Functionality to calculate excess statistics of time series."""

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
    Calculate total time, number of upwards crossings, average time and rms time
    above given thresholds of a sequence.

    Parameters
    ----------
    signal : 1D |array_like|
        Signal to be analyzed.

    thresholds : 1D |array_like|
        Threshold values.

    time_step : int
        Time step of ``signal``.

    pdf : bool, optional
        If set to True, the function estimates the PDF of the time above
        ``thresholds`` as well.

    bins : int, optional
        The number of bins in the estimation of the PDF above ``thresholds``.

    Returns
    -------
    `~collections.namedtuple`
        Contains the following attributes:

        - ``total_time_above_threshold``: 1D |array_like|
            Total time above threshold for each value in ``thresholds``.

        - ``number_of_crossings``: 1D |array_like|
            Total number of upwards crossings for each value in ``thresholds``.

        - ``average_times``: 1D |array_like|
            Average time above each value in ``thresholds``.

        - ``rms_times``: 1D |array_like|
            Root-mean-square values of time above each value in ``thresholds``.

        - ``hist``: 2D `~numpy.ndarray`, shape (``thresholds.size``, ``bins`` )
            For each value in ``thresholds``, returns the estimated PDF of time
            above threshold. Only returned if ``pdf`` set to True.

        - ``bin_centers``: 2D `~numpy.ndarray`, shape (``thresholds.size``, ``bins`` )
            Bin centers for ``hist``. Only returned if ``pdf`` set to True.

    Raises
    ------
    `TypeError`
        If ``bins`` is not of type `int`.

    `ValueError`
        If ``time_step`` <= 0.

    Example
    -------
    >>> from plasmapy.analysis.time_series.excess_statistics import excess_stat
    >>> excess_stat([0, 0, 2, 2, 0, 4], [1, 3, 5], 1)
    Excess_Statistics(total_time_above_threshold=[3, 1, 0],
                      number_of_crossings=[2, 1, 0],
                      average_times=[1.5, 1.0, 0],
                      rms_times=[0.5, 0.0, 0],
                      hist=None,
                      bin_centers=None)
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
        hist, bin_centers = _calculate_pdfs(events_per_threshold, bins)

    return _excess_stat_tuple(
        total_time_above_threshold=total_time_above_threshold,
        number_of_crossings=number_of_crossings,
        average_times=average_times,
        rms_times=rms_times,
        hist=hist,
        bin_centers=bin_centers,
    )


def _calculate_pdfs(events_per_threshold, bins):
    """
    Helper function for ``excess_stat``.
    Calculates the pdf P(time_step| thresholds).
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
