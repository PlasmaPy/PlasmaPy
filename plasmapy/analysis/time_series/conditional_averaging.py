"""Functionality to calculate the conditional average and conditional variance of a time series."""

__all__ = ["ConditionalEvents"]


import numpy as np
from scipy.signal import find_peaks
from astropy import units as u


class ConditionalEvents:
    """
    Calculate conditional average, conditional variance, peaks,
    arrival times and waiting times of events of a time series.

    Parameters
    ----------
    signal : 1D |array_like|
        Signal to be analyzed.
    time : 1D |array_like|
        Corresponding time values for ``signal``.
    lower_threshold : `float` or `~astropy.units.Quantity`
        Lower threshold for event detection.
    upper_threshold : `float` or `~astropy.units.Quantity`, default: `None`
        Upper threshold for event detection.
    reference_signal : 1D |array_like|, default: `None`
        Reference signal.
        If `None`, ``signal`` is the reference signal.
    length_of_return : float, default, `None`
        Desired length of returned data.
        If `None`, estimated as ``len(signal) / len(number_of_events) * time_step``.
    distance : float, default: ``0``
        Minimum distance between peaks, in units of time.

    Raises
    ------
    `ValueError`:
        If length of ``signal`` and ``time`` are not equal.
        If length of ``reference_signal`` and ``time`` are not equal (when reference_signal is provided).
        If ``length_of_return`` is greater than the length of the time span.
        If ``length_of_return`` is negative.
        If ``upper_threshold`` is less than or equal to ``lower_threshold``.

    `TypeError`:
        If ``signal``/``reference_signal`` and ``lower_threshold`` have different astropy units.
        If ``signal``/``reference_signal`` and ``upper_threshold`` have different astropy units.

    Notes
    -----
    A detailed analysis of the conditional averaging method is presented in
    Rolf Nilsen's master thesis: "Conditional averaging of overlapping pulses"
    https://munin.uit.no/handle/10037/29416

    Example
    -------
    >>> from plasmapy.analysis.time_series.conditional_averaging import ConditionalEvents
    >>> cond_events = ConditionalEvents(signal = [1, 2, 1, 1, 2, 1], time = [1, 2, 3, 4, 5, 6], lower_threshold = 1.5)
    >>> cond_events.time
    [-1.0, 0.0, 1.0]
    >>> cond_events.average
    [1.0, 2.0, 1.0]
    >>> cond_events.variance
    [1.0, 1.0, 1.0]
    >>> cond_events.peaks
    [2, 2]
    >>> cond_events.waiting_times
    [3]
    >>> cond_events.arrival_times
    [2, 5]
    >>> cond_events.number_of_events
    2
    """

    def __init__(
        self,
        signal,
        time,
        lower_threshold,
        upper_threshold=None,
        reference_signal=None,
        length_of_return=None,
        distance=0,
    ):
        # This astropy unit checks are quite ugly in my view.
        # If a code reviewer has a better idea how to handle this I would be very grateful.
        if reference_signal is not None:
            self._check_units_consistency(
                reference_signal, lower_threshold, upper_threshold
            )
        else:
            self._check_units_consistency(signal, lower_threshold, upper_threshold)

        self._astropy_unit = None

        if isinstance(signal, u.Quantity):
            signal, self._astropy_unit = self._strip_unit_from_variable(signal)

        if isinstance(lower_threshold, u.Quantity):
            lower_threshold = lower_threshold.value

        if isinstance(upper_threshold, u.Quantity):
            upper_threshold = upper_threshold.value

        if isinstance(reference_signal, u.Quantity):
            reference_signal = reference_signal.value

        if distance < 0:
            raise ValueError("distance can't be negative")

        if len(signal) != len(time):
            raise ValueError("length of signal and time must be equal")

        if reference_signal is not None:
            if len(reference_signal) != len(time):
                raise ValueError("length of reference_signal and time must be equal")

        if length_of_return is not None:
            if length_of_return > time[-1] - time[0]:
                raise ValueError(
                    "choose length_of_return shorter or euqal to time length"
                )
            if length_of_return < 0:
                raise ValueError("length_of_return must be bigger than 0")

        if upper_threshold:
            if upper_threshold <= lower_threshold:
                raise ValueError(
                    "upper_threshold higher than lower_threshold, no events will be found"
                )

        if reference_signal is None:
            reference_signal = signal.copy()

        signal = self._ensure_numpy_array(signal)
        time = self._ensure_numpy_array(time)
        reference_signal = self._ensure_numpy_array(reference_signal)

        time_step = np.diff(time).sum() / (len(time) - 1)

        peak_locations, _ = find_peaks(
            reference_signal,
            height=[lower_threshold, upper_threshold],
            distance=int(distance / time_step) + 1,
        )

        conditional_events_indices = self._separate_events(
            reference_signal, lower_threshold, upper_threshold
        )

        peak_indices = self._choose_largest_peak_per_event(
            reference_signal, conditional_events_indices, peak_locations
        )

        if length_of_return is None:
            length_of_return = len(signal) / len(conditional_events_indices) * time_step

        self._return_time = list(
            (
                np.arange(
                    -int(length_of_return / (time_step * 2)),
                    int(length_of_return / (time_step * 2)) + 1,
                )
                * time_step
            )
        )

        self._peaks = list(signal[peak_indices])
        self._number_of_events = len(self._peaks)

        self._arrival_times = list(time[peak_indices])
        self._waiting_times = list(np.diff(self._arrival_times))

        self._conditional_average, conditional_events = self._average_over_events(
            signal, peak_indices
        )

        self._conditional_variance = list(
            self._calculate_conditional_variance(conditional_events)
        )

        self._conditional_average = list(self._conditional_average)

        if self._astropy_unit is not None:
            self._peaks *= self._astropy_unit
            self._conditional_average *= self._astropy_unit

    @property
    def time(self):
        """
        Time values corresponding to the analysis window.

        Returns
        -------
        time : 1D |array_like|
            Time values representing the analysis window.
        """
        return self._return_time

    @property
    def average(self):
        """
        Conditional average over events.

        Returns
        -------
        average : 1D |array_like|
            Array representing the conditional average over events.

        """
        return self._conditional_average

    @property
    def variance(self):
        """
        Conditional variance over events.

        Returns
        -------
        variance : 1D |array_like|
            Array representing the conditional variance over events.

        """
        return self._conditional_variance

    @property
    def peaks(self):
        """
        Peak values of conditional events.

        Returns
        -------
        peaks : 1D |array_like|
            Peak values of conditional events.

        """
        return self._peaks

    @property
    def waiting_times(self):
        """
        Waiting times between consecutive peaks.

        Returns
        -------
        waiting_times : 1D |array_like|
            Waiting times between consecutive peaks of conditional events.

        """
        return self._waiting_times

    @property
    def arrival_times(self):
        """
        Arrival times corresponding to the conditional events.

        Returns
        -------
        arrival_times : 1D |array_like|
            Arrival times the conditional events.

        """
        return self._arrival_times

    @property
    def number_of_events(self):
        """
        Total number of conditional events.

        Returns
        -------
        number_of_events : int
            Total number of conditional events.

        """
        return self._number_of_events

    # This astropy unit checks are quite ugly in my view.
    # If a code reviewer has a better idea how to handle this I would be very grateful.
    def _check_units_consistency(self, signal, lower_threshold, upper_threshold):
        if isinstance(signal, u.Quantity):
            if not isinstance(lower_threshold, u.Quantity):
                raise TypeError(
                    "signal/reference_signal and lower_threshold must have same astropy unit"
                )
            if signal.unit != lower_threshold.unit:
                raise TypeError(
                    "signal/reference_signal and lower_threshold must have same astropy unit"
                )
            if upper_threshold is not None:
                if not isinstance(upper_threshold, u.Quantity):
                    raise TypeError(
                        "signal/reference_signal and lower_threshold must have same astropy unit"
                    )
                if signal.unit != upper_threshold.unit:
                    raise TypeError(
                        "signal/reference_signal and upper_threshold must have same astropy unit"
                    )

    def _strip_unit_from_variable(self, variable):
        return variable.value, variable.unit

    def _ensure_numpy_array(self, variable):
        if not isinstance(variable, np.ndarray):
            variable = np.array(variable)
        return variable

    def _separate_events(self, reference_signal, lower_threshold, upper_threshold):
        places = np.where(reference_signal > lower_threshold)[0]
        if upper_threshold:
            higher = np.where(reference_signal < upper_threshold)[0]
            places = np.intersect1d(places, higher)

        distance_between_places = np.diff(places)
        _split = np.where(distance_between_places != 1)[0]
        return np.split(places, _split + 1)

    def _choose_largest_peak_per_event(
        self, reference_signal, conditional_events_indices, peak_indices
    ):
        for event in conditional_events_indices:
            peaks_in_event = np.isin(peak_indices, event)
            if peaks_in_event.sum() > 1:
                peak_ind = peak_indices[peaks_in_event]
                highest_local_peak = reference_signal[peak_ind].argmax()
                not_highest_local_peaks = np.delete(peak_ind, highest_local_peak)
                peak_indices = np.delete(
                    peak_indices, np.isin(peak_indices, not_highest_local_peaks)
                )
        return peak_indices

    def _average_over_events(self, signal, peak_indices):

        t_half_len = int((len(self._return_time) - 1) / 2)
        conditional_events = np.zeros([len(self._return_time), len(peak_indices)])

        for i, global_peak_loc in enumerate(peak_indices):
            low_ind = int(max(0, global_peak_loc - t_half_len))
            high_ind = int(min(len(signal), global_peak_loc + t_half_len + 1))
            single_event = signal[low_ind:high_ind]
            if low_ind == 0:
                single_event = np.append(
                    np.zeros(-global_peak_loc + t_half_len), single_event
                )
            if high_ind == len(signal):
                single_event = np.append(
                    single_event,
                    np.zeros(global_peak_loc + t_half_len + 1 - len(signal)),
                )

            conditional_events[:, i] = single_event

        return np.mean(conditional_events, axis=1), conditional_events

    def _calculate_conditional_variance(self, conditional_events):
        return self._conditional_average**2 / np.mean(conditional_events**2, axis=1)
