"""Functionality to calculate the conditional average and conditional variance of a time series."""

__all__ = ["ConditionlEvents"]


import numpy as np
from scipy.signal import find_peaks


class ConditionlEvents:
    def __init__(
        self,
        signal,
        time,
        lower_threshold,
        upper_threshold=None,
        reference_signal=None,
        length_of_return=None,
        distance=0,
        weight="amplitude",
    ):
        assert weight in ("amplitude", "equal")
        assert distance >= 0, "distance can't be negative"
        if length_of_return:
            assert (
                length_of_return < time[-1] - time[0]
            ), "choose length_of_return shorter than time length"
            assert length_of_return > 0, "length_of_return must be positive"

        if upper_threshold:
            assert (
                upper_threshold > lower_threshold
            ), "upper_threshold higher than lower_threshold, no events will be found"

        if reference_signal is None:
            reference_signal = signal.copy()

        assert len(reference_signal) == len(signal) and len(signal) == len(time)

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

        self._return_time = (
            np.arange(
                -int(length_of_return / (time_step * 2)),
                int(length_of_return / (time_step * 2)) + 1,
            )
            * time_step
        )

        self._peaks = signal[peak_indices]
        self._number_of_events = len(self._peaks)

        self._arrival_times = time[peak_indices]
        self._waiting_times = np.diff(
            np.append(np.array([time[0]]), self._arrival_times)
        )

        self._conditional_average, conditional_events = self._average_over_events(
            signal, weight, peak_indices
        )

        self._conditional_variance = self._calculate_conditional_variance(
            conditional_events
        )

    @property
    def time(self):
        return self._return_time

    @property
    def average(self):
        return self._conditional_average

    @property
    def variance(self):
        return self._conditional_variance

    @property
    def peaks(self):
        return self._peaks

    @property
    def waiting_times(self):
        return self._waiting_times

    @property
    def arrival_times(self):
        return self._arrival_times

    @property
    def number_of_events(self):
        return self._number_of_events

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
        for i, event in enumerate(conditional_events_indices):
            peaks_in_event = np.isin(peak_indices, event)
            if peaks_in_event.sum() > 1:
                peak_ind = peak_indices[peaks_in_event]
                highest_local_peak = reference_signal[peak_ind].argmax()
                not_highest_local_peaks = np.delete(peak_ind, highest_local_peak)
                peak_indices = np.delete(
                    peak_indices, np.isin(peak_indices, not_highest_local_peaks)
                )
        return peak_indices

    def _average_over_events(self, signal, weight, peak_indices):

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
            if weight == "equal":
                conditional_events[:, i] /= single_event[t_half_len]

        return np.mean(conditional_events, axis=1), conditional_events

    def _calculate_conditional_variance(self, conditional_events):
        # The conditional variance of the conditional event f(t) is defined as
        # CV = <(f-<f>)^2>/<f^2> = 1 - <f>^2/<f^2>
        # at each time t.
        # For a highly reproducible signal, f~<f> and CV = 0.
        # For a completely random signal, <f^2> >> <f>^2 and CV = 1.
        # HINT: We return 1-CV = <f>^2/<f^2>.
        return self._conditional_average**2 / np.mean(conditional_events**2, axis=1)
