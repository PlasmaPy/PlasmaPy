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
        size_of_return=None,
        distance=1,  # minimum distance in index
        weight="amplitude",
    ):

        if reference_signal is None:
            reference_signal = signal.copy()
        assert len(reference_signal) == len(signal) and len(signal) == len(time)

        normalized_signal = (reference_signal - np.mean(reference_signal)) / np.std(
            reference_signal
        )
        time_step = np.diff(time).sum() / (len(time) - 1)

        # Find peak indices.
        self.peak_indices, _ = find_peaks(
            normalized_signal,
            height=[lower_threshold, upper_threshold],
            distance=distance,
        )

        places = np.where(normalized_signal > lower_threshold)[0]
        if upper_threshold:
            higher = np.where(normalized_signal < upper_threshold)[0]
            places = np.intersect1d(places, higher)

        distance_between_places = np.diff(places)
        _split = np.where(distance_between_places != 1)[0]
        self.conditional_events = np.split(places, _split + 1)


        self.peak_indices = self._choose_largest_peak_per_event(reference_signal, self.conditional_events, self.peak_indices)

        if size_of_return is None:
            size_of_return = len(normalized_signal) / len(self.conditional_events) * time_step

        return_time = (
            np.arange(
                -int(size_of_return / (time_step * 2)),
                int(size_of_return / (time_step * 2)) + 1,
            )
            * time_step
        )

        peaks = signal[self.peak_indices]
        waiting_times = np.append(np.array([time[0]]), time[self.peak_indices])
        waiting_times = np.diff(waiting_times)

        Svals = np.zeros(len(normalized_signal))
        Svals[:] = np.nan

        badcount = 0

        t_half_len = int((len(return_time) - 1) / 2)
        s_tmp = np.zeros([len(return_time), len(self.peak_indices)])

        for i, global_peak_loc in enumerate(self.peak_indices):

            # Find the average values and their variance
            low_ind = int(max(0, global_peak_loc - t_half_len))
            high_ind = int(
                min(len(normalized_signal), global_peak_loc + t_half_len + 1)
            )
            tmp_sn = signal[low_ind:high_ind]
            Svals[low_ind:high_ind] = signal[low_ind:high_ind]
            if low_ind == 0:
                tmp_sn = np.append(np.zeros(-global_peak_loc + t_half_len), tmp_sn)
            if high_ind == len(signal):
                tmp_sn = np.append(
                    tmp_sn, np.zeros(global_peak_loc + t_half_len + 1 - len(signal))
                )
            if tmp_sn.max() != tmp_sn[t_half_len]:
                badcount += 1

            s_tmp[:, i] = tmp_sn
            if weight == "equal":
                s_tmp[:, i] /= tmp_sn[t_half_len]

        conditional_average = np.mean(s_tmp, axis=1)

        # The conditional variance of the conditional event f(t) is defined as
        # CV = <(f-<f>)^2>/<f^2> = 1 - <f>^2/<f^2>
        # at each time t.
        # For a highly reproducible signal, f~<f> and CV = 0.
        # For a completely random signal, <f^2> >> <f>^2 and CV = 1.
        # HINT: We return 1-CV = <f>^2/<f^2>.
        conditional_variance = conditional_average**2 / np.mean(s_tmp**2, axis=1)

        self._return_time = return_time
        self._conditional_average = conditional_average
        self._conditional_variance = conditional_variance
        self._waiting_times = waiting_times
        self._peaks = peaks
        self._number_of_events = peaks
        self._number_of_unused_peaks = badcount
        # return Svals, s_av, s_var, t_av, peaks, wait

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
    def number_of_events(self):
        return self._number_of_events

    @property
    def number_of_unused_peaks(self):
        return self._number_of_unused_peaks


    def _choose_largest_peak_per_event(self, reference_signal, conditional_events, peak_indices):
        for i, event in enumerate(conditional_events):
            peaks_in_event = np.isin(peak_indices, event)
            if peaks_in_event.sum() > 1:
                peak_ind = peak_indices[peaks_in_event]
                highest_local_peak = reference_signal[peak_ind].argmax()
                not_highest_local_peaks = np.delete(peak_ind, highest_local_peak)
                peak_indices = np.delete(
                    peak_indices, np.isin(peak_indices, not_highest_local_peaks)
                )
        return peak_indices
