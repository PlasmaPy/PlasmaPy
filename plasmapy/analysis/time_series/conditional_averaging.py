"""Functionality to calculate the conditional average and conditional variance of a time series."""

__all__ = ["ConditionlEvents"]


import numpy as np


class ConditionlEvents:
    def __init__(
        self,
        signal,
        time,
        lower_threshold,
        upper_threshold=None,
        reference_signal=None,
        delta=None,
        window=False,
        print_verbose=True,
    ):
        if reference_signal is None:
            reference_signal = signal
        assert len(reference_signal) == len(signal) and len(signal) == len(time)

        normalizes_signal = (reference_signal - np.mean(reference_signal)) / np.std(
            reference_signal
        )
        time_step = sum(np.diff(time)) / (len(time) - 1)

        places = np.where(normalizes_signal > lower_threshold)[0]

        assert len(places) > 0, "No conditional events"

        if print_verbose:
            print("places to check:{}".format(len(places)), flush=True)

        distance_between_places = np.diff(places)
        _split = np.where(distance_between_places != 1)[0]
        # (+1 since dplaces is one ahead with respect to places)
        conditional_events = np.split(places, _split + 1)

        if delta is None:
            delta = len(normalizes_signal) / len(conditional_events) * time_step
        if upper_threshold is not None:
            too_high = []
            for i in range(len(conditional_events)):
                if max(normalizes_signal[conditional_events[i]]) > upper_threshold:
                    too_high.append(i)
            removeset = set(too_high)
            conditional_events = [v for i, v in enumerate(conditional_events) if i not in removeset]

        # Use arange instead of linspace to guarantee 0 in the middle of the array.
        t_av = (
            np.arange(-int(delta / (time_step * 2)), int(delta / (time_step * 2)) + 1)
            * time_step
        )

        # diagnostics
        lplmax = 0
        lpldiff = 0
        lplcount = 0

        gpl_array = np.array([])
        if window:
            conditional_events.sort(key=lambda t: max(normalizes_signal[t]), reverse=True)

        # Find local and global peak values
        for i in range(len(conditional_events)):
            local_peak_loc = np.where(
                reference_signal[conditional_events[i]] == max(reference_signal[conditional_events[i]])
            )[0]

            # Troubleshooting in case there are more than one unique peak
            if len(local_peak_loc) > 1:
                lplmax = max(lplmax, len(local_peak_loc))
                lpldiff = max(lpldiff, local_peak_loc[-1] - local_peak_loc[0])
                lplcount += 1
                # Prefer the peak closest to the mean of the candidates
                # (earliest time breaks ties)
                local_peak_loc = local_peak_loc[
                    abs(local_peak_loc - np.mean(local_peak_loc)).argmin()
                ]

            gpl_array = np.append(gpl_array, local_peak_loc + conditional_events[i][0])

        # Ensure distance delta between peaks.
        if window:
            index = 0
            while index < len(gpl_array):
                t_to_close = np.where(
                    abs(gpl_array[index + 1 :] - gpl_array[index])
                    < int(delta / time_step)
                )[0] + (index + 1)
                gpl_array = np.delete(gpl_array, t_to_close)
                removeset = set(t_to_close)
                conditional_events = [v for i, v in enumerate(conditional_events) if i not in removeset]
                index += 1
            conditional_events.sort(key=lambda t: max(t))
            gpl_array.sort()

        gpl_array = gpl_array.astype(int)
        peaks = signal[gpl_array]
        wait = np.append(np.array([time[0]]), time[gpl_array])
        wait = np.diff(wait)

        Svals = np.zeros(len(normalizes_signal))
        Svals[:] = np.nan

        badcount = 0

        t_half_len = int((len(t_av) - 1) / 2)
        s_tmp = np.zeros([len(t_av), len(conditional_events)])

        for i in range(len(conditional_events)):
            global_peak_loc = gpl_array[i]

            # Find the average values and their variance
            low_ind = int(max(0, global_peak_loc - t_half_len))
            high_ind = int(
                min(len(normalizes_signal), global_peak_loc + t_half_len + 1)
            )
            tmp_sn = signal[low_ind:high_ind]
            Svals[low_ind:high_ind] = signal[low_ind:high_ind]
            if low_ind == 0:
                tmp_sn = np.append(np.zeros(-global_peak_loc + t_half_len), tmp_sn)
            if high_ind == len(signal):
                tmp_sn = np.append(
                    tmp_sn, np.zeros(global_peak_loc + t_half_len + 1 - len(signal))
                )
            if max(tmp_sn) != tmp_sn[t_half_len]:
                badcount += 1

            s_tmp[:, i] = tmp_sn
        s_av = np.mean(s_tmp, axis=1)

        # The conditional variance of the conditional event f(t) is defined as
        # CV = <(f-<f>)^2>/<f^2> = 1 - <f>^2/<f^2>
        # at each time t.
        # For a highly reproducible signal, f~<f> and CV = 0.
        # For a completely random signal, <f^2> >> <f>^2 and CV = 1.
        # Note: We return 1-CV = <f>^2/<f^2>.
        s_var = s_av**2 / np.mean(s_tmp**2, axis=1)
        if print_verbose:
            print("conditional events:{}".format(len(peaks)), flush=True)
            if badcount > 0:
                print(
                    "bursts where the recorded peak is not the largest:" + str(badcount)
                )
            if lplcount > 0:
                print(
                    "There were problems locating unique peaks in {0} bursts, {1:.3f} of all bursts".format(
                        lplcount, lplcount / len(conditional_events)
                    )
                )
                print(
                    "Largest number of peaks in burst:{0}\nLargest difference (data points) between peaks:{1}".format(
                        lplmax, lpldiff
                    )
                )
                print("In all cases, the first peak per burst was used.")

        return Svals, s_av, s_var, t_av, peaks, wait

    @property
    def time(self):
        pass

    @property
    def average(self):
        pass

    @property
    def variance(self):
        pass

    @property
    def peaks(self):
        pass

    @property
    def waiting_times(self):
        pass

    @property
    def number_of_events(self):
        pass

    @property
    def number_of_unrecorded_events(self):
        pass
