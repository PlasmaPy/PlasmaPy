"""Functionality to calculate the conditional average and conditional variance of a time series."""

import numpy as np

class ConditionlEvents:
    pass

def cond_av(S, T, smin, smax=None, Sref=None, delta=None, window=False, print_verbose=True):
    """
    Use: cond_av(S, T, smin, smax=None, Sref=None, delta=None, window=False)
    Use the level crossing algorithm to compute the conditional average of
    a process.

    Inputs:
        S: Signal. Size N ................................. (1xN) np array
        T: Time base ...................................... (1xN) np array
        smin: Minimal peak amplitude
              in units of rms-value above mean value. ..... float
        smax: Maximal peak amplitude. ..................... float, def None
        Sref: Reference signal.
              If None, S is the reference. ................ (1xN) np array,
                                                            def None
        delta: The size of the conditionally averaged signal. If window = True,
               it is also the minimal distance between two peaks.
               If delta = None, it is estimated as
               delta = len(S)/(number of conditional events)*timestep.
               ............................................ float, def None
        window: If True, delta also gives the minimal distance between peaks.
                ........................................... bool, def False
        print_verbose: If True, information about conditional events is printed, 
                       If False, no information is printed. 
                       .................................... bool, def True 

    Outputs:
        Svals: Signal values used in the conditional average.
               S with unused values set to nan. ........... (1xN) np array
        s_av: conditionally averaged signal ............... np array
        s_var: conditional variance of events ............. np array
        t_av: time base of s_av ........................... np array
        peaks: max amplitudes of conditionally averaged events
        wait: waiting times between peaks
    """
    import numpy as np

    if Sref is None:
        Sref = S
    assert len(Sref) == len(S) and len(S) == len(T)

    sgnl = (Sref - np.mean(Sref)) / np.std(Sref)
    dt = sum(np.diff(T)) / (len(T) - 1)

    places = np.where(sgnl > smin)[0]
    assert len(places) > 0, "No conditional events"
    if print_verbose:
        print("places to check:{}".format(len(places)), flush=True)
    dplaces = np.diff(places)
    split = np.where(dplaces != 1)[0]
    # (+1 since dplaces is one ahead with respect to places)
    lT = np.split(places, split + 1)

    if delta is None:
        delta = len(sgnl) / len(lT) * dt
    if smax is not None:
        too_high = []
        for i in range(len(lT)):
            if max(sgnl[lT[i]]) > smax:
                too_high.append(i)
        removeset = set(too_high)
        lT = [v for i, v in enumerate(lT) if i not in removeset]

    # Use arange instead of linspace to guarantee 0 in the middle of the array.
    t_av = np.arange(-int(delta / (dt * 2)), int(delta / (dt * 2)) + 1) * dt

    # diagnostics
    lplmax = 0
    lpldiff = 0
    lplcount = 0

    gpl_array = np.array([])
    if window:
        lT.sort(key=lambda t: max(sgnl[t]), reverse=True)

    # Find local and global peak values
    for i in range(len(lT)):
        local_peak_loc = np.where(Sref[lT[i]] == max(Sref[lT[i]]))[0]

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

        gpl_array = np.append(gpl_array, local_peak_loc + lT[i][0])

    # Ensure distance delta between peaks.
    if window:
        index = 0
        while index < len(gpl_array):
            t_to_close = np.where(
                abs(gpl_array[index + 1 :] - gpl_array[index]) < int(delta / dt)
            )[0] + (index + 1)
            gpl_array = np.delete(gpl_array, t_to_close)
            removeset = set(t_to_close)
            lT = [v for i, v in enumerate(lT) if i not in removeset]
            index += 1
        lT.sort(key=lambda t: max(t))
        gpl_array.sort()

    gpl_array = gpl_array.astype(int)
    peaks = S[gpl_array]
    wait = np.append(np.array([T[0]]), T[gpl_array])
    wait = np.diff(wait)

    Svals = np.zeros(len(sgnl))
    Svals[:] = np.nan

    badcount = 0

    t_half_len = int((len(t_av) - 1) / 2)
    s_tmp = np.zeros([len(t_av), len(lT)])

    for i in range(len(lT)):
        global_peak_loc = gpl_array[i]

        # Find the average values and their variance
        low_ind = int(max(0, global_peak_loc - t_half_len))
        high_ind = int(min(len(sgnl), global_peak_loc + t_half_len + 1))
        tmp_sn = S[low_ind:high_ind]
        Svals[low_ind:high_ind] = S[low_ind:high_ind]
        if low_ind == 0:
            tmp_sn = np.append(np.zeros(-global_peak_loc + t_half_len), tmp_sn)
        if high_ind == len(S):
            tmp_sn = np.append(
                tmp_sn, np.zeros(global_peak_loc + t_half_len + 1 - len(S))
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
    # OBS: We return 1-CV = <f>^2/<f^2>.
    s_var = s_av ** 2 / np.mean(s_tmp ** 2, axis=1)
    if print_verbose:
        print("conditional events:{}".format(len(peaks)), flush=True)
        if badcount > 0:
            print("bursts where the recorded peak is not the largest:" + str(badcount))
        if lplcount > 0:
            print(
                "There were problems locating unique peaks in {0} bursts, {1:.3f} of all bursts".format(
                    lplcount, lplcount / len(lT)
                )
            )
            print(
                "Largest number of peaks in burst:{0}\nLargest difference (data points) between peaks:{1}".format(
                    lplmax, lpldiff
                )
            )
            print("In all cases, the first peak per burst was used.")

    return Svals, s_av, s_var, t_av, peaks, wait
