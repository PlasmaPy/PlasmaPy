"""Functionality to calculate the conditional average and conditional variance of a time series."""

__all__ = ["ConditionlEvents"]


import numpy as np
from scipy.signal import find_peaks, peak_prominences

class ConditionlEvents:
    def __init__(
        self, S, T, smin=None, smax=None, Sref=None, delta=None, window=False,  prominence=None, threshold = None, weight='amplitude'
    ):
    
    if all(i is None for i in[smin, prominence]):
        raise TypeError('Missing 1 required positional argument: \'smin\' '
                        'or \'prominence\'')
    
    if Sref is None:
        Sref = S.copy()
    assert len(Sref) == len(S) and len(S) == len(T)

    sgnl = (Sref - np.mean(Sref)) / np.std(Sref)
    dt = np.diff(T).sum() / (len(T) - 1)
    
    #Estimating delta.
    if delta is None:
        if smin is None:
            tmpmin = prominence
        else:
            tmpmin = smin
            
        places = np.where(sgnl > tmpmin)[0]
        dplaces = np.diff(places)
        split = np.where(dplaces != 1)[0]
        # (+1 since dplaces is one ahead with respect to places)
        lT = np.split(places, split + 1)
        delta = len(sgnl) / len(lT) * dt
    
    distance = None
    # Ensure distance delta between peaks.
    if window:
        distance = int(delta / dt)
        
    
    
    # Find peak indices.
    gpl_array, properties = find_peaks(sgnl, height = [smin, smax], distance = distance,
                              prominence = prominence, threshold = threshold)
    
    if not prominence and not window:
        places = np.where(sgnl > smin)[0]
        if smax:
            higher = np.where(sgnl<smax)[0]
            places = np.intersect1d(places, higher)
        assert len(places) > 0, "No conditional events"
   
        dplaces = np.diff(places)
        split = np.where(dplaces != 1)[0]
        lT = np.split(places, split + 1)
        #Ensure one peak for every threshold crossing. Largest peak is chosen.
        #Ties are broken my earliest time.
        for i, lTi in enumerate(lT): 
            peak_check = np.isin(gpl_array, lTi)
            if peak_check.sum()>1:
                peak_ind = gpl_array[peak_check]
                highest_local_peak = Sref[peak_ind].argmax()
                not_highest_local_peaks = np.delete(peak_ind, highest_local_peak)
                gpl_array = np.delete(gpl_array, np.isin(gpl_array, not_highest_local_peaks))
        
    
    # Use arange instead of linspace to guarantee 0 in the middle of the array.
    t_av = np.arange(-int(delta / (dt * 2)), int(delta / (dt * 2)) + 1) * dt

    
    peaks = S[gpl_array]
    prominences = peak_prominences(S, gpl_array)[0]
    wait = np.append(np.array([T[0]]), T[gpl_array])
    wait = np.diff(wait)

    Svals = np.zeros(len(sgnl))
    Svals[:] = np.nan

    badcount = 0

    t_half_len = int((len(t_av) - 1) / 2)
    s_tmp = np.zeros([len(t_av), len(gpl_array)])

    for i, global_peak_loc in enumerate(gpl_array):

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
        if tmp_sn.max() != tmp_sn[t_half_len]:
            badcount += 1
        
        s_tmp[:, i] = tmp_sn
        if weight == 'equal':
            if smin is None:
                s_tmp[:, i] /= prominences[i]
                s_tmp[:, i] += (1-s_tmp[t_half_len, i])
                
            else:
                s_tmp[:, i] /= tmp_sn[t_half_len]
    s_av = np.mean(s_tmp, axis=1)

    # The conditional variance of the conditional event f(t) is defined as
    # CV = <(f-<f>)^2>/<f^2> = 1 - <f>^2/<f^2>
    # at each time t.
    # For a highly reproducible signal, f~<f> and CV = 0.
    # For a completely random signal, <f^2> >> <f>^2 and CV = 1.
    # OBS: We return 1-CV = <f>^2/<f^2>.
    s_var = s_av ** 2 / np.mean(s_tmp ** 2, axis=1)
    print("conditional events:{}".format(len(peaks)), flush=True)
    if badcount > 0:
        print("bursts where the recorded peak is not the largest:" + str(badcount))

    self.return_time = t_av
    self.s_av = s_av
    # return Svals, s_av, s_var, t_av, peaks, wait, prominences
i
    @property
    def time(self):
        return self.return_time

    @property
    def average(self):
        return self.s_av

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
    def number_of_unused_peaks(self):
        pass
