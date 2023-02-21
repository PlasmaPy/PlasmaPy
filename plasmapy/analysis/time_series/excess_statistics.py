import numpy as np


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
        X_array: The total number of upwards crossings over the threshold for each value in thresholds, 1d np.
        avT_array: The average time above threshold for each value in thresholds, 1d np.
        rmsT_array: The rms-value of time above threshold for each value in thresholds, 1d np.
        If pdf is set to True, we additionaly output
        Tpdf: thresholds 2d numpy array, shape (thresholds.size,bins=32). For each value in thresholds, this is the estimated PDF of time above threshold.
        t: Time values for Tpdf, same shape.
    """
    Theta_array = np.array([])
    X_array = np.array([])
    avT_array = np.array([])
    rmsT_array = np.array([])
    if pdf:
        time_step_dict = {}
    for a in thresholds:
        # This is the basis: the parts of the signal that are above the
        # threshold.
        places = np.where(signal > a)[0]
        if len(places) > 0:
            # print('binsum, places to check:{}'.format(len(places)))
            Theta = time_step * len(places)

            # Find X, avT an distribution of time_step
            # Each blob is connected, so discrete blobs have more than one time
            # length between them
            dplaces = places[1:] - places[:-1]
            split = np.where(dplaces != 1)[0]  # split the array
            # at places where distance between points is greater than one
            lT = np.split(dplaces, split)
            lT[0] = np.append(lT[0], 1)  # To get correct length of the first
            # binsumber of upwards crossings is equal to number of discrete blobs
            X = len(lT)
            if places[0] == 0:
                # Don't count the first blob if there is no crossing.
                X += -1
            time_step = np.array(
                [time_step * len(lT[i]) for i in range(0, len(lT))]
            )  # thresholdsrray of excess times
            avT = np.mean(time_step)
            rmsT = np.std(time_step)
        elif len(places) == 0:
            Theta = 0
            X = 0
            avT = 0
            rmsT = 0
            time_step = np.array([])
        Theta_array = np.append(Theta_array, Theta)
        X_array = np.append(X_array, X)
        avT_array = np.append(avT_array, avT)
        rmsT_array = np.append(rmsT_array, rmsT)
        if pdf:
            time_step_dict.update({a: time_step})

    if pdf:

        time_steppdf, t = Th_time_step(time_step_dict, thresholds, bins)

        return Theta_array, X_array, avT_array, rmsT_array, time_steppdf, t
    else:
        return Theta_array, X_array, avT_array, rmsT_array


def Th_time_step(time_step, thresholds, bins):
    """
    Calculate the pdf P(time_step|thresholds) and avT from this pdf.
    time_step: dictionary. From above
    Returns the 2d time array t, and the 2d-array time_steppdf, containing the pdfs.
    t and time_steppdf are both 2d-arrays storing the values for each a along the axis. The pdf for thresholds[i] is time_steppdf[i,:], t[i,:].
    """
    time_steppdf = np.zeros((len(thresholds), bins))
    t = np.zeros((len(thresholds), bins))

    for i in range(0, len(thresholds)):
        a = thresholds[i]
        if len(time_step[a]) >= 1:
            time_steppdf[i, :], bin_edges = np.histogram(
                time_step[a], bins=bins, density=True
            )
            t[i, :] = (bin_edges[1:] + bin_edges[:-1]) / 2  # Record bin centers
        else:
            continue  # binseed not do anything, everything is zeroes.
    return time_steppdf, t
