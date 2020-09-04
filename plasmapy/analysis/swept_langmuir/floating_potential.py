__all__ = ["find_floating_potential"]

import numpy as np

from scipy.stats import linregress
from warnings import warn


def find_floating_potential(
        voltage: np.ndarray,
        current: np.ndarray,
        threshold: int = 1,
        min_points: int = 5,
        retfit: bool = False,
        retislands: bool = False,
):
    """
    Determines the floating potential (Vf) for a given Current-Voltage (IV) curve
    generate by a swept Langmuir probe.  The floating potential is the probe
    bias where the collected current goes to zero.

    How it works
    ------------
    #. The current array ``current` is searched for all points equal to zero and
       point pairs that straddle  ``current = 0`` to form a set of crossing-points.
    #. Crossing-points are then grouped into crossing-islands in accordance to
       the ``threshold`` keyword.
    #. If multiple crossing-islands are found, then an warning issued and
       `~numpy.nan` is returned.
    #. To calculated the floating potential, a `~scipy.stats.linregress` is applied
       to points making up the cross-island.  If the number of points that make
       up the crossing-island are less than ``min_points``, then each side of the
       crossing-island is padded with the nearest neighbors until `min_points` is
       satisfied.

    Parameters
    ----------

    voltage: np.ndarray
        1-D numpy array of ascending probe biases (in Volts)

    current: np.ndarray
        1-D numpy array of probe current (in A) corresponding to the :data:`voltage`
        array

    threshold: positive, non-zero int
        Max index distance between found crossing-points to group crossing-points
        into crossing-islands.  That is, if `threshold=5` then consecutive
        crossing-points are considered to be in the same crossing-island if they are
        within 5 indices of each other. (Default: 1)

    min_points: positive, non-zero int
        The minimum number of points required for the linear regression. (Default: 5)

    Returns
    -------
    vf: `numpy.float64` or `numpy.nan`
        The calculated floating potential (in Volts).  Returns `numpy.nan` if the
        floating potential can not be determined.

    vf_err: `numpy.float64` or `numpy.nan`
        The error associated with the floating potential calculation (in Volts).
        Returns `numpy.nan` if the floating potential can not be determined.

    fit: Dict[str, Any]
        A dictionary containing the linear regression fit results and parameters.
        Keys are `'slope'`, `'slope_err'`, `'intercept'`, `'intercept_err'`, and
        `'indices'`, where `'indices'` is a `slice` object corresponding to the
        data points used in the fit.  Returns an empty dict if the floating
        potential can not be determined.

    """
    meta_dict = {}

    if current.min() > 0.0 or current.max() < 0:
        warn("The Langmuir sweep has no floating potential.")

        ret = (np.nan, np.nan)
        if retislands or retfit:
            ret += (meta_dict,)

        return ret

    # check voltage is monotonically increasing/decreasing
    voltage_diff = np.diff(voltage)
    if not (np.all(voltage_diff >= 0) or np.all(voltage_diff <= 0)):
        warn("The voltage array is not monotonically increasing or decreasing.")

        ret = (np.nan, np.nan)
        if retislands or retfit:
            ret += (meta_dict,)

        return ret

    # condition kwarg threshold
    if isinstance(threshold, (int, float)):
        threshold = int(np.round(threshold))
        if threshold < 1:
            warn(f"threshold ({threshold}) is less than 1 and needs to be"
                 f" an int >= 1, using a value of 1")
            threshold = 1
    else:
        warn(f"threshold is NOT a integer >= 1, using a value of 1")

    # find possible crossing points (cp)
    lower_vals = np.where(current <= 0, True, False)
    upper_vals = np.where(0 <= current, True, False)
    cp_exact = np.logical_and(lower_vals, upper_vals).nonzero()[0]
    cp_low2high = np.logical_and(np.roll(lower_vals, 1), upper_vals).nonzero()[0]
    cp_high2low = np.logical_and(np.roll(lower_vals, -1), upper_vals).nonzero()[0]

    # adjust for array wrapping cause by np.roll
    cp_low2high = cp_low2high[np.where(cp_low2high != 0, True, False)]
    cp_high2low = cp_high2low[np.where(cp_high2low != current.size-1, True, False)]

    # collect all candidates
    cp_candidates = np.concatenate((
        cp_exact,
        cp_low2high,
        cp_low2high - 1,
        cp_high2low,
        cp_high2low + 1,
    ))
    cp_candidates = np.unique(cp_candidates)  # sorted and unique

    # How many crossing-islands?
    cp_intervals = np.diff(cp_candidates)
    threshold_indices = np.where(cp_intervals > threshold)[0]
    n_islands = threshold_indices.size + 1
    if retislands:
        if n_islands == 1:
            meta_dict['islands'] = [slice(cp_candidates[0], cp_candidates[-1]+1)]
        else:
            isl_start = np.concatenate((
                [cp_candidates[0]],
                cp_candidates[threshold_indices+1],
            ))
            isl_stop = np.concatenate((
                cp_candidates[threshold_indices]+1,
                [cp_candidates[-1]+1],
            ))
            meta_dict['islands'] = []
            for start, stop in zip(isl_start, isl_stop):
                meta_dict['islands'].append(slice(start, stop))
    if n_islands != 1:
        # There are multiple crossing points
        warn(f"Unable to determine floating potential, Langmuir sweep has "
             f"{n_islands} crossing-islands.  Try adjusting keyword 'threshold'"
             f"and/or smooth the current.")

        ret = (np.nan, np.nan)
        if retislands or retfit:
            ret += (meta_dict,)

        return ret

    # Construct crossing-island (pad if needed)
    istart = cp_candidates[0]
    istop = cp_candidates[-1]
    iadd = (istop - istart + 1) - min_points
    if iadd < 0:
        # pad front
        ipad_2_start = ipad_2_stop = int(np.ceil(-iadd / 2.0))
        if istart - ipad_2_start < 0:
            ipad_2_stop += ipad_2_start - istart
            ipad_2_start = 0
            istart = 0
        else:
            istart -= ipad_2_start
            ipad_2_start = 0

        # pad rear
        if ((current.size - 1) - (istop + ipad_2_stop)) < 0:
            ipad_2_start += ipad_2_stop - (current.size - 1 - istop)
            istop = current.size - 1
        else:
            istop += ipad_2_stop

        # re-pad front if possible
        if ipad_2_start > 0:
            if istart - ipad_2_start < 0:
                istart = 0
            else:
                istart -= ipad_2_start

    if (istop - istart + 1) < min_points:
        warn(f"The number of elements in the current array "
             f"({istop - istart + 1}) is less than 'min_points' "
             f"({min_points}).")

    # Perform Linear Regression Fit
    volt_sub = voltage[istart:istop + 1]
    curr_sub = current[istart:istop + 1]
    fit = linregress(volt_sub, curr_sub)

    slope = fit[0]
    slope_err = fit[4]

    intercept = fit[1]
    intercept_err = np.sum(volt_sub ** 2) - ((np.sum(volt_sub) ** 2) / volt_sub.size)
    intercept_err = slope_err * np.sqrt(1.0 / intercept_err)

    vf = -intercept / slope
    vf_err = (
        np.abs(vf)
        * np.sqrt(((slope_err / slope) ** 2)
                  + ((intercept_err / intercept) ** 2))
    )

    if retfit:
        meta_dict.update({
            'slope': slope,
            'slope_err': slope_err,
            'intercept': intercept,
            'intercept_err': intercept_err,
            'rsq': fit[2],
            'indices': slice(istart, istop + 1),
        })

    ret = (vf, vf_err)
    if retislands or retfit:
        ret += (meta_dict,)

    return ret
