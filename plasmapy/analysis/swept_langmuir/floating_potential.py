"""Functionality for determining the floating potential of a Langmuir sweep."""
__all__ = ["find_floating_potential"]

import numpy as np

from collections import namedtuple
from warnings import warn
from typing import Union

from plasmapy.analysis.swept_langmuir.fit_functions import (
    ExponentialOffsetFitFunction,
    LinearFitFunction,
)

FloatingPotentialResults = namedtuple("FloatingPotentialResults",
                                      ("vf", "vf_err", "info"))


def find_floating_potential(
        voltage: np.ndarray,
        current: np.ndarray,
        threshold: int = 1,
        min_points: Union[int, float] = None,
        fit_type: str = "exponential",
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
    fit_funcs = {
        "linear": {
            "func": LinearFitFunction(),
            "min_point_factor": 0.1,
        },
        "exponential": {
            "func": ExponentialOffsetFitFunction(),
            "min_point_factor": 0.2,
        },
    }
    try:
        fit_func = fit_funcs[fit_type]["func"]
        min_point_factor = fit_funcs[fit_type]["min_point_factor"]
        meta_dict = {"func": fit_func}
    except KeyError:
        raise KeyError(
            f"Requested fit function '{fit_type}' is  not a valid option.  "
            f"Examine kwarg 'fit_curve' for valid options.")

    if current.min() > 0.0 or current.max() < 0:
        warn("The Langmuir sweep has no floating potential.")

        return FloatingPotentialResults(np.nan, np.nan, meta_dict)

    # check voltage is monotonically increasing/decreasing
    voltage_diff = np.diff(voltage)
    if not (np.all(voltage_diff >= 0) or np.all(voltage_diff <= 0)):
        warn("The voltage array is not monotonically increasing or decreasing.")

        return FloatingPotentialResults(np.nan, np.nan, meta_dict)

    # condition kwarg threshold
    if isinstance(threshold, (int, float)):
        threshold = int(np.round(threshold))
        if threshold < 1:
            warn(f"threshold ({threshold}) is less than 1 and needs to be"
                 f" an int >= 1, using a value of 1")
            threshold = 1
    else:
        warn(f"threshold is NOT a integer >= 1, using a value of 1")

    # condition min_points
    if min_points is None:
        min_points = int(np.max([5, np.around(min_point_factor * voltage.size)]))
    elif min_points == 0:
        # this signals to use all points
        pass
    elif 0 < min_points < 1:
        min_points = int(np.round(min_points * voltage.size))
    elif min_points >= 1:
        min_points = int(np.round(min_points))
    else:
        raise ValueError(f"Got {min_points}, but 'min_points' must be an int or float "
                         f"greater than or equal to 0.")

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

    if min_points == 0:
        meta_dict['islands'] = [slice(cp_candidates[0], cp_candidates[-1]+1)]
    elif n_islands == 1:
        meta_dict['islands'] = [slice(cp_candidates[0], cp_candidates[-1]+1)]
    else:
        # There are multiple crossing points
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

        # do islands fall within min_points window
        isl_window = np.abs(np.r_[meta_dict["islands"][-1]][-1]
                            - np.r_[meta_dict["islands"][0]][0]) + 1
        if isl_window > min_points:
            warn(f"Unable to determine floating potential, Langmuir sweep has "
                 f"{n_islands} crossing-islands.  Try adjusting keyword 'threshold' "
                 f"and/or smooth the current.")

            return FloatingPotentialResults(np.nan, np.nan, meta_dict)

    # Construct crossing-island (pad if needed)
    if min_points == 0:
        # us all points
        istart = 0
        istop = voltage.size
    else:
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
    fit_func.curve_fit(volt_sub, curr_sub)

    vf, vf_err = fit_func.root_solve()
    meta_dict.update({
        "rsq": fit_func.rsq,
        "indices": slice(istart, istop+1)
    })

    return FloatingPotentialResults(vf, vf_err, meta_dict)
