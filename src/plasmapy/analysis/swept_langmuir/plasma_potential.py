"""
Functionality for determining the plasma potential of a Langmuir sweep.
"""

__all__ = ["find_didv_peak_location", "dIdVExtras"]
__aliases__ = []
__all__ += __aliases__

import numbers
from collections.abc import Sequence
from typing import NamedTuple

import numpy as np
from scipy import signal

from plasmapy.analysis.swept_langmuir.helpers import check_sweep, merge_voltage_clusters


class dIdVExtras(NamedTuple):  # noqa: N801
    std: float | None
    data_slice: slice | None
    savgol_windows: list[int] | None
    savgol_peaks: list[float] | None


def _condition_voltage_window(voltage, voltage_window) -> slice:
    """
    Condition ``voltage_window`` and return resulting `slice` object to
    index ``voltage``.
    """
    if isinstance(voltage_window, np.ndarray):
        voltage_window = voltage_window.tolist()

    if voltage_window is None:
        voltage_window = [None, None]
    elif not isinstance(voltage_window, Sequence):
        raise TypeError(
            f"Expected a 2-element list of floats or None for 'voltage_window', "
            f"but got type {type(voltage_window)}."
        )
    elif len(voltage_window) != 2:
        raise ValueError(
            f"Expected a 2-element list of floats or None for 'voltage_window', "
            f"but got type {len(voltage_window)} elements."
        )
    elif not all(
        isinstance(element, numbers.Real) or element is None
        for element in voltage_window
    ):
        raise TypeError(f"Not all elements of 'voltage_window' are floats or None.")
    elif None not in voltage_window:
        voltage_window = np.sort(voltage_window).tolist()

    # determine data window
    if voltage_window[0] is not None and voltage_window[0] >= voltage[-1]:
        raise ValueError(
            f"The min value for the voltage window ({voltage_window[0]}) "
            f"is larger than the max value of the langmuir trace "
            f"({voltage[-1]})."
        )

    if voltage_window[1] is not None and voltage_window[1] <= voltage[0]:
        raise ValueError(
            f"The max value for the voltage window ({voltage_window[1]}) "
            f"is smaller than the min value of the langmuir trace "
            f"({voltage[0]})."
        )

    if all(isinstance(element, numbers.Real) for element in voltage_window):
        voltage_window = np.sort(voltage_window).tolist()

    first_index = (
        None
        if voltage_window[0] is None
        else int(np.where(voltage >= voltage_window[0])[0][0])
    )

    last_index = (
        None
        if voltage_window[1] is None
        else int(np.where(voltage <= voltage_window[1])[0][-1])
    )

    return slice(first_index, last_index, 1)


def _condition_smooth_fractions(smooth_fractions, data_size):
    """
    Condition ``smooth_fractions`` and return the resulting
    Savitzky-Golay filter windows sizes, based on the ``data_size``.
    """
    if smooth_fractions is None:
        # Note: If this default value is changed, then the docstring entry
        #       for smooth_fractions in the find_didv_peak() docstring
        #       needs to be updated accordingly.
        smooth_fractions = np.linspace(0.01, 0.25, num=30)
    elif (
        isinstance(smooth_fractions, Sequence)
        and not isinstance(smooth_fractions, np.ndarray)
    ):
        smooth_fractions = np.array(smooth_fractions)

    if not isinstance(smooth_fractions, np.ndarray):
        raise TypeError(
            "Expected a 1-D list of floats in the interval (0, 1] for argument "
            f"'smooth_fractions', but got type {type(smooth_fractions)}."
        )
    elif smooth_fractions.ndim != 1:
        raise ValueError(
            "Expected a 1-D list of floats in the interval (0, 1] for argument "
            f"'smooth_fractions', but got a {smooth_fractions.ndim}-D list.")
    elif not np.issubdtype(smooth_fractions.dtype, np.floating):
        raise ValueError(
            "Expected a 1-D list of floats in the interval (0, 1] for argument "
            f"'smooth_fractions', not all elements are floats."
        )

    smooth_fractions = np.unique(np.sort(smooth_fractions))
    mask1 = smooth_fractions > 0
    mask2 = smooth_fractions <= 1
    mask = np.logical_and(mask1, mask2)
    if np.count_nonzero(mask) == 0:
        raise ValueError(
            "Expected a 1-D list of floats in the interval (0, 1] for argument "
            f"'smooth_fractions', no elements are within this interval "
            f"{smooth_fractions.tolist()}."
        )

    # create bin sizes (savgol_windows) for the savgol_filter
    savgol_windows = np.unique(np.rint(smooth_fractions * data_size).astype(int))

    # windows need to have at least 2 points
    mask = savgol_windows > 2
    savgol_windows = savgol_windows[mask]

    # force windows sizes to be odd
    mask = savgol_windows % 2 == 0
    if np.count_nonzero(mask) > 0:
        savgol_windows[mask] = savgol_windows[mask] + 1
        savgol_windows = np.unique(savgol_windows)

    # do not let windows sizes exceed data_size
    mask = savgol_windows <= data_size
    savgol_windows = savgol_windows[mask]

    # check savgol_windows is not null
    if savgol_windows.size == 0:
        raise ValueError(
            f"The given smooth_fractions ({smooth_fractions}) and "
            f"window size ({data_size}) resulted in no valid Savitzky-Golay "
            f"filter windows.  Computed windows must be odd, greater than 3, "
            f"and less than or equal to the windows size."
        )

    return savgol_windows


def find_didv_peak_location(  # noqa: C901, PLR0912
    voltage: np.ndarray,
    current: np.ndarray,
    *,
    voltage_window: list[float | None] | None = None,
    smooth_fractions: list[float] | None = None,
) -> tuple[float, dIdVExtras]:
    """
    Find the bias voltage at which the peak slope (:math:`dI/dV_{max}`)
    of the swept Langmuir trace occurs.

    The peak slope bias voltage is often used as a rough estimate of
    the plasma potential.  However, the value will always be slightly
    less than the actual plasma potential.

    Parameters
    ----------
    voltage : `numpy.ndarray`
        1-D numpy array of monotonically increasing probe biases
        (should be in volts).

    current : `numpy.ndarray`
        1-D numpy array of probe current (should be in amperes)
        corresponding to the ``voltage`` array.

    voltage_window : `list[float | None]` | `None`, default: `None`
        A two-element list ``[v_min, v_max]`` that specifies the voltage
        range in which the peak slope will be looked for.  Specifying
        `None` for either the first or second element will result in a
        window using the lower or upper bound of ``voltage``
        respectively.  If set to `None` (default), then the whole
        ``voltage`` window will be used.

    smooth_fractions : `list[float]` | `None`, default: `None`
        An order list of fractions in the interval :math:`(0, 1]` used
        to compute the Savitzky-Golay filter window sizes.  For example,
        if the ``voltage_windows`` had a size of 50, then a
        ``smooth_fraction`` value of 0.5 would result in a
        Savitzky-Golay window size of 25.  If `None` (default), then
        ``smooth_fractions`` will default to
        ``numpy.linspace(0.01, 0.25, num=30)``.

    Notes
    -----

    Add details about algorithm.
    """
    rtn_extras = dIdVExtras(
        std=None,
        data_slice=None,
        savgol_windows=None,
        savgol_peaks=None,
    )._asdict()

    # check voltage and current arrays
    voltage, current = check_sweep(voltage, current, strip_units=True)
    voltage, current = merge_voltage_clusters(voltage, current, voltage_step_size=0)

    # condition voltage_window
    _slice = _condition_voltage_window(voltage, voltage_window)
    rtn_extras["data_slice"] = _slice
    data_size = voltage[_slice].size
    if data_size < 3:
        raise ValueError(
            f"The specified voltage_window ({voltage_window}) would result "
            f"in a null window or a window with fewer than 3 elements."
        )

    # define starting savgol windows
    savgol_windows = _condition_smooth_fractions(smooth_fractions, data_size)

    voltage_slice = voltage[_slice]
    current_slice = current[_slice]
    plasma_potentials = []
    rtn_extras["savgol_windows"] = []
    for _window in savgol_windows:
        v_smooth = signal.savgol_filter(voltage_slice, _window, 1)
        c_smooth = signal.savgol_filter(current_slice, _window, 1)

        didv = np.gradient(c_smooth, v_smooth)
        imax = np.argmax(didv)
        if imax.size > 1:
            if np.all(np.diff(imax) == 1):
                vp = np.average(voltage_slice[imax])
            else:
                continue
        elif np.isscalar(imax):
            vp = voltage_slice[imax]
        else:
            vp = voltage_slice[imax[0]]

        plasma_potentials.append(float(vp))
        rtn_extras["savgol_windows"].append(int(_window))

    rtn_extras["savgol_peaks"] = plasma_potentials
    rtn_extras["std"] = float(np.std(plasma_potentials))

    vp = float(np.average(plasma_potentials))
    return vp, dIdVExtras(**rtn_extras)
