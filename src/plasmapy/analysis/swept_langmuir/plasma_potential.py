"""
Functionality for determining the plasma potential of a Langmuir sweep.
"""

__all__ = ["find_didv_peak"]
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
        voltage_window = np.sort(voltage).tolist()

    # determine data window
    if (
        voltage_window[0] is None
        or voltage_window[0] <= voltage[0]
    ):
        first_index = None
    elif voltage_window[0] >= voltage[-1]:
        raise ValueError(
            f"The min value for the voltage window ({voltage_window[0]}) "
            f"is larger than the max value of the langmuir trace "
            f"({voltage[-1]})."
        )
    else:
        first_index = int(np.where(voltage >= voltage_window[0])[0][0])

    if (
        voltage_window[1] is None
        or voltage_window[1] >= voltage[-1]
    ):
        last_index = None
    elif voltage_window[1] <= voltage[0]:
        raise ValueError(
            f"The max value for the voltage window ({voltage_window[1]}) "
            f"is smaller than the min value of the langmuir trace "
            f"({voltage[0]})."
        )
    else:
        last_index = int(np.where(voltage < voltage_window[1])[0][0])

    return slice(first_index, last_index, 1)


def find_didv_peak(  # noqa: C901, PLR0912
    voltage: np.ndarray,
    current: np.ndarray,
    *,
    voltage_window: list[float | None] | None = None,
    smooth_fractions=None,
) -> tuple[float, dIdVExtras]:
    """
    Find the peak slope (:math:`dI/dV_{max}`) of the swept Langmuir
    trace.

    This is often used as a rough estimate of the plasma potential.
    However, it will always be slightly less than the actual plasma
    potential.

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
    data_size = len(_slice.indices(voltage.size))
    if data_size <= 1:
        raise ValueError(
            f"The specified voltage_window ({voltage_window}) would result "
            f"in a null window or a 1-element window."
        )

    # define smooth_fractions
    # TODO: add better description
    if smooth_fractions is None:
        smooth_fractions = np.linspace(0.01, 0.25, num=30)

    # create bin sizes (smooth_windows) for the savgol_filter
    smooth_windows = np.unique(np.rint(smooth_fractions * data_size).astype(int))
    mask = smooth_windows > 2
    smooth_windows = smooth_windows[mask]
    mask = smooth_windows % 2 == 0
    if np.count_nonzero(mask) > 0:
        smooth_windows[mask] = smooth_windows[mask] + 1
        smooth_windows = np.unique(smooth_windows)

    voltage_slice = voltage[_slice]
    current_slice = current[_slice]
    plasma_potentials = []
    rtn_extras["savgol_windows"] = []
    for _window in smooth_windows:
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
