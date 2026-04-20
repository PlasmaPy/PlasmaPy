"""
Functionality for determining the plasma potential of a Langmuir sweep.
"""

__all__ = ["find_didv_peak_location", "dIdVExtras"]
__aliases__ = []
__all__ += __aliases__

from collections.abc import Sequence
from typing import NamedTuple

import numpy as np
from scipy import signal

from plasmapy.analysis.swept_langmuir.helpers import (
    _condition_voltage_window,
    check_sweep,
    merge_voltage_clusters,
)


class dIdVExtras(NamedTuple):  # noqa: N801
    """
    `~typing.NamedTuple` structured to contain the extra parameters
    calculated by
    `~plasmapy.analysis.swept_langmuir.plasma_potential.find_didv_peak_location`.
    """

    std: float | None
    """
    Alias for field number 0, standard deviation of all the computed
    peak slope bias locations for each Savitzky-Golay filtered Langmuir
    trace.  Standard deviation of :attr:`savgol_peaks`.
    """

    data_slice: slice | None
    """
    Alias for field number 1, `slice` objected corresponding to the
    sub-arrays of the Langmuir trace used for the calculation.
    """

    savgol_windows: list[int] | None
    """
    Alias for field number 2, list of windows sizes used for each
    Savitzky-Golay filtered Langmuir trace.
    """

    savgol_peaks: list[float] | None
    """
    Alias for field number 3, list of computed peak slope bias locations
    for each Savitzky-Golay filtered Langmuir trace.
    """


def _condition_smooth_fractions(  # noqa: C901
    smooth_fractions, data_size: int
):
    """
    Condition ``smooth_fractions`` and return the resulting
    Savitzky-Golay filter windows sizes, based on the ``data_size``.
    """
    # condition smooth_fractions
    if smooth_fractions is None:
        # Note: If this default value is changed, then the docstring entry
        #       for smooth_fractions in the find_didv_peak() docstring
        #       needs to be updated accordingly.
        smooth_fractions = np.linspace(0.01, 0.25, num=30)
    elif isinstance(smooth_fractions, np.ndarray):
        pass
    elif not isinstance(smooth_fractions, Sequence) or isinstance(
        smooth_fractions, str
    ):
        raise TypeError(
            f"Expected a list-like object for 'smooth_fractions', "
            f"but got type {type(smooth_fractions)}."
        )
    else:
        smooth_fractions = np.array(smooth_fractions)

    if smooth_fractions.ndim != 1:
        raise ValueError(
            "Expected a 1-D list of floats in the interval [0, 1] for argument "
            f"'smooth_fractions', but got a {smooth_fractions.ndim}-D list."
        )
    elif not np.issubdtype(smooth_fractions.dtype, np.floating):
        raise ValueError(
            "Expected a 1-D list of floats in the interval [0, 1] for argument "
            "'smooth_fractions', not all elements are floats."
        )

    smooth_fractions = np.unique(np.sort(smooth_fractions))
    mask1 = smooth_fractions < 0
    mask2 = smooth_fractions > 1
    mask = np.logical_or(mask1, mask2)
    if np.count_nonzero(mask) > 0:
        raise ValueError(
            "Expected a 1-D list of floats in the interval [0, 1] for argument "
            f"'smooth_fractions', no elements are within this interval "
            f"{smooth_fractions.tolist()}."
        )

    # condition data_size
    if not isinstance(data_size, int):
        raise TypeError(
            "Expected a positive, non-zero integer for argument 'data_size', "
            f"but got type {type(data_size)}."
        )
    elif data_size < 1:
        raise ValueError(
            "Expected a positive, non-zero integer for argument 'data_size', "
            f"but got {data_size}."
        )

    # create bin sizes (savgol_windows) for the savgol_filter
    savgol_windows = np.unique(np.rint(smooth_fractions * data_size).astype(int))

    # windows need to have at least 2 points
    mask = savgol_windows > 2
    savgol_windows = savgol_windows[mask]

    if 0.0 in smooth_fractions:
        savgol_windows = np.sort(np.append(savgol_windows, 0))

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


def find_didv_peak_location(
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
        An order list of fractions in the interval :math:`[0, 1]` used
        to compute the Savitzky-Golay filter window sizes.  For example,
        if the ``voltage_windows`` had a size of 50, then a
        ``smooth_fraction`` value of 0.5 would result in a
        Savitzky-Golay window size of 25.  For ``smooth_fraction``
        value of ``0.0`` no smoother will be performed.  If `None`
        (default), then ``smooth_fractions`` will default to
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
        if _window == 0:
            v_smooth = voltage_slice
            c_smooth = current_slice
        else:
            v_smooth = signal.savgol_filter(voltage_slice, _window, 1)
            c_smooth = signal.savgol_filter(current_slice, _window, 1)

        didv = np.gradient(c_smooth, v_smooth)
        didv_max = np.max(didv)
        imax = np.where(np.isclose(didv, didv_max, atol=0.0))[0]
        if imax.size > 1:
            if np.all(np.diff(imax) == 1):
                # neighboring indices have the same peak value, take
                # voltage average
                vp = np.average(voltage_slice[imax])
            else:
                # peaks occur at separated indices (i.e. NOT neighboring
                # indices) ... unable to determine vp just go to the next
                # savgol window
                continue
        elif np.isscalar(imax):
            vp = voltage_slice[imax]
        else:
            vp = voltage_slice[imax[0]]

        plasma_potentials.append(float(vp))
        rtn_extras["savgol_windows"].append(int(_window))

    if len(plasma_potentials) == 0:
        raise RuntimeError("Unable to identify any dI/dV peaks in the Langmuir trace.")

    rtn_extras["savgol_peaks"] = plasma_potentials
    rtn_extras["std"] = float(np.std(plasma_potentials))

    vp = float(np.average(plasma_potentials))
    return vp, dIdVExtras(**rtn_extras)
