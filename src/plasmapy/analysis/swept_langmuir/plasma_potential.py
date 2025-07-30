"""
Functionality for determining the plasma potential of a Langmuir sweep.
"""

__all__ = ["find_didv_peak"]
__aliases__ = []
__all__ += __aliases__

from typing import NamedTuple

import numpy as np

from scipy import signal

from plasmapy.analysis.swept_langmuir.helpers import check_sweep, merge_voltage_clusters


class dIdVExtras(NamedTuple):
    std: float | None
    data_slice: slice | None
    savgol_windows: list[int] | None
    savgol_peaks: list[float] | None


def find_didv_peak(
    voltage: np.ndarray,
    current: np.ndarray,
    # *,
    floating_potential: float | None = None,
    smooth_fractions=None,
):
    rtn_extras = dIdVExtras(
        std=None,
        data_slice=None,
        savgol_windows=None,
        savgol_peaks=None,
    )._asdict()

    # check voltage and current arrays
    voltage, current = check_sweep(voltage, current, strip_units=True)
    voltage, current = merge_voltage_clusters(voltage, current, voltage_step_size=0)

    # condition floating potential
    if floating_potential is None:
        pass
    elif not isinstance(floating_potential, float | int | np.floating | np.integer):
        raise TypeError(
            f"Expected None or float for floating_potential, got type "
            f"{type(floating_potential)}."
        )
    elif not (np.min(voltage) <= floating_potential <= np.max(voltage)):
        raise ValueError(
            f"Given floating_potential ({floating_potential}) is outside range "
            f"[{np.min(voltage)}, {np.max(voltage)}]."
        )

    # determine data window
    first_index = 0
    if floating_potential is not None:
        mask = voltage > floating_potential
        if np.count_nonzero(mask) > 0:
            first_index = np.where(mask)[0][0]
    data_size = voltage.size - first_index
    rtn_extras["data_slice"] = slice(int(first_index), None, 1)

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

    voltage_slice = voltage[first_index:]
    current_slice = current[first_index:]
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
