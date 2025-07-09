"""
Functionality for determining the plasma potential of a Langmuir sweep.
"""

__all__ = ["find_plasma_potential_via_didv_peak"]

import numpy as np

from scipy import signal

from plasmapy.analysis.swept_langmuir.helpers import check_sweep, merge_voltage_clusters


def find_plasma_potential_via_didv_peak(
    voltage: np.ndarray,
    current: np.ndarray,
    # *,
    floating_potential: float | None = None,
    smooth_fractions=None,
):
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

    # define smooth_fractions
    # TODO: add better description
    if smooth_fractions is None:
        smooth_fractions = np.linspace(0.01, 0.25, num=25)

    # create bin sizes (smooth_windows) for the savgol_filter
    smooth_windows = np.unique(np.rint(smooth_fractions * data_size, dtype=int))
    mask = smooth_windows > 2
    smooth_windows = smooth_windows[mask]
    mask = smooth_windows % 2 == 0
    if np.count_nonzero(mask) > 0:
        smooth_windows[mask] = smooth_windows[mask] + 1
        smooth_windows = np.unique(smooth_windows)

    voltage_slice = voltage[first_index:]
    current_slice = current[first_index:]
    plasma_potentials = []
    for _window in smooth_windows:
        v_smooth = signal.savgol_filter(voltage_slice, _window, 1)
        c_smooth = signal.savgol_filter(current_slice, _window, 1)

        didv = np.gradient(c_smooth, v_smooth)
        imax = np.argmax(didv)
        if imax.size > 1:
            ...
        else:
            vp = voltage_slice[imax[0]]

        plasma_potentials.append(vp)




find_vp_ = find_plasma_potential_via_didv_peak
"""
Alias to
:func:`~plasmapy.analysis.swept_langmuir.plasma_potential.find_plasma_potential_via_didv_peak`.
"""
