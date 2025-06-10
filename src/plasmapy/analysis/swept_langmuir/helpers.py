"""Helper functions for analyzing swept Langmuir traces."""

__all__ = ["check_sweep", "merge_voltage_clusters", "sort_sweep_arrays"]

from typing import Literal

import astropy.units as u
import numpy as np


def check_sweep(  # noqa: C901, PLR0912
    voltage: np.ndarray,
    current: np.ndarray,
    strip_units: bool = True,
    allow_unsorted: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Function for checking that the voltage and current arrays are
    properly formatted for analysis by
    `plasmapy.analysis.swept_langmuir`.

    Parameters
    ----------
    voltage: `numpy.ndarray`
        1D `numpy.ndarray` representing the voltage of the swept
        Langmuir trace. Voltage should be monotonically increasing.
        *No units are assumed or checked, but values should be in
        volts.*

    current: `numpy.ndarray`
        1D `numpy.ndarray` representing the current of the swept
        Langmuir trace.  Values should start from a negative
        ion-saturation current and increase to a positive
        electron-saturation current.  *No units are assumed or checked,
        but values should be in amperes.*

    strip_units: `bool`
        (Default: `True`) If `True`, then the units on ``voltage``
        and/or ``current`` will be stripped if either are passed in as
        an Astropy `~astropy.units.Quantity`.

    allow_unsorted: `bool`
        (Default: `False`) If `True`, then the supplied ``voltage``
        array must be monotonically increasing.

    Returns
    -------
    voltage : `numpy.ndarray`
        Input argument ``voltage`` after it goes through all of its checks
        and conditioning.

    current : `numpy.ndarray`
        Input argument ``current`` after it goes through all of its checks
        and conditioning.

    Raises
    ------
    `TypeError`
        If either the ``voltage`` or ``current`` arrays are not instances of a
        `numpy.ndarray`.

    `ValueError`:
        If either the ``voltage`` or ``current`` arrays are not 1D.

    `ValueError`
        If the ``voltage`` array is not monotonically increasing.

    `ValueError`
        If the ``current`` array never crosses zero (i.e. has no floating
        potential).

    `ValueError`
        If the ``current`` array does not start form a negative ion-saturation
        current and increases to a positive electron-saturation current.

    `ValueError`
        If either the ``voltage`` or ``current`` array does not have a
        `numpy.dtype` of either `numpy.integer` or `numpy.floating`.

    """
    # -- examine voltage array --
    # check type
    if isinstance(voltage, np.ndarray):
        pass
    elif isinstance(voltage, list | tuple):
        voltage = np.array(voltage)
    else:
        raise TypeError(
            f"Expected 1D numpy array for voltage, but got {type(voltage)}.",
        )

    # check array structure
    if not (
        np.issubdtype(voltage.dtype, np.floating)
        or np.issubdtype(voltage.dtype, np.integer)
    ):
        raise ValueError(
            f"Expected 1D numpy array of floats or integers for voltage, but"
            f" got an array with dtype '{voltage.dtype}'."
        )
    elif voltage.ndim != 1:
        raise ValueError(
            f"Expected 1D numpy array for voltage, but got array with "
            f"{voltage.ndim} dimensions.",
        )
    elif not np.all(np.diff(voltage) >= 0) and not allow_unsorted:
        raise ValueError("The voltage array is not monotonically increasing.")

    if isinstance(voltage, u.Quantity) and strip_units:
        voltage = voltage.value

    # -- examine current array --
    # check type
    if isinstance(current, np.ndarray):
        pass
    elif isinstance(current, list | tuple):
        current = np.array(current)
    else:
        raise TypeError(
            f"Expected 1D numpy array for current, but got {type(current)}.",
        )

    # check array structure
    if not (
        np.issubdtype(current.dtype, np.floating)
        or np.issubdtype(current.dtype, np.integer)
    ):
        raise ValueError(
            f"Expected 1D numpy array of floats or integers for current, but"
            f" got an array with dtype '{current.dtype}'."
        )
    elif current.ndim != 1:
        raise ValueError(
            f"Expected 1D numpy array for current, but got array with "
            f"{current.ndim} dimensions.",
        )
    elif current.min() > 0.0 or current.max() < 0:
        raise ValueError(
            "Invalid swept Langmuir trace, the current never crosses zero "
            "'current = 0'."
        )
    elif current[0] > 0.0 or current[-1] < 0.0:
        raise ValueError(
            "The current array needs to start from a negative ion-saturation "
            "current to a positive electron-saturation current."
        )

    if voltage.size != current.size:
        raise ValueError(
            f"Incompatible arrays, 'voltage' size {voltage.size} must be the same"
            f" as the 'current' size {current.size}."
        )

    if isinstance(current, u.Quantity) and strip_units:
        current = current.value

    return voltage, current


def sort_sweep_arrays(
    voltage: np.ndarray,
    current: np.ndarray,
    voltage_order: Literal["ascending", "descending"] = "ascending",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sort the swept langmuir ``voltage`` and ``current`` traces to
    ensure the ``voltage`` array is either monotonically increasing or
    decreasing.

    Parameters
    ----------
    voltage: `numpy.ndarray`
        1D `numpy.ndarray` representing the voltage of the swept
        Langmuir trace.  *No units are assumed or checked, but values
        should be in volts.*

    current: `numpy.ndarray`
        1D `numpy.ndarray` representing the current of the swept
        Langmuir trace.  Values should start from a negative
        ion-saturation current and increase to a positive
        electron-saturation current.  *No units are assumed or checked,
        but values should be in amperes.*

    voltage_order: `str`
        Either 'ascending' or 'descending' to indicate how the
        ``voltage`` array should be sorted.  (DEFAULT: ``'ascending'``)

    Returns
    -------
    voltage : `numpy.ndarray`
        Sorted ``voltage`` array.

    current : `numpy.ndarray`
        Matched ``current`` array to the sorted ``voltage`` array.
    """
    if not isinstance(voltage_order, str):
        raise TypeError(
            "Expected 'voltage_order' to be a string equal to 'ascending' "
            f"or 'descending', but got type {type(voltage_order)}."
        )
    elif voltage_order not in ["ascending", "descending"]:
        raise ValueError(
            "Expected 'voltage_order' to be a string equal to 'ascending' "
            f"or 'descending', but got '{voltage_order}'."
        )

    voltage, current = check_sweep(
        voltage, current, strip_units=True, allow_unsorted=True
    )

    # determine order
    voltage_diff = np.diff(voltage)
    if np.all(voltage_diff >= 0):
        _order = "ascending"
    elif np.all(voltage_diff <= 0):
        _order = "descending"
    else:
        _order = None

    if _order == voltage_order:
        # already ordered
        return voltage, current

    # perform sorting
    if _order is None:
        index_sort = np.argsort(voltage)
        if voltage_order == "descending":
            index_sort = index_sort[::-1]

        voltage = voltage[index_sort]
        current = current[index_sort]
    else:
        voltage = voltage[::-1]
        current = current[::-1]

    return voltage, current


def merge_voltage_clusters(
    voltage: np.ndarray,
    current: np.ndarray,
    voltage_step_size: float = None,
    force_regular_grid: bool = False,
):
    # step_size = 0 (only merge identical voltages)
    # None : calculate steps size as average of array step sizes

    voltage, current = check_sweep(voltage, current)

    if not isinstance(force_regular_grid, bool):
        raise TypeError(
            "Expected 'force_regular_grid' to be a bool, but got type "
            f"{type(force_regular_grid)}."
        )

    if voltage_step_size is None:
        voltage_step_size = np.abs(np.average(np.diff(voltage)))
    elif voltage_step_size == 0:
        voltage_step_size = 0.0
    elif not isinstance(voltage_step_size, float):
        raise TypeError(
            "Expected 'voltage_step_size' to be a float or None, got type "
            f"{type(voltage_step_size)}."
        )
    elif voltage_step_size < 0:
        voltage_step_size = -voltage_step_size

    voltage_ascending = bool(np.all(np.diff(voltage) >= 0))
    if not voltage_ascending:
        voltage, current = sort_sweep_arrays(voltage, current, voltage_order="ascending")

    # now merge clusters
    new_voltage = np.unique(voltage)
    if voltage_step_size == 0 and new_voltage.size == voltage.size:
        new_current = current.copy()
    elif voltage_step_size == 0:
        new_current = np.empty_like(new_voltage, dtype=current.dtype)

        for ii in range(new_voltage.size):
            mask = voltage == new_voltage[ii]
            if np.count_nonzero(mask) == 1:
                new_current[ii] = current[mask]
            else:
                new_current[ii] = np.average(current[mask])
    else:
        new_voltage = np.empty_like(voltage, dtype=voltage.dtype)
        new_voltage[...] = np.nan
        new_current = np.empty_like(current, dtype=current.dtype)
        new_current[...] = np.nan

        start_voltage = voltage[0]
        stop_voltage = start_voltage + voltage_step_size

        ii = 0
        while start_voltage <= voltage[-1]:
            mask1 = voltage >= start_voltage
            mask2 = voltage < stop_voltage
            mask = np.logical_and(mask1, mask2)

            if np.count_nonzero(mask) > 0:
                new_voltage[ii] = (
                    start_voltage + 0.5 * voltage_step_size if force_regular_grid
                    else np.average(voltage[mask])
                )
                new_current[ii] = np.average(current[mask])

                ii += 1

            start_voltage = stop_voltage
            stop_voltage = start_voltage + voltage_step_size

        # filter out NaN values
        new_voltage = new_voltage[:ii]
        new_current = new_current[:ii]

        # force regular spacing
        if force_regular_grid:
            # Need to fill array with NaN values to force the regular spacing
            ...

    if not voltage_ascending:
        voltage, current = sort_sweep_arrays(voltage, current, voltage_order="descending")
        new_voltage, new_current = sort_sweep_arrays(
            new_voltage, new_current, voltage_order="descending"
        )
    else:
        new_voltage, new_current = sort_sweep_arrays(
            new_voltage, new_current, voltage_order="ascending"
        )

    return new_voltage, new_current
