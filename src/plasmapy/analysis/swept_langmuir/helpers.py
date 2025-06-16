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


def _force_regular_spacing(
    voltage: np.ndarray,
    current: np.ndarray,
    voltage_step_size: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Take a ``voltage`` array that already has voltage steps at some
    integer spacing of ``voltage_step_size`` and its paired ``current``
    array to generate a voltage array with a fixed step size
    ``voltage_step_size`` along with an appropriately `~numpy.nan`
    stuffed ``current`` array.

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

    voltage_step_size: `float`
        A non-zero, positive step size for the ``voltage`` array
        spacing.

    Returns
    -------
    voltage : `numpy.ndarray`
        An array with step size ``voltage_step_size`` spanning the same
        voltages given in the input ``voltage`` array.

    current : `numpy.ndarray`
        A `~numpy.nan` stuffed copy of the input ``current`` array to
        match the size of the returned ``voltage`` array.
    """
    size = int(np.round((voltage[-1] - voltage[0]) / voltage_step_size)) + 1

    if voltage.size == size:
        # all voltages steps are already equal
        return voltage, current

    reg_voltage = np.linspace(voltage[0], voltage[-1], num=size, dtype=voltage.dtype)

    _, reg_indices, new_indices = np.intersect1d(
        np.round(reg_voltage, decimals=5),
        np.round(voltage, decimals=5),
        return_indices=True,
    )

    mask = np.zeros_like(reg_voltage, dtype=bool)
    mask[reg_indices] = True

    reg_current = np.full(size, np.nan, dtype=current.dtype)
    reg_current[mask] = current[...]

    return reg_voltage, reg_current


def _is_voltage_regularly_spaced(
    voltage_diff: np.ndarray, mask_zero_diff: np.ndarray,
) -> bool:
    is_regular_grid = False
    if np.count_nonzero(mask_zero_diff) > 0:
        # is_regular_grid = False
        pass
    elif np.allclose(voltage_diff, voltage_diff[0]):
        # grid is already regularly spaced
        is_regular_grid = True
    elif np.allclose(
            np.rint(voltage_diff / np.min(voltage_diff))
            - voltage_diff / np.min(voltage_diff),
            0,
    ):
        # grid has a common dV, but at times jumps N * dV times
        is_regular_grid = True

    return is_regular_grid


def merge_voltage_clusters(  # noqa: C901, PLR0912, PLR0915
    voltage: np.ndarray,
    current: np.ndarray,
    voltage_step_size: float | None = None,
    force_regular_spacing: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Search the ``voltage`` array for closely spaced voltage clusters
    based on the ``voltage_step_size`` parameter and merge cluster, and
    associated ``current`` values, into a single point.

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

    voltage_step_size: `float` or `None`
        A non-zero, positive step size for the ``voltage`` array
        cluster identification.  A value of ``0`` will merge only
        duplicate voltage values.  A value of `None` will default to the
        average step size of the ``voltage`` array.  (DEFAULT:`None`)

    force_regular_spacing: `bool`
        If `False`, then the new voltage value of a cluster is the
        average of the contained points.  If `True`, then the new
        ``voltage`` array will be regularly spaced with a delta of
        ``voltage_step_size``.  (DEFAULT: `False`)

    Returns
    -------
    voltage : `numpy.ndarray`
        The new ``voltage`` array.

    current : `numpy.ndarray`
        The new ``current`` array.

    Notes
    -----
    ``voltage_step_size = 0`` :
        blah

    ``voltage_step_size > 0`` or `None`
        blah

    ``force_regular_spacing = True``
        blah
    """
    # condition force_regular_grid
    if not isinstance(force_regular_spacing, bool):
        raise TypeError(
            "Expected 'force_regular_spacing' to be a bool, but got type "
            f"{type(force_regular_spacing)}."
        )

    # condition voltage_step_size
    if voltage_step_size is not None and not isinstance(
        voltage_step_size, float | np.floating | int | np.integer
    ):
        raise TypeError(
            "Expected 'voltage_step_size' to be a float or None, got type "
            f"{type(voltage_step_size)}."
        )
    elif isinstance(voltage_step_size, int | np.integer):
        voltage_step_size = float(voltage_step_size)

    voltage, current = check_sweep(voltage, current)

    # check if grid is regularly spaced
    voltage_diff = np.diff(voltage)
    mask_zero_diff = np.isclose(voltage_diff, 0.0)
    is_regular_grid = _is_voltage_regularly_spaced(voltage_diff, mask_zero_diff)

    # return if voltage is already regularly spaced and no voltage merging is
    # requested
    if is_regular_grid and (voltage_step_size is None or voltage_step_size == 0):
        return voltage.copy(), current.copy()

    # condition voltage_step_size ... Round 2
    if voltage_step_size is None:
        voltage_step_size = np.abs(
            np.average(voltage_diff[np.logical_not(mask_zero_diff)])
        )
    elif voltage_step_size == 0:
        voltage_step_size = 0.0
    elif voltage_step_size < 0:
        voltage_step_size = -voltage_step_size

    # ensure voltage is ascending for calculation
    voltage_ascending = bool(np.all(voltage_diff >= 0))
    if not voltage_ascending:
        voltage, current = sort_sweep_arrays(
            voltage, current, voltage_order="ascending"
        )

    # now merge clusters
    voltage_diff = np.diff(voltage)
    if voltage_step_size != 0 and np.all(voltage_diff >= voltage_step_size):
        new_voltage = voltage.copy()
        new_current = current.copy()

        if force_regular_spacing:
            size = (
                int(np.round((new_voltage[-1] - new_voltage[0]) / voltage_step_size))
                + 1
            )
            reg_voltage = np.linspace(
                new_voltage[0], new_voltage[-1], num=size, dtype=new_voltage.dtype
            )
            reg_current = np.interp(reg_voltage, new_voltage, new_current)

            new_voltage = reg_voltage
            new_current = reg_current
    elif voltage_step_size == 0:
        voltage_diff = np.diff(voltage)
        mask_zero_diff = np.isclose(voltage_diff, 0.0)
        mask_zero_diff = np.append(mask_zero_diff, [mask_zero_diff[-1]])

        new_voltage = np.full(voltage.shape, np.nan, dtype=voltage.dtype)
        new_current = np.full(current.shape, np.nan, dtype=current.dtype)

        mask = np.logical_not(mask_zero_diff)
        new_voltage[mask] = voltage[mask]
        new_current[mask] = current[mask]

        while np.any(mask_zero_diff):
            volt = voltage[mask_zero_diff][0]
            index = np.where(new_voltage == volt)[0]
            if len(index) == 0:
                index = np.where(mask_zero_diff)[0]

            index = index[0]
            volt_mask = np.isclose(voltage, volt)

            new_voltage[index] = volt
            new_current[index] = np.average(current[volt_mask])

            mask_zero_diff[volt_mask] = False

        # remove nan entries
        mask = np.logical_not(np.isnan(new_voltage))
        new_voltage = new_voltage[mask]
        new_current = new_current[mask]

        if force_regular_spacing:
            new_voltage, new_current = merge_voltage_clusters(
                voltage=new_voltage,
                current=new_current,
                voltage_step_size=voltage_step_size,
            )

    else:
        new_voltage = np.full(voltage.shape, np.nan, dtype=voltage.dtype)
        new_current = np.full(current.shape, np.nan, dtype=current.dtype)

        start_voltage = voltage[0]
        stop_voltage = start_voltage + voltage_step_size

        # populate
        ii = 0
        while start_voltage <= voltage[-1]:
            mask1 = voltage >= start_voltage
            mask2 = voltage < stop_voltage
            mask = np.logical_and(mask1, mask2)

            if np.count_nonzero(mask) > 0:
                new_voltage[ii] = (
                    start_voltage + 0.5 * voltage_step_size
                    if force_regular_spacing
                    else np.average(voltage[mask])
                )
                new_current[ii] = np.average(current[mask])

                ii += 1

            start_voltage = stop_voltage
            stop_voltage += voltage_step_size

        # crop new arrays
        new_voltage = new_voltage[:ii]
        new_current = new_current[:ii]

        if force_regular_spacing:
            # Need to fill array with NaN values to force the regular spacing
            new_voltage, new_current = _force_regular_spacing(
                voltage=new_voltage,
                current=new_current,
                voltage_step_size=voltage_step_size,
            )

    if not voltage_ascending:
        voltage = voltage[::-1]
        current = current[::-1]
        new_voltage = new_voltage[::-1]
        new_current = new_current[::-1]

    return new_voltage, new_current
