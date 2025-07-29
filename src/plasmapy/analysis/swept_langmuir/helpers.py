"""Helper functions for analyzing swept Langmuir traces."""

__all__ = ["check_sweep", "merge_voltage_clusters", "sort_sweep_arrays"]

import numbers
import warnings
from typing import Literal

import astropy.units as u
import numpy as np

from plasmapy.utils.exceptions import PlasmaPyWarning


def check_sweep(  # noqa: C901, PLR0912
    voltage: np.ndarray,
    current: np.ndarray,
    strip_units: bool = True,
    allow_unsorted: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Check that the voltage and current arrays are properly formatted
    for analysis by `plasmapy.analysis.swept_langmuir`.

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

    strip_units: `bool`, default: `True`
        If `True`, then the units on ``voltage`` and/or ``current``
        will be stripped if either are passed in as an Astropy
        `~astropy.units.Quantity`.

    allow_unsorted: `bool`, default: `False`
        If `True`, then the supplied ``voltage`` array must be
        monotonically increasing.

    Returns
    -------
    voltage : `numpy.ndarray`
        Input argument ``voltage`` after it goes through all of its
        checks and conditioning.

    current : `numpy.ndarray`
        Input argument ``current`` after it goes through all of its
        checks and conditioning.

    Raises
    ------
    `TypeError`
        If either the ``voltage`` or ``current`` arrays are not
        instances of a `numpy.ndarray`.

    `ValueError`:
        If either the ``voltage`` or ``current`` arrays are not 1D.

    `ValueError`
        If the ``voltage`` array is not monotonically increasing.

    `ValueError`
        If the ``current`` array never crosses zero (i.e. has no
        floating potential).

    `ValueError`
        If the ``current`` array does not start form a negative
        ion-saturation current and increases to a positive
        electron-saturation current.

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
        Either ``'ascending'`` or ``'descending'`` to indicate how the
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


def _is_voltage_regularly_spaced(
    voltage_diff: np.ndarray,
    mask_zero_diff: np.ndarray,
) -> bool:
    """
    Determine if the voltage difference array ``voltage_diff`` is
    regularly spaced; that is the differences are all equal or some
    integer multiple of the smallest difference.
    """
    is_regular_grid = False
    if np.count_nonzero(mask_zero_diff) > 0:
        # is_regular_grid = False
        return False
    elif np.allclose(voltage_diff, voltage_diff[0]):
        # grid is already regularly spaced
        is_regular_grid = True
    elif np.allclose(
        np.rint(voltage_diff / np.min(voltage_diff))
        - voltage_diff / np.min(voltage_diff),
        0,
    ):
        # grid has a common dV, but at times jumps N * dV times
        min_dV = np.min(voltage_diff)
        ndV = np.rint(voltage_diff / min_dV)

        is_regular_grid = not np.any(ndV > 10)

    return is_regular_grid


def _merge_voltage_clusters__zero_diff_neighbors(
    voltage: np.ndarray, current: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Take the ``voltage`` and ``current`` arrays associated with a swept
    langmuir trace and average together clusters of identical voltage
    values.
    """
    # generate boolean mask for zero diff locations
    voltage_diff = np.diff(voltage)
    mask_zero_diff = np.isclose(voltage_diff, 0.0)
    mask_zero_diff = np.append(mask_zero_diff, [mask_zero_diff[-1]])

    # initialize new voltage and current arrays
    new_voltage = np.full(voltage.shape, np.nan, dtype=voltage.dtype)
    new_current = np.full(current.shape, np.nan, dtype=current.dtype)

    mask = np.logical_not(mask_zero_diff)
    new_voltage[mask] = voltage[mask]
    new_current[mask] = current[mask]

    # merge zero diff clusters
    while np.any(mask_zero_diff):
        leading_cluster_volt = voltage[mask_zero_diff][0]
        index = np.where(new_voltage == leading_cluster_volt)[0]
        if len(index) == 0:
            index = np.where(mask_zero_diff)[0]

        index = index[0]
        volt_mask = np.isclose(voltage, leading_cluster_volt)

        new_voltage[index] = leading_cluster_volt
        new_current[index] = np.average(current[volt_mask])

        mask_zero_diff[volt_mask] = False

    # remove nan entries
    mask = np.logical_not(np.isnan(new_voltage))
    new_voltage = new_voltage[mask]
    new_current = new_current[mask]

    return new_voltage, new_current


def _merge_voltage_clusters__within_dv(  # noqa: C901
    voltage: np.ndarray,
    current: np.ndarray,
    voltage_step_size: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Take the ``voltage`` and ``current`` arrays associated with a swept
    langmuir trace and average together clusters of voltage values
    within step size ``voltage_step_size``.
    """
    # initialize new voltage and current arrays
    new_voltage = voltage.copy()
    new_current = current.copy()

    voltage_diff = np.diff(voltage)

    # generate global cluster mask
    mask1 = voltage_diff < voltage_step_size
    mask2 = np.isclose(voltage_diff, voltage_step_size)
    cluster_mask = np.logical_or(mask1, mask2)

    # determine cluster locations
    indices = np.where(np.diff(cluster_mask))[0]
    if cluster_mask[0]:
        # 1st element in voltage is in a cluster
        indices = np.append([0], indices)
        indices[1] = indices[1] + 1
    else:
        indices[0:2] = indices[0:2] + 1
    if indices.size > 1:
        indices[2:] = indices[2:] + 1
    if cluster_mask[-1] and indices[-1] != voltage.size - 1:
        # last element in voltage is in a cluster
        indices = np.append(indices, [voltage.size - 1])
    indices = np.reshape(indices, newshape=(int(indices.size / 2), 2))
    # ^ using newshape kwarg for backwards compatibility, newshape kwarg
    #   has been deprecated since numpy v2.1

    # merge clusters and update new_voltage and new_current
    for ii in range(indices.shape[0]):
        start_index = indices[ii][0]
        stop_index = indices[ii][1] + 1

        new_voltage[start_index:stop_index] = np.nan
        new_current[start_index:stop_index] = np.nan

        sub_voltage = voltage[start_index:stop_index]
        sub_current = current[start_index:stop_index]

        v_range = sub_voltage[-1] - sub_voltage[0]
        nbins = int(np.rint(v_range / voltage_step_size))
        if not np.isclose(nbins - (v_range / voltage_step_size), 0):
            nbins = int(
                np.floor((sub_voltage[-1] - sub_voltage[0]) / voltage_step_size)
            )

        if nbins == 0:
            nbins = 1

        if nbins == 1:
            new_voltage[start_index] = np.average(sub_voltage)
            new_current[start_index] = np.average(sub_current)
            continue

        range_array = np.linspace(sub_voltage[0], sub_voltage[-1], nbins + 1)
        for jj in range(nbins):
            start_voltage = range_array[jj]
            stop_voltage = range_array[jj + 1]

            mask1 = sub_voltage >= start_voltage
            if jj == nbins - 1:
                mask2 = sub_voltage <= stop_voltage
            else:
                mask2 = sub_voltage < stop_voltage
            mask = np.logical_and(mask1, mask2)

            if np.count_nonzero(mask) > 0:
                new_voltage[start_index + jj] = np.average(sub_voltage[mask])
                new_current[start_index + jj] = np.average(sub_current[mask])

    # filter out NaN values
    nan_mask = np.logical_not(np.isnan(new_voltage))
    new_voltage = new_voltage[nan_mask]
    new_current = new_current[nan_mask]

    return new_voltage, new_current


def merge_voltage_clusters(  # noqa: C901, PLR0912
    voltage: np.ndarray,
    current: np.ndarray,
    voltage_step_size: float | None = None,
    filter_nan: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Search the ``voltage`` array for closely spaced voltage clusters
    based on the ``voltage_step_size`` parameter and merge those
    clusters, and associated ``current`` values, into a single point.

    This function is intended to merge together identical, or close,
    voltage points to remove small step sizes in the swept langmuir
    trace before differentiation.

    Parameters
    ----------
    voltage: `numpy.ndarray`
        1D `numpy.ndarray` representing the voltage of the swept
        Langmuir trace.  Voltage needs to be monotonically ascending.
        *No units are assumed or checked, but values should be in
        volts.*

    current: `numpy.ndarray`
        1D `numpy.ndarray` representing the current of the swept
        Langmuir trace.  Values should start from a negative
        ion-saturation current and increase to a positive
        electron-saturation current.  *No units are assumed or checked,
        but values should be in amperes.*

    voltage_step_size: `float` | `None`, default: `None`
        A non-zero, positive step size for the ``voltage`` array
        cluster identification.  A value of ``0`` will merge only
        duplicate voltage values.  A value of `None` will default to
        95% of the average step size of the ``voltage`` array (only
        counting duplicate voltages once).

    filter_nan: `bool`, default: `False`
        Set `True` to automatically filter `~numpy.nan` values from the
        ``voltage`` and ``current`` arrays.

    Returns
    -------
    voltage : `numpy.ndarray`
        The new ``voltage`` array.

    current : `numpy.ndarray`
        The new ``current`` array.

    Notes
    -----
    An identified voltage cluster can span a voltage range larger than
    ``voltage_step_size`` and still have every voltage step being
    smaller than ``voltage_step_size``.  In this scenario, the voltage
    cluster will be divided up into :math:`N`-sections for averaging,
    where :math:`N` is given by

    .. math::

        N = \texttt{floor}\left(
            \frac{V_{cluster,max} - V_{cluster,min}}{\texttt{voltage_step_size}}
        \right)
    """
    # condition voltage_step_size
    if voltage_step_size is not None and not isinstance(
        voltage_step_size, numbers.Real
    ):
        raise TypeError(
            "Expected 'voltage_step_size' to be a float or None, got type "
            f"{type(voltage_step_size)}."
        )
    elif isinstance(voltage_step_size, numbers.Integral):
        voltage_step_size = float(voltage_step_size)

    try:
        voltage, current = check_sweep(voltage, current)
    except ValueError as err:
        # check if the ValueError was a result of voltage containing nan
        # values
        if not np.any(np.isnan(voltage)):
            # voltage array has no nan values
            raise ValueError(*err.args) from err

        if not filter_nan:
            raise ValueError(
                "The voltage array contains NaN values.  If you want NaN "
                "values to be automatically filtered, then set argument "
                "filter_nan=True."
            ) from err

        # filter (voltage) NaN values and check sweep again
        nan_mask = np.logical_not(np.isnan(voltage))
        voltage = voltage[nan_mask]
        current = current[nan_mask]

        voltage, current = check_sweep(voltage, current)

    # filter (current) NaN values
    if not np.any(np.isnan(current)):
        # voltage array has no nan values
        pass
    elif filter_nan:
        # mask out current NaN values
        nan_mask = np.logical_not(np.isnan(current))
        voltage = voltage[nan_mask]
        current = current[nan_mask]
    else:
        raise ValueError(
            "The current array contains NaN values.  If you want NaN "
            "values to be automatically filtered, then set argument "
            "filter_nan=True."
        )

    # check if grid is regularly spaced
    voltage_diff = np.diff(voltage)
    mask_zero_diff = np.isclose(voltage_diff, 0.0)
    is_regular_grid = _is_voltage_regularly_spaced(voltage_diff, mask_zero_diff)

    # return if voltage is already regularly spaced and no voltage merging is
    # requested
    if is_regular_grid and voltage_step_size is None:
        warnings.warn(
            "The supplied ``voltage`` array is already regularly spaced. If "
            "you want to re-bin the arrays to a different voltage_step_size, "
            "then use something like numpy.interp.  No merging performed.",
            PlasmaPyWarning,
        )

        return voltage.copy(), current.copy()

    # condition voltage_step_size ... Round 2
    if voltage_step_size is None:
        voltage_step_size = 0.95 * np.abs(
            np.average(voltage_diff[np.logical_not(mask_zero_diff)])
        )
    elif voltage_step_size == 0:
        voltage_step_size = 0.0
    elif voltage_step_size < 0:
        voltage_step_size = -voltage_step_size

    if np.all(voltage_step_size <= np.abs(voltage_diff)) and voltage_step_size != 0:
        warnings.warn(
            f"The supplied voltage_step_size ({voltage_step_size}) is smaller than "
            f"any of the voltage steps in the voltage array.  Supply a step size "
            f"greater than the smallest voltage step ({np.min(voltage_diff)}).  "
            f"No merging performed.",
            PlasmaPyWarning,
        )

        return voltage.copy(), current.copy()

    # now merge clusters
    if voltage_step_size == 0:
        new_voltage, new_current = _merge_voltage_clusters__zero_diff_neighbors(
            voltage, current
        )

    else:
        new_voltage, new_current = _merge_voltage_clusters__within_dv(
            voltage, current, voltage_step_size
        )

    return new_voltage, new_current
