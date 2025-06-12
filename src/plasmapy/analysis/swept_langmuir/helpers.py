"""Helper functions for analyzing swept Langmuir traces."""

__all__ = ["check_sweep", "sort_sweep_arrays"]

import astropy.units as u
import numpy as np

from typing import Literal


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

    index_sort = np.argsort(voltage)
    if voltage_order == "descending":
        index_sort = index_sort[::-1]

    voltage = voltage[index_sort]
    current = current[index_sort]

    return voltage, current
