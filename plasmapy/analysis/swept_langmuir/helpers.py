"""Helper functions for analyzing swept Langmuir traces."""
__all__ = ["check_sweep"]

import numpy as np


def check_sweep(voltage: np.ndarray, current: np.ndarray) -> None:
    """
    Function for checking that the voltage and current arrays are properly
    formatted for analysis by `plasmapy.analysis.swept_langmuir`.

    Parameters
    ----------
    voltage: `numpy.ndarray`
        1D `numpy.ndarray` representing the voltage of the swept Langmuir trace.
        Voltage should be monotonically increasing.

    current: `numpy.ndarray`
        1D `numpy.ndarray` representing the current of the swept Langmuir trace.
        Values should start from a negative ion-saturation current and increase
        to a positive electron-saturation current.

    Raises
    ------
    `TypeError`
        If either the ``voltage`` or ``current`` arrays are not instances of a
        `numpy.ndarray`.

    `ValueError`
        If either the ``voltage`` or ``current`` arrays are not 1D.

        If the ``voltage`` array is not monotonically increasing.

        If the ``current`` array never cross zero (i.e. has not floating
        potential).

        If the ``current`` array does not start form a negative ion-saturation
        current and increases to a positive electron-saturation current.

    """
    # examine voltage array
    if not isinstance(voltage, np.ndarray):
        raise TypeError(
            f"Expected 1D numpy array for voltage, but got {type(voltage)}.",
        )
    elif voltage.ndim != 1:
        raise ValueError(
            f"Expected 1D numpy array for voltage, but got array with "
            f"{voltage.ndim} dimensions.",
        )
    elif not np.all(np.diff(voltage) >= 0):
        raise ValueError(f"The voltage array is not monotonically increasing.")

    # examine current array
    if not isinstance(current, np.ndarray):
        raise TypeError(
            f"Expected 1D numpy array for current, but got {type(current)}.",
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
            "The current array needs to start from a negative ion-saturation current"
            " to a positive electron-saturation current."
        )

    if voltage.size != current.size:
        raise ValueError(
            f"Incompatible arrays, 'voltage' size {voltage.size} must be the same"
            f" as the 'current' size {current.size}."
        )
