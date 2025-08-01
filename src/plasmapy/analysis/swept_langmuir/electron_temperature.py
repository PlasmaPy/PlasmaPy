"""
Functionality for determining the electron temperature of a given
Langmuir sweep.
"""
__all__ = ["TeExtras", "find_electron_temperature"]
__aliases__ = ["find_te_"]
__all__ += __aliases__

import numbers
from typing import Callable, NamedTuple, Sequence

import numpy as np

from plasmapy.analysis import fit_functions as ffuncs
from plasmapy.analysis.swept_langmuir.helpers import check_sweep
from plasmapy.analysis.swept_langmuir.plasma_potential import (
    _condition_voltage_window
)


class TeExtras(NamedTuple):
    err: float | None
    rsq: float | None
    fitted_func:  ffuncs.AbstractFitFunction | None
    fitted_indices: slice | None


def _condition_current_threshold_factors(
    current_threshold_factors: list[float | None] | None
) -> tuple[float, float]:

    thresholds = current_threshold_factors

    if thresholds is None:
        return 0.0, 1.0

    if isinstance(thresholds, np.ndarray):
        thresholds = thresholds.tolist()

    if not isinstance(thresholds, Sequence):
        raise TypeError(
            f"Expected a 2-element list of floats or None for "
            f"'current_threshold_factors', but got type "
            f"{type(thresholds)}."
        )
    elif len(thresholds) != 2:
        raise ValueError(
            f"Expected a 2-element list of floats or None for "
            f"'current_threshold_factors', but got type "
            f"{len(thresholds)} elements."
        )
    elif not all(
        isinstance(element, numbers.Real) or element is None
        for element in thresholds
    ):
        raise TypeError(
            f"Not all elements of 'current_threshold_factors' are floats or None."
        )
    elif None not in thresholds:
        thresholds = np.sort(thresholds).tolist()
    else:
        thresholds = list(thresholds)

    if thresholds[0] is None:
        thresholds[0] = 0.0
    elif not (0 <= thresholds[0] <= 1):
        raise ValueError(
            "The lower current_threshold_factors needs to be in the "
            f"interval [0, 1], got {thresholds[0]}."
        )

    if thresholds[1] is None:
        thresholds[1] = 1.0
    elif not (0 <= thresholds[1] <= 1):
        raise ValueError(
            "The upper current_threshold_factors needs to be in the "
            f"interval [0, 1], got {thresholds[1]}."
        )

    if thresholds[1] < thresholds[0]:
        thresholds = np.sort(thresholds).tolist()
    elif thresholds[0] == thresholds[1]:
        raise ValueError(
            "The upper and lower current_threshold_factors can not be equal."
        )

    return thresholds[0], thresholds[1]


def find_electron_temperature(
    voltage: np.ndarray,
    current: np.ndarray,
    *,
    isat: numbers.Real | Callable | None,
    voltage_window: list[numbers.Real | None] | None = None,
    current_threshold_factors: list[float | None] | None = (0.05, 1.0),
):
    """threshhold"""
    rtn_extras = TeExtras(
        err=None,
        rsq=None,
        fitted_func=None,
        fitted_indices=None,
    )

    # check voltage and current arrays
    voltage, current = check_sweep(voltage, current, strip_units=True)

    # condition voltage_window
    voltage_window = _condition_voltage_window(voltage, voltage_window)
    data_size = voltage[voltage_window].size
    if data_size < 3:
        raise ValueError(
            f"The specified voltage_window ({voltage_window}) would result "
            f"in a null window or a window with fewer than 3 elements."
        )

    # condition isat
    if isat is None:
        # this indicates the user has already removed isat from the
        # given current array
        current_adjusted = current.copy()
    elif isinstance(isat, numbers.Real):
        current_adjusted = current - isat
    elif callable(isat):
        try:
            current_adjusted = current - isat(voltage)
        except Exception as err:
            raise RuntimeError(
                "Supplied isat callable does not have appropriate "
                "signature.  isat should take a voltage array and "
                "produce an isat current array that can be subtracted "
                "from the probe current array, i.e. "
                "`current - isat(voltage)`."
            ) from err
    else:
        raise TypeError(
            f"Expected None, float, or a callable for argument isat, "
            f"but got {type(isat)}."
        )

    # condition current_threshold_factors
    thresholds = _condition_current_threshold_factors(current_threshold_factors)

    # generate data slice

    # fit data

    # calculated Te and error


find_te_ = find_electron_temperature
"""
Alias to
:func:`~plasmapy.analysis.swept_langmuir.electron_temperature.find_electron_temperature`.
"""
