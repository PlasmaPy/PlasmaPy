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

    # generate data slice

    # fit data

    # calculated Te and error


find_te_ = find_electron_temperature
"""
Alias to
:func:`~plasmapy.analysis.swept_langmuir.electron_temperature.find_electron_temperature`.
"""
