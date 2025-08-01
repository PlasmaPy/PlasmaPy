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
    ...


find_te_ = find_electron_temperature
"""
Alias to
:func:`~plasmapy.analysis.swept_langmuir.electron_temperature.find_electron_temperature`.
"""
