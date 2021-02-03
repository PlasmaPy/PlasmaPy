"""
Functionality for determining the ion-saturation current of a Langmuir sweep.
"""
__all__ = ["find_ion_saturation_current", "find_isat_"]

import numbers
import numpy as np

from collections import namedtuple

from plasmapy.analysis import fit_functions as ffuncs
from plasmapy.analysis.swept_langmuir.helpers import check_sweep

IonSaturationCurrentResults = namedtuple(
    "FloatingPotentialResults",
    ("isat_func", "rsq", "func", "indices"),
)


def find_ion_saturation_current(
    voltage: np.ndarray,
    current: np.ndarray,
    upper_bound: float = None,
    fit_type: str = "exp_plus_linear",
):
    rtn = IonSaturationCurrentResults(
        isat_func=None, rsq=None, func=None, indices=None,
    )._asdict()

    _settings = {
        "linear": {
            "func": ffuncs.Linear, "default_ubound_frac": 0.4,
        },
        "exp_plus_linear": {
            "func": ffuncs.ExponentialPlusLinear, "default_ubound_frac": 1.0,
        },
        "exp_plus_offset": {
            "func": ffuncs.ExponentialPlusOffset, "default_ubound_frac": 1.0,
        },
    }
    try:
        default_ubound_frac = _settings[fit_type]["default_ubound_frac"]
        fit_func = _settings[fit_type]["func"]()
        rtn["func"] = fit_func
    except KeyError:
        raise ValueError(
            f"Requested fit '{fit_type}' is not a valid option.  Valid options "
            f"are {list(_settings.keys())}."
        )

    # check voltage and current arrays
    voltage, current = check_sweep(voltage, current, strip_units=True)

    # condition kwarg upper_bound
    if upper_bound is None:
        current_min = current.min()
        current_bound = (1.0 - default_ubound_frac) * current_min
        mask = np.where(current <= current_bound)[0]

    elif not isinstance(upper_bound, numbers.Real):
        raise TypeError(
            f"Keyword 'upper_bound' is of type {type(upper_bound)}, expected an "
            f"int or float."
        )
    else:
        mask = np.where(voltage <= upper_bound)[0]
        if mask.size == 0:
            raise ValueError(
                f"The value for keyword 'upper_bound' ({upper_bound}) resulted in "
                f"identifying a fit window of size 0."
            )
    mask = slice(0, mask[-1] + 1)
    rtn["indices"] = mask

    volt_sub = voltage[mask]
    curr_sub = current[mask]
    fit_func.curve_fit(volt_sub, curr_sub)

    rtn["rsq"] = fit_func.rsq

    m = getattr(fit_func.params, "m", 0.0)
    b = getattr(fit_func.params, "b", 0.0)
    m_err = getattr(fit_func.param_errors, "m", 0.0)
    b_err = getattr(fit_func.param_errors, "b", 0.0)
    isat = ffuncs.Linear(params=(m, b), param_errors=(m_err, b_err))
    rtn["isat_func"] = isat

    return IonSaturationCurrentResults(**rtn)


find_isat_ = find_ion_saturation_current
"""Alias to :func:`find_ion_saturation_current`."""
