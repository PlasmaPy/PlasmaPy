"""
Functionality for determining the ion-saturation current of a Langmuir sweep.
"""
__all__ = ["find_ion_saturation_current"]

import numpy as np

from collections import namedtuple

from plasmapy.analysis.swept_langmuir.fit_functions import (
    ExponentialLinearFitFunction,
    LinearFitFunction,
)


IonSaturationCurrentResults = namedtuple(
    "FloatingPotentialResults",
    ("isat_func", "rsq", "func", "indices"),
)


def find_ion_saturation_current(
    voltage: np.ndarray,
    current: np.ndarray,
    upper_bound=None,
    fit_type: str = "explinear",
):
    rtn = IonSaturationCurrentResults(
        isat_func=None,
        rsq=None,
        func=None,
        indices=None,
    )._asdict()

    fit_funcs = {
        "linear": LinearFitFunction(),
        "explinear": ExponentialLinearFitFunction()
    }
    try:
        fit_func = fit_funcs[fit_type]
        rtn["func"] = fit_func
    except KeyError:
        raise KeyError(
            f"Requested fit function '{fit_type}' is  not a valid option, "
            f"expecting on of '{list(fit_funcs)}'.  Examine kwarg 'fit_curve' "
            f"for valid options."
        )

    # mask = np.where(current <= 0)[0]
    if upper_bound is None:
        upper_bound = 0
    mask = np.where(voltage <= upper_bound)[0]
    rtn["indices"] = mask

    volt_sub = voltage[mask]
    curr_sub = current[mask]
    fit_func.curve_fit(volt_sub, curr_sub)

    rtn["rsq"] = fit_func.rsq

    isat = LinearFitFunction()
    isat._parameters = [fit_func.parameters.m, fit_func.parameters.b]
    isat._parameters_err = [fit_func.parameters_err.m, fit_func.parameters_err.b]
    rtn["isat_func"] = isat

    return IonSaturationCurrentResults(**rtn)
