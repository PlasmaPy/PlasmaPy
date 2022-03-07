"""
Functionality for determining the ion-saturation current of a Langmuir sweep.
"""
__all__ = ["find_ion_saturation_current", "ISatExtras"]
__aliases__ = ["find_isat_"]

import numbers
import numpy as np

from typing import Tuple, NamedTuple, Union

from plasmapy.analysis import fit_functions as ffuncs
from plasmapy.analysis.swept_langmuir.helpers import check_sweep

__all__ += __aliases__


class ISatExtras(NamedTuple):
    """
    Create a tuple containing the extra parameters calculated by
    `find_ion_saturation_current`.
    """
    rsq: Union[float, None]
    """
    Alias for field number 0, the r-squared value of the ion-saturation
    curve fit.
    """

    fitted_func:  Union[ffuncs.AbstractFitFunction, None]
    """
    Alias for field number 1, the :term:`fit-function` fitted during
    the ion-saturation curve fit.
    """

    fitted_indices: Union[slice, None]
    """
    Alias for field number 2, the indices used in the ion-saturation
    curve fit.
    """


def find_ion_saturation_current(
    voltage: np.ndarray,
    current: np.ndarray,
    *,
    fit_type: str = "exp_plus_linear",
    upper_bound: float = None,
) -> Tuple[ffuncs.Linear, ISatExtras]:
    """
    Determines the ion-saturation current (:math:`I_{sat}`) for a given
    current-voltage (IV) curve obtained from a swept Langmuir probe.
    The current collected by a Langmuir probe reaches ion-saturation
    when the probe is sufficiently biased so the influx of electrons is
    completely repelled leading to only the collection of ions.  (For
    additional details see the **Notes** section below.)

    **Aliases:** `find_isat_`

    Parameters
    ----------

    voltage: `numpy.ndarray`
        1-D numpy array of monotonically ascending/descending probe
        biases (should be in volts).

    current: `numpy.ndarray`
        1-D numpy array of probe current (should be in amperes)
        corresponding to the ``voltage`` array.

    fit_type: `str`
        The type of curve (:term:`fit-function`) to be fitted to the
        Langmuir trace, valid options are listed below.
        (DEFAULT ``"exp_plus_linear"``)

        +-----------------------+----------------------------------------------------------+
        | ``"linear"``          | `~plasmapy.analysis.fit_functions.Linear`                |
        +-----------------------+----------------------------------------------------------+
        | ``"exponential"``     | `~plasmapy.analysis.fit_functions.ExponentialPlusOffset` |
        +-----------------------+----------------------------------------------------------+
        | ``"exp_plus_linear"`` | `~plasmapy.analysis.fit_functions.ExponentialPlusLinear` |
        +-----------------------+----------------------------------------------------------+

    upper_bound: `float`
        A bias voltage (in volts) that specifies an upper bound used to
        collect the points for the curve fit.  That is, points that
        satisfy ``voltage <= upper_bound`` are used in the fit.
        (DEFAULT ``None``)

        |

        If ``upper_bound`` is not specified, then the route will collect
        indices based on a current bound relative to ``np.min(current)``
        that is dependent on the ``fit_type`` specified.

    Returns
    -------
    isat: `~plasmapy.analysis.fit_functions.Linear`
        A :term:`fit-function` representing the linear portion of the
        fitter curve.

    extras: `ISatExtras`
        Additional information from the curve fit:

        * ``extras.fitted_func`` is the :term:`fit-function` (specified
          by ``fit_type``) fitted to the IV-curve
        * ``extras.rsq`` is the coefficient of determination
          (r-squared) value of the ``extras.fitted_func`` to the IV-curve
        * ``extras.fitted_indices`` is a `slice` object representing the
          points used in the curve fit (i.e.
          ``(voltage[extras.fitted_indices], current[extras.fitted_indices])``).
    """
    rtn_extras = ISatExtras(rsq=None, fitted_func=None, fitted_indices=None)._asdict()

    _settings = {
        "linear": {
            "func": ffuncs.Linear,
            "default_ubound_frac": 0.4,
        },
        "exp_plus_linear": {
            "func": ffuncs.ExponentialPlusLinear,
            "default_ubound_frac": 1.0,
        },
        "exp_plus_offset": {
            "func": ffuncs.ExponentialPlusOffset,
            "default_ubound_frac": 1.0,
        },
    }
    try:
        default_ubound_frac = _settings[fit_type]["default_ubound_frac"]
        fit_func = _settings[fit_type]["func"]()
        rtn_extras["fitted_func"] = fit_func
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
    rtn_extras["fitted_indices"] = mask

    volt_sub = voltage[mask]
    curr_sub = current[mask]
    fit_func.curve_fit(volt_sub, curr_sub)

    rtn_extras["rsq"] = fit_func.rsq

    m = getattr(fit_func.params, "m", 0.0)
    b = getattr(fit_func.params, "b", 0.0)
    m_err = getattr(fit_func.param_errors, "m", 0.0)
    b_err = getattr(fit_func.param_errors, "b", 0.0)
    isat = ffuncs.Linear(params=(m, b), param_errors=(m_err, b_err))

    return isat, ISatExtras(**rtn_extras)


find_isat_ = find_ion_saturation_current
"""Alias to :func:`find_ion_saturation_current`."""
