"""
Various decorators to validate input/output arguments to functions.
"""
__all__ = ['validate_values']

import functools
import inspect
import numpy as np

from astropy import units as u
from plasmapy.utils.decorators.helpers import preserve_signature
from typing import (Any, Dict, List, Union)


class ValidateValues:
    """
    Class for defining a decorator that validates input values to a function.
    """

    @classmethod
    def as_decorator(cls,
                     func=None,
                     **validations: Dict[str, Union[str, List, None, u.Quantity]]):
        """
        The decorator for the class.

        Parameters
        ----------
        func:
            The function to be decorated

        **validations: Dict[str, Union[str, List, None, :class:`astropy.units.Quantity`]]
            Input arguments to the wrapped function whose values are to be validated.
            The name of each keyword should be the name of the input argument to the
            wrapped function and its value is a dictionary specifying how the value
            should be validated.  For example, `mass={'can_be_negative': False}` would
            specify the `mass` argument to the wrapped function can not be negative.
            The following keys are allows in the validations:

            ================ ===== ================================================
            Key              Type  Description
            ================ ===== ================================================
            can_be_negative  bool  `True` (DEFAULT) values can be negative
            can_be_complex   bool  `False` (DEFAULT) values can be complex numbers
            can_be_inf       bool  `True` (DEFAULT) values can be infinite
            can_be_nan       bool  `True` (DEFAULT) values can be NaN
            none_shall_pass  bool  `False` (DEFAULT) values can be python `None`
            ================ ===== ================================================

        """
        # if func is not None and not validations:
        if func is not None:
            return cls(**validations)(func)
        else:
            return cls(**validations)

    def __init__(self, **validations: Dict[str, Union[bool, List, None, u.Quantity]]):
        self.validations = validations

    def __call__(self, f):
        self.f = f
        wrapped_sign = inspect.signature(f)

        @preserve_signature
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            # combine args and kwargs into dictionary
            bound_args = wrapped_sign.bind(*args, **kwargs)
            bound_args.apply_defaults()
            given_args = bound_args.arguments

            # get validations
            validations = self._get_validations(bound_args)

            # Does `validations` indicate arguments not used by f?
            missing_params = [
                param
                for param in set(validations.keys()) - set(given_args.keys())
            ]
            if len(missing_params) > 0:
                params_str = ", ".join(missing_params)
                raise TypeError(
                    f"Call to {self.f.__name__} is missing validated "
                    f"params {params_str}")

            # Review validations and check argument
            for arg_name in validations:
                self._validate_value(given_args[arg_name],
                                     arg_name,
                                     **validations[arg_name])

            return f(**given_args)

        return wrapper

    def _get_validations(self, bound_args: inspect.BoundArguments) -> Dict[str, Any]:
        # initialize validation dictionary
        out_validations = {}

        # Iterate through parameters, determine validation keys, and build validations
        # dictionary `value_validations`
        for param in bound_args.signature.parameters.values():
            # variable arguments (*args, **kwargs) not validted
            if param.kind in (inspect.Parameter.VAR_KEYWORD,
                              inspect.Parameter.VAR_POSITIONAL):
                continue

            # catch the (never triggered) case where bind relied on a default value
            if param.name not in bound_args.arguments \
                    and param.default is not param.empty:
                bound_args.arguments[param.name] = param.default

            # grab the validations dictionary for the desired parameter
            try:
                param_in_validations = self.validations[param.name]
            except KeyError:
                # validations for parameter not specified
                continue

            # build `value_validations`
            # read validations and/or apply defaults values
            out_validations[param.name] = {}
            for v_name, v_default in self._validation_item_defaults.items():
                out_validations[param.name][v_name] = \
                    param_in_validations.get(v_name, v_default)

        return out_validations

    @property
    def _validation_item_defaults(self) -> Dict[str, bool]:
        # adding a new validation key and default value here will automatically
        # update the checks of the `validations` passed in
        _defaults = {
            'can_be_negative': True,
            'can_be_complex': False,
            'can_be_inf': True,
            'can_be_nan': True,
            'none_shall_pass': False,
        }
        return _defaults

    def _validate_value(self, arg, arg_name, **validations):
        valueerror_msg = (f"The argument '{arg_name}'' to function "
                          f"{self.f.__name__}() can not contain")

        # check values
        if arg is None and validations['none_shall_pass']:
            return
        elif arg is None:
            raise ValueError(f"{valueerror_msg} Nones.")
        elif not validations['can_be_negative']:
            # Allow NaNs through without raising a warning
            with np.errstate(invalid='ignore'):
                isneg = np.any(arg < 0)
            if isneg:
                raise ValueError(f"{valueerror_msg} negative numbers.")
        elif not validations['can_be_complex'] and np.any(np.iscomplexobj(arg)):
            raise ValueError(f"{valueerror_msg} complex numbers.")
        elif not validations['can_be_inf'] and np.any(np.isinf(arg)):
            raise ValueError(f"{valueerror_msg} infs.")
        elif not validations['can_be_nan'] and np.any(np.isnan(arg)):
            raise ValueError(f"{valueerror_msg} NaNs.")


validate_values = ValidateValues.as_decorator
