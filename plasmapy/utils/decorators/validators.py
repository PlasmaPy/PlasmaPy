"""
Various decorators to validate input/output arguments to functions.
"""
__all__ = ['validate_quantities', 'ValidateQuantities']

import functools
import inspect
import warnings

from astropy import units as u
from plasmapy.utils.decorators.checks import (CheckValues, CheckUnits)
from plasmapy.utils.decorators.helpers import preserve_signature
from plasmapy.utils.exceptions import ImplicitUnitConversionWarning
from typing import (Any, Dict, List, Union)


class ValidateQuantities(CheckUnits, CheckValues):
    """
    A decorator class to "validate" (i.e. control convert) the units (and values) of
    input/output arguments to a function.  (Validating of function arguments `*args`
    and `**kwargs` is not supported.)

    """

    def __init__(self, validations_on_return=None, **validations: Dict[str, Any]):

        if 'checks_on_return' in validations:
            raise TypeError(f"keyword argument 'checks_on_return' is not allowed, "
                            f"use 'validations_on_return' to set validations "
                            f"on the return variable")

        self._validations = validations

        checks = validations.copy()
        if validations_on_return is not None:
            self._validations['validations_on_return'] = validations_on_return
            checks['checks_on_return'] = validations_on_return

        super().__init__(**checks)

    def __call__(self, f):
        self.f = f
        wrapped_sign = inspect.signature(f)

        @preserve_signature
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            # combine args and kwargs into dictionary
            bound_args = wrapped_sign.bind(*args, **kwargs)
            bound_args.apply_defaults()

            # get conditioned validations
            validations = self._get_validations(bound_args)

            # validate (input) argument units and values
            for arg_name in validations:
                # skip check of output/return
                if arg_name == 'validations_on_return':
                    continue

                # validate argument & update for conversion
                arg = self._validate_quantity(bound_args.arguments[arg_name],
                                              arg_name,
                                              **validations[arg_name])
                bound_args.arguments[arg_name] = arg

            # call function
            _return = f(**bound_args.arguments)

            # validate output
            if 'validations_on_return' in validations:
                _return = self._validate_quantity(_return, 'validations_on_return',
                                                  **validations['validations_on_return'])

            return _return
        return wrapper

    def _get_validations(
            self,
            bound_args: inspect.BoundArguments) -> Dict[str, Dict[str, Any]]:
        """
        Blah
        """
        unit_checks = self._get_unit_checks(bound_args)
        value_checks = self._get_value_checks(bound_args)

        # combine all validations
        # * `unit_checks` will encompass all argument "checks" defined either by
        #   function annotations or **validations.
        # * `value_checks` may miss some arguments if **validations only defines
        #   unit validations or some validations come from function annotations
        #
        validations = unit_checks.copy()
        for arg_name in validations:
            # augment 'none_shall_pass' (if needed)
            try:
                # if 'none_shall_pass' was in the original passed-in validations,
                # then override the value determined by CheckUnits
                validations[arg_name]['none_shall_pass'] = \
                    self.validations[arg_name]['none_shall_pass']
            except (KeyError, TypeError):
                # 'none_shall_pass' was not in the original passed-in validations, so
                # rely on the value determined by CheckUnits
                pass
            finally:
                try:
                    del value_checks[arg_name]['none_shall_pass']
                except KeyError:
                    dvc = self._CheckValues__check_defaults.copy()
                    del dvc['none_shall_pass']
                    value_checks[arg_name] = dvc

            # update the validations dictionary
            validations[arg_name].update(value_checks[arg_name])

        if 'checks_on_return' in validations:
            validations['validations_on_return'] = validations.pop('checks_on_return')

        return validations

    def _validate_quantity(self, arg, arg_name, **validations):
        # rename to work with "check" methods
        if arg_name == 'validations_on_return':
            arg_name = 'checks_on_return'

        # initialize str for error message
        if arg_name == 'checks_on_return':
            err_msg = f"The return value  "
        else:
            err_msg = f"The argument '{arg_name}' "
        err_msg += f"to function {self.f.__name__}()"

        # initialize TypeError message
        typeerror_msg = (f"{err_msg} should be an astropy Quantity with units"
                         f" equivalent to one of [")
        for ii, unit in enumerate(validations['units']):
            typeerror_msg += f"{unit}"

            if ii != len(validations['units']) - 1:
                typeerror_msg += f", "
        typeerror_msg += f"]"

        # add units to arg if possible
        # * a None value will be taken care of by `_check_unit_core`
        #
        if arg is None or hasattr(arg, 'unit'):
            pass
        elif len(validations['units']) != 1:
            raise TypeError(typeerror_msg)
        else:
            try:
                arg = arg * validations['units'][0]
            except (TypeError, ValueError):
                raise TypeError(typeerror_msg)
            else:
                warnings.warn(u.UnitsWarning(
                    f"{err_msg} has no specified units. Assuming units of "
                    f"{validations['units'][0]}. To silence this warning, "
                    f"explicitly pass in an astropy Quantity "
                    f"(e.g. 5. * astropy.units.cm) "
                    f"(see http://docs.astropy.org/en/stable/units/)"
                ))

        # check units
        arg, unit, equiv, err = self._check_unit_core(arg, arg_name, validations)

        # convert quantity
        if arg is not None and unit is not None:
            # if non-standard conversion then warn
            if not arg.unit.is_equivalent(unit, equivalencies=None):
                warnings.warn(ImplicitUnitConversionWarning(
                    f"{err_msg} has a non-standard unit conversion..."
                    f"converting {arg.unit} to {unit}"))

            arg = arg.to(unit, equivalencies=equiv)
        elif err is not None:
            raise err

        # check value
        self._check_value(arg, arg_name, validations)

        return arg

    @property
    def validations(self):
        """
        Requested validations on the decorated function's input arguments and
        return variable.
        """
        return self._validations


def validate_quantities(func=None, validations_on_return=None, **validations):
    """
    Parameters
    ----------
    func
    validations_on_return
    validations

    Returns
    -------

    """
    if validations_on_return is not None:
        validations['validations_on_return'] = validations_on_return

    if func is not None:
        return ValidateQuantities(**validations)(func)
    else:
        return ValidateQuantities(**validations)
