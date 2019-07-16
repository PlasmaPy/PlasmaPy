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
from plasmapy.utils.exceptions import (ImplicitUnitConversionWarning, PlasmaPyWarning)
from typing import (Any, Dict, List, Union)


class ValidateQuantities(CheckUnits, CheckValues):
    """
    A decorator class to "validate" (i.e. control convert) the units (and values) of
    input/output arguments to a function.  (Validating of function arguments `*args`
    and `**kwargs` is not supported.)

    """

    def __init__(self, validations_on_return=None, **validations: Dict[str, Any]):

        checks_on_return = validations.pop('checks_on_return', None)
        if checks_on_return is not None:
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
                if arg_name == 'checks_on_return':
                    continue

                # validate argument & update for conversion
                arg = self._validate_quantity(bound_args.arguments[arg_name],
                                              arg_name,
                                              **validations[arg_name])
                bound_args.arguments[arg_name] = arg

            # call function
            _return = f(**bound_args.arguments)

            # validate output
            if 'checks_on_return' in validations:
                _return = self._validate_quantity(_return, 'checks_on_return',
                                                  **validations['checks_on_return'])

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

        # ensure both dictionaries have the exact same keys
        # key_set = set(list(value_validations.keys()) + list(validations))
        # if len(key_set - set(value_validations.keys())) != 0 \
        #         or len(key_set - set(validations.keys())) != 0:
        #     raise TypeError

        # combine all validations
        key_set = set(list(unit_checks.keys()) + list(value_checks.keys()))
        validations = {}
        for arg_name in key_set:
            # add unit checks to validations
            try:
                validations[arg_name] = unit_checks[arg_name].copy()
            except KeyError:
                # no unit checks were defined, but value checks were
                validations[arg_name] = value_checks[arg_name]

                # TODO: improve warning string
                warnings.warn(u.UnitsWarning(f"No units defined for argument"))

                # jump to next for-loop iteration
                continue

            # augment 'none_shall_pass' (if needed)
            try:
                # if 'none_shall_pass' was in the original passed-in validations,
                # then override the value determined by CheckUnits
                validations[arg_name]['none_shall_pass'] = \
                    self._validations[arg_name]['none_shall_pass']
            except KeyError:
                # 'none_shall_pass' was not in the original passed-in validations, so
                # rely on the value determined by CheckUnits
                pass
            else:
                del value_checks[arg_name]['none_shall_pass']

            # update the validations dictionary
            try:
                validations[arg_name].update(value_checks[arg_name])
            except KeyError:
                # no value checks were defined
                pass

        return validations

    def _validate_quantity(self, arg, arg_name, **validations: Dict[str, Any]):

        # initialize TypeError message
        typeerror_msg = (
            f"The argument {arg_name} to {self.f.__name__} should be an astropy"
            f"Quantity with units equivalent to one of ["
        )
        for ii, unit in enumerate(validations['units']):
            typeerror_msg += f"{unit}"

            if ii != len(validations['units']) - 1:
                typeerror_msg += f", "
        typeerror_msg += f"]"

        # add units to arg if possible
        # * a None value will be taken care of by `_check_unit`
        #
        if 'units' not in validations:
            pass
        elif arg is not None and not isinstance(arg, u.Quantity):
            if len(validations['units']) != 1:
                raise TypeError(typeerror_msg)
            else:
                try:
                    arg = arg * validations['units'][0]
                except (TypeError, ValueError):
                    raise TypeError(typeerror_msg)
                else:
                    if not isinstance(arg, u.Quantity):  # pragma: no cover
                        # this should never be reached...if so, except
                        # is not setup correctly
                        #
                        raise TypeError(typeerror_msg)

                    warnings.warn(u.UnitsWarning(
                        f"No units are specified for {arg_name} = {arg.value} "
                        f"in {self.f.__name__}(). Assuming units of "
                        f"{validations['units'][0]}. To silence this warning, "
                        f"explicitly pass in an astropy Quantity "
                        f"(e.g. 5. * astropy.units.cm) "
                        f"(see http://docs.astropy.org/en/stable/units/)"
                    ))

        if 'units' in validations:
            # check units
            arg, unit, equiv = self._check_unit(arg, arg_name, **validations)

            # convert quantity
            if arg is None:
                pass
            elif unit is not None and unit != arg.unit:
                arg = arg.to(unit, equivalencies=equiv)

                # if non-standard conversion then warn
                if equiv is not None:
                    warnings.warn(ImplicitUnitConversionWarning)

        # check value
        if all([key in validations for key in self._CheckValues__check_defaults]):
            self._check_value(arg, arg_name, **validations)

        return arg


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
