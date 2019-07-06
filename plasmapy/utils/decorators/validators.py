"""
Various decorators to validate input/output arguments to functions.
"""
__all__ = ['ValidateQuantities']

import functools
import inspect
import warnings

from astropy import units as u
from plasmapy.utils.decorators.checks import (CheckValues, CheckUnits)
from plasmapy.utils.decorators.helpers import preserve_signature
from plasmapy.utils.exceptions import ImplicitUnitConversionWarning
from typing import (Any, Dict, List, Union)


class ValidateQuantities(CheckUnits):
    """
    A decorator class to "validate" (i.e. control convert) the units (and values) of
    input/output arguments to a function.  (Validating of function arguments `*args`
    and `**kwargs` is not supported.)

    """

    def __init__(self, equivalencies: Union[None, List] = None,
                 **validations: Dict[str, Any]):

        super().__init__(equivalencies, **validations)
        self._validations = validations

    def __call__(self, f):
        self.f = f
        wrapped_sign = inspect.signature(f)

        @preserve_signature
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            # combine args and kwargs into dictionary
            bound_args = wrapped_sign.bind(*args, **kwargs)
            bound_args.apply_defaults()

            # check and convert units
            unit_validations = self._get_unit_validations(bound_args)

            for arg_name in unit_validations:
                arg = self._validate_unit(bound_args.arguments[arg_name],
                                          arg_name,
                                          **unit_validations[arg_name])
                bound_args.arguments[arg_name] = arg

            return CheckValues(**self._validations)(f)(**bound_args.arguments)
        return wrapper

    def _get_unit_validations(
            self,
            bound_args: inspect.BoundArguments) -> Dict[str, Dict[str, Any]]:
        """
        Blah
        """
        unit_validations = self._get_checks(bound_args)

        # some extra conditioning
        for arg_name in unit_validations:
            # augment 'none_shall_pass' (if needed)
            try:
                unit_validations[arg_name]['none_shall_pass'] = \
                    self._validations[arg_name]['none_shall_pass']
            except KeyError:
                pass

        return unit_validations

    def _validate_unit(self, arg, arg_name, **validations: Dict[str, Any]):

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
        if arg is not None and not isinstance(arg, u.Quantity):
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
                        f"No units are specified for {arg_name} = {arg} "
                        f"in {self.f.__name__}. Assuming units of "
                        f"{validations['units'][0]}. To silence this warning, "
                        f"explicitly pass in an astropy Quantity "
                        f"(e.g. 5. * astropy.units.cm) "
                        f"(see http://docs.astropy.org/en/stable/units/)"
                    ))

        # check units
        arg, unit, equiv = self._check_unit(arg, arg_name, **validations)

        # convert quantity
        if not validations['pass_equivalent_units']:
            if arg.unit != unit:
                arg.to(unit, equivalencies=equiv)

                # if non-standard conversion then warn
                if equiv is not None:
                    warnings.warn(ImplicitUnitConversionWarning)

        return arg
