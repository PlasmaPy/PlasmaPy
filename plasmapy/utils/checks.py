import functools
import inspect

import numpy as np
from astropy import units as u
from ..constants import c
import warnings
from plasmapy.utils.exceptions import RelativityWarning, RelativityError


def check_quantity(validations):
    """
    Raise an exception if an annotated argument in a decorated function
    is not a `~astropy.units.Quantity` with correct units and valid
    numerical values.

    Parameters
    ----------
    validations : `dict`
        Validation dictionary.

    Raises
    ------
    `TypeError`
        If the argument is not a `~astropy.units.Quantity`, units is not
        entirely units or `argname` does not have a type annotation.

    `~astropy.units.UnitConversionError`
        If the argument is not in acceptable units.

    `ValueError`
        If the argument contains `~numpy.nan` or other invalid values as
        determined by the keywords.

    Returns
    -------
    function
        Decorated function.

    Examples
    --------
    >>> from astropy import units as u
    >>> @check_quantity({
    ... "x": {"units": u.m},
    ... "y": {"units": u.s,
    ...       "can_be_negative": False,
    ...       "can_be_complex": True,
    ...       "can_be_inf": False}
    ... })
    ... def func(x: u.m, y: u.s=1 * u.s):
    ...     return x
    ...
    >>> func(1 * u.m)
    <Quantity 1. m>
    >>> func(1 * u.m, 2 * u.m)
    Traceback (most recent call last):
      ...
    astropy.units.core.UnitConversionError: The argument y to func should be a Quantity with the following units: s

    """
    def decorator(f):
        wrapped_sign = inspect.signature(f)
        fname = f.__name__

        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            # combine args and kwargs into dictionary
            bound_args = wrapped_sign.bind(*args, **kwargs)
            bound_args.apply_defaults()
            given_params_values = bound_args.arguments
            given_params = set(given_params_values.keys())

            # names of params to check
            validated_params = set(validations.keys())

            missing_params = [
                param for param in (validated_params - given_params)
            ]

            if len(missing_params) > 0:
                params_str = ", ".join(missing_params)
                raise TypeError(
                    f"Call to {fname} is missing "
                    f"validated params {params_str}")

            for param_to_check, validation_settings in validations.items():
                value_to_check = given_params_values[param_to_check]

                can_be_negative = validation_settings.get(
                    'can_be_negative', True)
                can_be_complex = validation_settings.get(
                    'can_be_complex', False)
                can_be_inf = validation_settings.get(
                    'can_be_inf', True)

                _check_quantity(value_to_check,
                                param_to_check,
                                fname,
                                validation_settings['units'],
                                can_be_negative=can_be_negative,
                                can_be_complex=can_be_complex,
                                can_be_inf=can_be_inf)

            return f(*args, **kwargs)
        return wrapper
    return decorator


def check_relativistic(func=None, betafrac=0.1):
    r"""
    Warns or raises an exception when the output of the decorated
    function is greater than `betafrac` times the speed of light.

    Parameters
    ----------
    func : `function`, optional
        The function to decorate.

    betafrac : float, optional
        The minimum fraction of the speed of light that will raise a
        `~plasmapy.utils.RelativityWarning`. Defaults to 0.1.

    Returns
    -------
    function
        Decorated function.

    Raises
    ------
    TypeError
        If `V` is not a `~astropy.units.Quantity`.

    ~astropy.units.UnitConversionError
        If `V` is not in units of velocity.

    ValueError
        If `V` contains any `~numpy.nan` values.

    ~plasmapy.utils.RelativityError
        If `V` is greater than or equal to the speed of light.

    Warns
    -----
    ~plasmapy.utils.RelativityWarning
        If `V` is greater than or equal to `betafrac` times the speed of light,
        but less than the speed of light.

    Examples
    --------
    >>> from astropy import units as u
    >>> @check_relativistic
    ... def speed():
    ...     return 1 * u.m / u.s

    Passing in a custom `betafrac`:

    >>> @check_relativistic(betafrac=0.01)
    ... def speed():
    ...     return 1 * u.m / u.s

    """
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            return_ = f(*args, **kwargs)
            _check_relativistic(return_, f.__name__, betafrac=betafrac)
            return return_
        return wrapper
    if func:
        return decorator(func)
    return decorator


def _check_quantity(arg, argname, funcname, units, can_be_negative=True,
                    can_be_complex=False, can_be_inf=True):
    """
    Raise an exception if an object is not a `~astropy.units.Quantity`
    with correct units and valid numerical values.

    Parameters
    ----------
    arg : ~astropy.units.Quantity
        The object to be tested.

    argname : str
        The name of the argument to be printed in error messages.

    funcname : str
        The name of the original function to be printed in error messages.

    units : `~astropy.units.Unit` or list of `~astropy.unit.Unit`
        Acceptable units for `arg`.

    can_be_negative : bool, optional
        `True` if the `~astropy.units.Quantity` can be negative,
        `False` otherwise.  Defaults to `True`.

    can_be_complex : bool, optional
        `True` if the `~astropy.units.Quantity` can be a complex number,
        `False` otherwise.  Defaults to `False`.

    can_be_inf : bool, optional
        `True` if the `~astropy.units.Quantity` can contain infinite
        values, `False` otherwise.  Defaults to `True`.

    Raises
    ------
    TypeError
        If the argument is not a `~astropy.units.Quantity` or units is
        not entirely units.

    ~astropy.units.UnitConversionError
        If the argument is not in acceptable units.

    ValueError
        If the argument contains any `~numpy.nan` or other invalid
        values as determined by the keywords.

    Warns
    -----
    UserWarning
        If a `~astropy.units.Quantity` is not provided and unique units
        are provided, a `UserWarning` will be raised and the inputted
        units will be assumed.

    Examples
    --------
    >>> from astropy import units as u
    >>> _check_quantity(4*u.T, 'B', 'f', u.T)

    """

    # TODO: Replace `funcname` with func.__name__?

    if not isinstance(units, list):
        units = [units]

    for unit in units:
        if not isinstance(unit, (u.Unit, u.CompositeUnit, u.IrreducibleUnit)):
            raise TypeError(
                "The keyword 'units' to check_quantity must be "
                "a unit or a list/tuple containing only units.")

    # Create a generic error message

    typeerror_message = (
        f"The argument {argname} to {funcname} should be a Quantity with "
    )

    if len(units) == 1:
        typeerror_message += f"the following units: {str(units[0])}"
    else:
        typeerror_message += "one of the following units: "
        for unit in units:
            typeerror_message += str(unit)
            if unit != units[-1]:
                typeerror_message += ", "

    if isinstance(arg, (u.Unit, u.CompositeUnit, u.IrreducibleUnit)):
        raise TypeError(typeerror_message)

    # Make sure arg is a quantity with correct units

    if not isinstance(arg, (u.Quantity)):
        if len(units) != 1:
            raise TypeError(typeerror_message)
        else:
            try:
                arg = arg * units[0]
            except Exception:
                raise TypeError(typeerror_message)
            else:
                raise UserWarning(
                    f"No units are specified for {argname} in {funcname}. "
                    f"Assuming units of {str(units[0])}."
                )

    in_acceptable_units = []

    for unit in units:
        try:
            arg.unit.to(unit, equivalencies=u.temperature_energy())
        except Exception:
            in_acceptable_units.append(False)
        else:
            in_acceptable_units.append(True)

    if not np.any(in_acceptable_units):
        raise u.UnitConversionError(typeerror_message)

    # Make sure that the quantity has valid numerical values

    valueerror_message = (
        f"The argument {argname} to function {funcname} cannot contain"
    )

    if np.any(np.isnan(arg.value)):
        raise ValueError(f"{valueerror_message} NaNs.")
    elif np.any(np.iscomplex(arg.value)) and not can_be_complex:
        raise ValueError(f"{valueerror_message} complex numbers.")
    elif not can_be_negative and np.any(arg.value < 0):
        raise ValueError(f"{valueerror_message} negative numbers.")
    elif not can_be_inf and np.any(np.isinf(arg.value)):
        raise ValueError(f"{valueerror_message} infs.")


def _check_relativistic(V, funcname, betafrac=0.1):
    r"""
    Warn or raise error for relativistic or superrelativistic
    velocities.

    Parameters
    ----------
    V : ~astropy.units.Quantity
        A velocity.

    funcname : str
        The name of the original function to be printed in the error
        messages.

    betafrac : float
        The minimum fraction of the speed of light that will generate
        a warning.

    Raises
    ------
    TypeError
        If `V` is not a `~astropy.units.Quantity`.

    ~astropy.units.UnitConversionError
        If `V` is not in units of velocity.

    ValueError
        If `V` contains any `~numpy.nan` values.

    RelativityError
        If `V` is greater than or equal to the speed of light.

    Warns
    -----
    ~plasmapy.utils.RelativityWarning
        If `V` is greater than or equal to the specified fraction of the
        speed of light.

    Examples
    --------
    >>> from astropy import units as u
    >>> _check_relativistic(1*u.m/u.s, 'function_calling_this')

    """

    # TODO: Replace `funcname` with func.__name__?

    errmsg = ("V must be a Quantity with units of velocity in"
              "_check_relativistic")

    if not isinstance(V, u.Quantity):
        raise TypeError(errmsg)

    if V.si.unit != u.m / u.s:
        raise u.UnitConversionError(errmsg)

    if np.any(np.isnan(V.value)):
        raise ValueError(f"V includes NaNs in {funcname}")

    beta = np.max(np.abs((V.si / c).value))

    if beta == np.inf:
        raise RelativityError(f"{funcname} is yielding an infinite velocity.")
    elif beta >= 1:
        raise RelativityError(
            f"{funcname} is yielding a velocity that is {str(round(beta, 3))} "
            f"times the speed of light.")
    elif beta >= betafrac:
        warnings.warn(
            f"{funcname} is yielding a velocity that is "
            f"{str(round(beta * 100, 3))}% of the speed of "
            f"light. Relativistic effects may be important.",
        RelativityWarning)
