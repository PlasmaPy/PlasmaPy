import numpy as np
from astropy import units as u
from ..constants import c
import inspect
from functools import wraps


def check_quantity(validations):
    def decorator(f):
        @wraps(f)
        def wrapped_function(*args, **kwargs):
            fname = f.__name__
            fsig = inspect.signature(f)

            # names of params to check
            validated_params = set(validations.keys())

            # map of given param names to given values
            args_params_values = {k:v for (k,v) in zip(fsig.parameters, args)}

            # combine args and kwargs
            given_params_values = {**args_params_values, **kwargs}

            # given param names, including kwargs
            given_params = set(given_params_values.keys())

            default_params_values = {name:param.default for (name, param) in fsig.parameters.items() if param.default != inspect.Parameter.empty}
            default_params = set(default_params_values.keys())

            missing_params = [param for param in (validated_params - given_params - default_params)\
                                    if not validations[param].get('ignore', False)]
            if len(missing_params) > 0:
                raise TypeError("Call to {} is missing validated params {}".format(fname, ", ".join(missing_params)))

            for param_to_check in validated_params:
                if validations[param_to_check].get('ignore', False):
                    continue

                value_to_check = given_params_values.get(param_to_check, default_params_values.get(param_to_check))

                if value_to_check is None:
                    raise TypeError("Parameter {} in call to {} was expected to be validated by check_quantity, "
                                        "but not given".format(param_to_check, fname))

                validation_settings = validations[param_to_check]

                can_be_negative = validation_settings.get('can_be_negative', True)
                can_be_complex = validation_settings.get('can_be_complex', False)
                can_be_inf = validation_settings.get('can_be_inf', True)

                _check_quantity(value_to_check,
                                param_to_check,
                                fname,
                                validation_settings['units'],
                                can_be_negative=can_be_negative,
                                can_be_complex=can_be_complex,
                                can_be_inf=can_be_inf)

            # if we are validating our params, make sure we mention all Quantity params in our validation
            for param in given_params:
                if given_params_values[param] is not None and\
                        issubclass(type(given_params_values[param]), u.quantity.Quantity) and\
                        param not in validated_params:
                    raise ValueError("Parameter {} in call to {} is a Quantity "
                                        "but is not validated".format(param, fname))

            return f(*args, **kwargs)
        return wrapped_function
    return decorator


def _check_quantity(arg, argname, funcname, units, can_be_negative=True,
                    can_be_complex=False, can_be_inf=True):
    """Raises exceptions if an object is not an astropy Quantity with
    correct units and valid numerical values.

    Parameters
    ----------
    arg : Quantity
        The object to be tested.

    argname : string
        The name of the argument to be printed in error messages.

    funcname : string
        The name of the original function to be printed in error messages.

    units : Unit or list of Units
        Acceptable units for arg.

    can_be_negative : boolean, optional
        True if the Quantity can be negative, False otherwise.
        Defaults to True.

    can_be_complex : boolean, optional
        True if the Quantity can be a complex number, False otherwise.
        Defaults to False.

    can_be_inf : boolean, optional
        True if the Quantity can contain infinite values, False
        otherwise.  Defaults to True.

    Raises
    ------
    TypeError
        If the argument is not a Quantity or units is not entirely units.

    UnitConversionError
        If the argument is not in acceptable units.

    ValueError
        If the argument contains NaNs or other invalid values as
        determined by the keywords.

    Examples
    --------
    >>> from astropy import units as u
    >>> check_quantity(4*u.T, 'B', 'f', u.T)

    """

    if not isinstance(units, list):
        units = [units]

    for unit in units:
        if not isinstance(unit, (u.Unit, u.CompositeUnit, u.IrreducibleUnit)):
            raise TypeError("The keyword 'units' to check_quantity must be "
                            "a unit or a list/tuple containing only units.")

    # Create a generic error message

    typeerror_message = ("The argument " + argname + " to " + funcname +
                         " should be a Quantity with ")

    if len(units) == 1:
        typeerror_message += "the following units: " + str(units[0])
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
                arg = arg*units[0]
            except Exception:
                raise TypeError(typeerror_message)
            else:
                raise UserWarning("No units are specified for " + argname +
                                  " in " + funcname + ". Assuming units of " +
                                  str(units[0]) + ".")

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

    valueerror_message = ("The argument " + argname + " to function " +
                          funcname + " cannot contain ")

    if np.any(np.isnan(arg.value)):
        raise ValueError(valueerror_message + "NaNs.")
    elif np.any(np.iscomplex(arg.value)) and not can_be_complex:
        raise ValueError(valueerror_message + "complex numbers.")
    elif not can_be_negative and np.any(arg.value < 0):
        raise ValueError(valueerror_message + "negative numbers.")
    elif not can_be_inf and np.any(np.isinf(arg.value)):
        raise ValueError(valueerror_message + "infs.")


def check_relativistic(betafrac=0.1):
    def _check_relativistic(V, funcname, betafrac=0.1):
        r"""Raise UserWarnings if a velocity is relativistic or superrelativistic

        Parameters
        ----------
        V : Quantity
            A velocity

        funcname : string
            The name of the original function to be printed in the error messages.

        betafrac : float
            The minimum fraction of the speed of light that will raise a
            UserWarning

        Raises
        ------
        TypeError
            If V is not a Quantity

        UnitConversionError
            If V is not in units of velocity

        ValueError
            If V contains any NaNs

        UserWarning
            If V is greater than betafrac times the speed of light

        Examples
        --------
        >>> from astropy import units as u
        >>> _check_relativistic(1*u.m/u.s, 'function_calling_this')

        """

        errmsg = ("V must be a Quantity with units of velocity in"
                  "_check_relativistic")

        if not isinstance(V, u.Quantity):
            raise TypeError(errmsg)

        if V.si.unit != u.m/u.s:
            raise u.UnitConversionError(errmsg)

        if np.any(np.isnan(V.value)):
            raise ValueError("V includes NaNs in _check_relativistic")

        beta = np.max(np.abs((V/c).value))

        if beta == np.inf:
            raise UserWarning(funcname + " is yielding an infinite velocity.")
        elif beta >= 1:
            raise UserWarning(funcname + " is yielding a velocity that is " +
                              str(round(beta, 3)) + " times the speed of light.")
        elif beta >= betafrac:
            raise UserWarning(funcname + " is yielding a velocity that is " +
                              str(round(beta*100, 3)) + "% of the speed of " +
                              "light. Relativistic effects may be important.")

    def decorator(f):
        @wraps(f)
        def wrapped_function(*args, **kwargs):
            result = f(*args, **kwargs)
            _check_relativistic(result, f.__name__, betafrac)
            return result
        return wrapped_function
    return decorator


