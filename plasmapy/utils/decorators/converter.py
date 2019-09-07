"""
Decorator to convert units of functions in /physics methods
"""
__all__ = ["angular_freq_to_hz"]

from astropy import units as u
import functools
import inspect
from plasmapy.utils.decorators import preserve_signature


def angular_freq_to_hz(fn):
    """Decorator to convert angular frequencies from radians per second
    to  Hz.

    Other Parameters
    ----------------
    fn: callable
        The function being wrapped.

    Returns
    -------
    callable
        Wrapped version of the function.
    -------
    """
    wrapped_sign = inspect.signature(fn)
    fname = fn.__name__
    @preserve_signature
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
            # Obtain to_hz value from Original function
            bound_args = wrapped_sign.bind(*args, **kwargs)
            bound_args.apply_defaults()
            given_params_values = bound_args.arguments
            given_params = set(given_params_values.keys())
            # Getting the original value from the original function
            to_hz_value = given_params_values["to_hz"]
            # result of the function
            func_result = fn(*args, **kwargs)
            # if true radians / s will be converted to hz
            if to_hz_value:
                return func_result.to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])
            else:
                return func_result
    return wrapper
