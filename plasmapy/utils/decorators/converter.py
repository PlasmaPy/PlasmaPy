"""
Decorator to convert units of functions in /physics methods
"""
from astropy import units as u
import functools
import inspect
from plasmapy.utils.decorators import preserve_signature
__all__ = ["from_radians_to_hz"]

def from_radians_to_hz(fn):
    """Decorator to convert angular frequencies from units of radians per second
     to  Hz.
        Parameters
    ----------
    """
    wrapped_sign = inspect.signature(fn)
    fname = fn.__name__
    @preserve_signature
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
            #Obtain to_hz value from Original function
            bound_args = wrapped_sign.bind(*args, **kwargs)
            bound_args.apply_defaults()
            given_params_values = bound_args.arguments
            given_params = set(given_params_values.keys())
            #Getting the original value from the origina function
            to_hz_value = given_params_values["to_hz"]
            #result of the function
            func_result = fn(*args, **kwargs)
            #if true radians / s will be converted to hz
            if to_hz_value:
                return func_result.to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])
            else:
                return func_result
    return wrapper
