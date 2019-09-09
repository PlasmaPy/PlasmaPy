"""
Decorator to convert units of functions in /physics methods
"""
__all__ = ["angular_freq_to_hz"]

import functools
import inspect
from astropy import units as u
from plasmapy.utils.decorators import preserve_signature


def angular_freq_to_hz(fn):
    # raise exception if fn uses the 'to_hz' kwarg
    sig = inspect.signature(fn)
    if 'to_hz' in sig.parameters:
        raise raise ValueError(f"Wrapped function '{fn.__name__}'
        can not use keyword 'to_hz'.  Keyword reserved for decorator functionality.")

    # make new signature for fn
    new_params = sig.parameters.copy()
    new_params['to_hz'] = inspect.Parameter('to_hz', inspect.Parameter.POSITIONAL_OR_KEYWORD,
                                            default=False)
    new_sig = inspect.Signature(parameters=new_params.values(),
                                return_annotation=sig.return_annotation)
    fn.__signature__ = new_sig

    # add 'to_hz' to fn docstring
    fn.__doc__ = """
    Other Parameters
    ----------------
    to_hz: bool
        Set `True` to to convert function output from angular frequency to Hz
    """
    fn.__doc__ += """
    Additional parameters
    ------------
    to_hz: bool
        Set `True` to to convert function output from angular frequency to Hz
    """

    def wrapper(*args, to_hz=False, **kwargs):
        _result = fn(*args, **kwargs)
        if to_hz:
            return _result.to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])
        return _result
    return wrapper
