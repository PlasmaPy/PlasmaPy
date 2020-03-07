"""
Miscellaneous decorators for various package uses.
"""
__all__ = ["preserve_signature"]

import functools
import inspect


def preserve_signature(f):
    """
    A decorator for decorators, which preserves the signature of the function
    being wrapped. This preservation allows IDE function parameter hints to work
    on the wrapped function. To do this, the `__signature__` dunder is defined, or
    inherited, from the function being wrapped to the resulting wrapped function.

    Parameters
    ----------
    f: callable
        The function being wrapped.

    Returns
    -------
    callable
        Wrapped version of the function.

    Examples
    --------
    >>> def a_decorator(f):
    ...     @preserve_signature
    ...     @functools.wraps(f)
    ...     def wrapper(*args, **kwargs):
    ...         return wrapper(*args, **kwargs)
    ...
    ...     return wrapper
    """
    # add '__signature__' to methods that are copied from
    # f onto wrapper
    assigned = list(functools.WRAPPER_ASSIGNMENTS)
    assigned.append("__signature__")

    @functools.wraps(f, assigned=assigned)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    # add '__signature__' if it does not exist
    # - this will preserve parameter hints in IDE's
    if not hasattr(wrapper, "__signature__"):
        wrapper.__signature__ = inspect.signature(f)

    return wrapper
