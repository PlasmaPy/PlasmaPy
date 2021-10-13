"""
Module for defining functionality that marks and handle Lite-Function
creation.
"""
__all__ = ["bind_lite_func"]

import functools
import inspect

from typing import Callable, List, Tuple
from warnings import warn

from plasmapy.utils.exceptions import PlasmaPyWarning

_litefunc_registry = {}


def bind_lite_func(lite_func, attrs: List[Tuple[str, Callable]] = None):
    """
    Decorator to bind a lightweight "lite" version of a formulary
    function to the full formulary function, as well as any supporting
    attributes.

    Parameters
    ----------
    lite_func: Callable
        The lightweight function to be bound as the ``lite`` attribute
        to the function being decorated.

    attrs: List[Tuple[str, Callable]]
        A list of 2-element tuples where the second element is a
        function to be bound and the first element is a string giving
        the name for binding.

    Examples
    --------

    .. code-block:: python

        def foo_lite(x)
            return x

        def bar():
            print("Supporting function.")

        @bind_lite_func(foo_lite, attrs=[("bar", bar),])
        def foo(x):
            if not isinstance(x, float):
                raise TypeError("Argument x can only be a float.")
            return x

        >>> foo(5)
        5
        >>> foo.lite(5)
        5
        >>> foo.bar()
        Supporting function.

    """
    if attrs is None:
        attrs = []

    if not inspect.isfunction(lite_func) or inspect.isbuiltin(lite_func):
        raise ValueError(f"The given lite-function is not a user-defined function.")

    def decorator(f):
        parent_qualname = f"{f.__module__}.{f.__name__}"

        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            return f(*args, **kwargs)

        __bound_lite_func__ = {}

        attrs.append(("lite", lite_func))
        for bound_name, attr in attrs:
            # skip objects that are not allowed
            # - only allow functions
            if not inspect.isfunction(attr):
                warn(
                    f"Can not bind obj '{attr}' to function '{wrapper.__name__}'."
                    f"  Only functions are allowed to be bound. Skipping.",
                    PlasmaPyWarning
                )
                continue

            # build origin name
            if hasattr(attr, "__module__"):
                modname = attr.__module__
            else:
                # assume attr is defined in the module the function being
                # decorated is defined in
                modname = wrapper.__module__
            origin = f"{modname}.{attr.__name__}"
            if bound_name in __bound_lite_func__:
                raise ValueError(
                    f"Can NOT bind the same attribute name '{bound_name}' more"
                    f" than once."
                )
            __bound_lite_func__[bound_name] = origin

            # bind
            setattr(wrapper, bound_name, attr)

            # add to lite function registry
            reg_entry = {
                f"{parent_qualname}.{bound_name}" : {
                    "is_parent": False,
                    "parent": parent_qualname,
                    "shortname": bound_name,
                    "origin": origin,
                },
            }
            _litefunc_registry.update(reg_entry)

        setattr(wrapper, "__bound_lite_func__", __bound_lite_func__)

        reg_entry = {
            f"{parent_qualname}": {
                "is_parent": True,
                "parent": None,
                "shortname": parent_qualname.split(".")[-1],
                "origin": parent_qualname,
            },
        }
        _litefunc_registry.update(reg_entry)

        return wrapper

    return decorator
