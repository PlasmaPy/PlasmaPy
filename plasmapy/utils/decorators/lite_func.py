"""
Module for defining functionality that marks and handle Lite-Function creation.
"""
__all__ = ["bind_lite_func"]

import functools
import inspect

from typing import Callable, List, Tuple
from warnings import warn


_litefunc_registry = {}


class LiteFuncList(list):
    """
    A list containing which attributions have been bound as part of the
    "Lite-Function" functionality.  Each entry is a 2-element tuple where the
    first element is the bound name and the second element is the fully
    qualified name of the original object.
    """


def bind_lite_func(flite, attrs: List[Tuple[str, Callable]] = None):
    """
    Decorator to bind lightweight "lite" versions of formulary functions to the full
    formulary function.
    """
    if attrs is None:
        attrs = []

    if not inspect.isfunction(flite) or inspect.isbuiltin(flite):
        raise ValueError(f"The lite-function passed is not a user-defined function.")

    def decorator(f):
        parent_qualname = f"{f.__module__}.{f.__name__}"

        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            return f(*args, **kwargs)

        __has_lite_func__ = LiteFuncList()
        attrs.append(("lite", flite))
        for bound_name, attr in attrs:

            # skip objects that are not allowed
            # - only allow functions
            if not inspect.isfunction(attr):
                warn(
                    f"Can not bind obj '{attr}' to function '{wrapper.__name__}'."
                    f"  Only functions are allowed to be bound. Skipping."
                )
                continue

            # build origin name
            if hasattr(attr, "__module__"):
                modname = attr.__module__
            else:
                modname = wrapper.__module___
            origin = f"{modname}.{attr.__name__}"
            __has_lite_func__.append((bound_name, origin))

            # bind
            setattr(wrapper, bound_name, attr)

            reg_entry = {
                f"{parent_qualname}.{bound_name}" : {
                    "is_parent": False,
                    "parent": parent_qualname,
                    "shortname": bound_name,
                    "origin": origin,
                },
            }
            _litefunc_registry.update(reg_entry)

        if len(__has_lite_func__) == 0:
            raise ValueError(
                f"Lite-function marking of '{wrapper.__name__}' resulting in NO"
                f" attributes being bound."
            )

        setattr(wrapper, "__has_lite_func__", __has_lite_func__)

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
