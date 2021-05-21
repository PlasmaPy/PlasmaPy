"""
Module for defining functionality that marks and handle Lite-Function creation.
"""
__all__ = ["mark_has_lite_func"]

import functools
import inspect

from collections import namedtuple
from typing import List, Tuple
from warnings import warn


LiteFuncTupleEntry = namedtuple("LiteFuncTupleEntry", ("name", "origin"))

_litefunc_registry = {}


def mark_has_lite_func(flite, attrs: List[Tuple[str, str]] = None, scope=None):
    """
    Decorator to bind lightweight "lite" versions of formulary functions to the full
    formulary function.
    """
    if attrs is None:
        attrs = []
    elif scope is None:
        raise ValueError(
            f"If 'attrs' are being bound to the decorate function, then 'scope'"
            f" needs to be defined and contain the objects to be bound.  "
            f"Setting 'scope=globals()' will typically suffice."
        )

    if not inspect.isfunction(flite) or inspect.isbuiltin(flite):
        raise ValueError(f"The lite-function passed is not a user-defined function.")

    def decorator(f):
        parent_qualname = f"{f.__module__}.{f.__name__}"

        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            return f(*args, **kwargs)

        __has_litefuncs__ = []
        attrs.append(("lite", ""))
        for bound_name, attr_name in attrs:
            # get object to bind
            if bound_name == "lite":
                attr = flite
            else:
                attr = scope[attr_name]

            # skip objects that are not allowed
            _warn = None
            if inspect.isclass(attr):
                _warn = "class"
            elif inspect.ismodule(attr):
                _warn = "module"
            elif inspect.ismethod(attr):
                _warn = "method"
            if _warn is not None:
                warn(
                    f"Can not bind {_warn} '{attr_name}' to function"
                    f" '{wrapper.__name__}'.  Only functions and instances are "
                    f"allowed. Skipping."
                )
                continue

            # build origin name
            if hasattr(attr, "__module__"):
                modname = attr.__module__
            else:
                modname = wrapper.__module__

            if hasattr(attr, "__name__"):
                shortname = attr.__name__
            else:
                shortname = attr_name

            origin = f"{modname}.{shortname}"
            __has_litefuncs__.append(
                LiteFuncTupleEntry(name=bound_name, origin=origin)
            )

            # bind
            setattr(wrapper, bound_name, attr)

            reg_entry = {
                f"{parent_qualname}.{bound_name}" : {
                    "is_parent": False,
                    "is_child": True,
                    "parent": parent_qualname,
                    "shortname": bound_name,
                    "origin": origin,
                },
            }
            _litefunc_registry.update(reg_entry)

        if len(__has_litefuncs__) == 0:
            raise ValueError(
                f"Lite-function marking of '{wrapper.__name__}' resulting in NO"
                f" attributes being bound."
            )

        setattr(wrapper, "__has_litefunc__", __has_litefuncs__)

        reg_entry = {
            f"{parent_qualname}": {
                "is_parent": True,
                "is_child": False,
                "parent": None,
                "shortname": parent_qualname.split(".")[-1],
                "origin": parent_qualname,
            },
        }
        _litefunc_registry.update(reg_entry)

        return wrapper

    return decorator
