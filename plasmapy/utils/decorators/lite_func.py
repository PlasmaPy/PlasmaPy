"""
Module for defining functionality that marks and handle Lite-Function
creation.
"""
__all__ = ["bind_lite_func", "rust_sanitize"]

import astropy.units as u
import functools
import inspect
import numpy as np

from numba.extending import is_jitted
from typing import Callable


class _LiteFuncDict(dict):
    """
    Dictionary of Lite-Function functionality bound to the parent
    function.  The dictionary key is a string indicating the name
    the functionality was bound with and the dictionary value is a
    string containing the fully qualified path of the original
    functionality.
    """

    # This is only to give __bound_lite_func__ a docstring.


def bind_lite_func(lite_func, attrs: dict[str, Callable] = None):
    """
    Decorator to bind a lightweight "lite" version of a formulary
    function to the full formulary function, as well as any supporting
    attributes.

    Parameters
    ----------
    lite_func: Callable
        The lightweight function to be bound as the ``lite`` attribute
        to the function being decorated.

    attrs: Dict[str, Callable]
        A dictionary where the key is a string defining the bound name
        and the associated value is the functionality to be bound.

    Examples
    --------
    .. code-block:: python

        def foo_lite(x):
            return x


        def bar():
            print("Supporting function.")


        @bind_lite_func(
            foo_lite,
            attrs=[
                ("bar", bar),
            ],
        )
        def foo(x):
            if not isinstance(x, float):
                raise TypeError("Argument x can only be a float.")
            return x

    .. code-block:: pycon

        >>> foo(5)  # doctest: +SKIP
        5
        >>> foo.lite(5)  # doctest: +SKIP
        5
        >>> foo.bar()  # doctest: +SKIP
        Supporting function.

    Notes
    -----
    In addition to binding the functionality defined by the inputs, a
    ``__bound_lite_func__`` dunder is bound.  This dunder is a
    dictionary where a key is a string representing the bound name of
    the bound functionality and the associated value is a string
    representing the fully qualified path of the original bound
    functionality.
    """
    if attrs is None:
        attrs = {}
    elif not isinstance(attrs, dict):
        raise TypeError(
            f"Argument 'attrs' is a type '{type(attrs)}', expected a dictionary."
        )
    elif "lite" in attrs:
        raise ValueError(
            "Argument 'attr' can NOT define key 'lite', this is reserved for"
            " the 'lite_func' argument."
        )

    if inspect.isbuiltin(lite_func) or not (
        is_jitted(lite_func) or inspect.isfunction(lite_func)
    ):
        raise ValueError("The given lite-function is not a user-defined function.")

    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            return f(*args, **kwargs)

        __bound_lite_func__ = _LiteFuncDict()

        attrs["lite"] = lite_func
        for bound_name, attr in attrs.items():
            # only allow functions or jitted functions
            if not (inspect.isfunction(attr) or is_jitted(attr)):
                raise ValueError(
                    f"Can not bind obj '{attr}' to function '{wrapper.__name__}'."
                    f"  Only functions are allowed to be bound. Skipping."
                )

            # build origin name
            if hasattr(attr, "__module__"):
                modname = attr.__module__
            else:  # coverage: ignore
                # assume attr is defined in the module the function being
                # decorated in
                modname = wrapper.__module__
            origin = f"{modname}.{attr.__name__}"
            __bound_lite_func__[bound_name] = origin

            # bind
            setattr(wrapper, bound_name, attr)

        wrapper.__bound_lite_func__ = __bound_lite_func__

        return wrapper

    return decorator


def rust_sanitize(ndarray_type_dictionary: dict[str, object]):
    def decorator(target_function):
        signature = inspect.signature(target_function)

        for parameter in ndarray_type_dictionary:
            if parameter not in signature.parameters:
                raise TypeError(
                    f"Sanitize decorator got unexpected keyword argument '{parameter}'"
                )

        @functools.wraps(target_function)
        def wrapper(*args, **kwargs):
            bound_args = signature.bind(*args, **kwargs)
            bound_args.apply_defaults()

            sanitized_arguments = {}

            for arg_name, arg_value in bound_args.arguments.items():
                new_value = arg_value
                if isinstance(arg_value, u.Quantity):
                    new_value = arg_value.value

                if arg_name in ndarray_type_dictionary:
                    array_type = ndarray_type_dictionary[arg_name]

                    if (
                        isinstance(arg_value, np.ndarray)
                        and arg_value.dtype is not array_type
                    ):
                        new_value = arg_value.astype(array_type)
                    elif not isinstance(arg_value, np.ndarray):
                        new_value = np.asarray(arg_value, dtype=array_type)

                sanitized_arguments[arg_name] = np.atleast_1d(new_value)

            result = target_function(**sanitized_arguments)

            return result if result.shape else result[()]

        return wrapper

    return decorator
