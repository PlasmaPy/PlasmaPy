"""Tools for formatting strings, including for error messages."""

__all__ = ["call_string", "attribute_call_string", "method_call_string"]

import inspect
import numbers
import numpy as np

from astropy import units as u
from typing import Any, Callable, Dict, List, Optional, Tuple, Union


def _code_repr_of_ndarray(array: np.ndarray, max_items=np.inf) -> str:
    """
    Transform a `~numpy.ndarray` into a format as would appear in a
    function call, after having done ``import numpy as np``.

    If ``max_items`` is less than ``array.size``, then include only the
    first ``max_items`` elements of the array in the resulting string
    and replace the remaining items with ``"..."``.
    """

    def remove_excess_spaces(string: str) -> str:
        string = " ".join(string.split())
        string = string.replace(" ,", ",")
        string = string.replace("[ ", "[")
        return string

    def put_np_before_infs_and_nans(string: str) -> str:
        string = string.replace("inf", "np.inf")
        string = string.replace("nan", "np.nan")
        return string

    def replace_excess_items_with_ellipsis(s, max_items):
        substrings_between_commas = s.split(",")
        to_comma_before_ellipsis = ",".join(substrings_between_commas[0:max_items])
        substring_after_last_comma = substrings_between_commas[-1]
        closing_brackets = "]" * substring_after_last_comma.count("]")
        closing = f", ...{closing_brackets})"
        return to_comma_before_ellipsis + closing

    if not isinstance(array, np.ndarray):
        raise TypeError("Expecting an ndarray.")

    s = np.array_repr(array, max_line_width=np.inf, suppress_small=False)
    s = remove_excess_spaces(s)
    s = put_np_before_infs_and_nans(s)

    if array.size > max_items:
        s = replace_excess_items_with_ellipsis(s, max_items)

    return f"np.{s}"


def _code_repr_of_quantity(arg: u.Quantity) -> str:
    """
    Transform a `~astropy.units.Quantity` into a format as would appear
    in a function call.
    """
    if not isinstance(arg, u.Quantity):
        raise TypeError("Expecting a Quantity.")

    formatted = f"{repr(arg.value)}"

    if isinstance(arg.value, np.ndarray):
        formatted = "np." + formatted

    if arg.unit == u.dimensionless_unscaled:
        formatted += "*u.dimensionless_unscaled"
    else:
        for base, power in zip(arg.unit.bases, arg.unit.powers):
            if power == -1:
                formatted += f"/u.{base}"
            elif power == 1:
                formatted += f"*u.{base}"
            else:
                formatted += f"*u.{base}**{power}"

    return formatted


def _code_repr_of_arg(arg) -> str:
    """Transform an argument into a format as would appear in a function call."""
    if hasattr(arg, "__name__"):
        return arg.__name__

    function_for_type = {
        u.Quantity: _code_repr_of_quantity,
        np.ndarray: _code_repr_of_ndarray,
    }

    for arg_type in function_for_type.keys():
        if isinstance(arg, arg_type):
            code_repr_func = function_for_type[arg_type]
            return code_repr_func(arg)

    return repr(arg)


def _code_repr_of_args_and_kwargs(args: Any = tuple(), kwargs: Dict = {}) -> str:
    """
    Take positional and keyword arguments, and format them into a
    string as they would appear in a function call.
    """
    args_collection = args if isinstance(args, tuple) else (args,)

    args_and_kwargs = ""

    for arg in args_collection:
        args_and_kwargs += f"{_code_repr_of_arg(arg)}, "

    for kwarg in kwargs.keys():
        args_and_kwargs += f"{kwarg}={_code_repr_of_arg(kwargs[kwarg])}, "

    if args_and_kwargs[-2:] == ", ":
        args_and_kwargs = args_and_kwargs[:-2]

    return args_and_kwargs


def _name_with_article(ex: Exception) -> str:
    """
    Return a string with an indefinite article and the name of
    exception ``ex``.
    """
    exception_name = ex.__name__
    use_an = exception_name[0] in "aeiouAEIOU" and exception_name[0:4] != "User"
    indefinite_article = "an" if use_an else "a"
    return f"{indefinite_article} {exception_name}"


def _object_name(obj: Any, showmodule=False) -> str:
    """
    Return the name of an `object`.  If the `object` has a `__name__`
    attribute and ``showmodule`` is `True`, then prepend the module
    name if not in `builtins`.
    """
    obj_name = obj.__name__ if hasattr(obj, "__name__") else repr(obj)

    if hasattr(obj, "__name__") and showmodule is True:
        module_name = inspect.getmodule(obj).__name__  # ADD TEST
        if module_name != "builtins":
            obj_name = f"{module_name}.{obj_name}"  # ADD TEST

    return obj_name


def _string_together_warnings_for_printing(
    warning_types: List[Warning], warning_messages: List[str]
):
    """
    Take a list of warning types with a list of corresponding warning
    messages, and create a string that prints out each warning type
    followed by the corresponding message, separated by a full line.
    """
    warnings_with_messages = [
        _object_name(warning, showmodule=False) + ": " + message
        for warning, message in zip(warning_types, warning_messages)
    ]
    return "\n\n".join(warnings_with_messages)


def call_string(f: Callable, args: Any = tuple(), kwargs: Dict = {}) -> str:
    """
    Return a string with the equivalent call of a callable object such
    as a function or class.

    Parameters
    ----------
    f : callable
        The callable object to be used in the string representation

    args : `tuple`, `list`, or any; optional
        A `tuple` or `list` containing positional arguments, or any
        other `object` if there is only one positional argument

    kwargs : `dict`, optional
        A `dict` containing keyword arguments

    Examples
    --------
    >>> call_string(int, 3.14159)
    'int(3.14159)'
    >>> call_string(int, args=(9.2,), kwargs={'base': 2})
    'int(9.2, base=2)'

    See Also
    --------
    attribute_call_string
    method_call_string
    """
    args_and_kwargs = _code_repr_of_args_and_kwargs(args, kwargs)
    return f"{f.__name__}({args_and_kwargs})"


def attribute_call_string(
    cls,
    attr: str,
    cls_args: Optional[Union[Tuple, Any]] = None,
    cls_kwargs: Optional[Dict[str, Any]] = None,
) -> str:
    """
    Return a string to represent accessing an attribute of an object,
    such as an attribute of an instance in a class.

    Parameters
    ----------
    cls : `class`
        The class to be used in the string representation

    attr: `str`
        The name of the attribute of class ``cls``

    cls_args : `tuple`, `list`, or any; optional
        A `tuple` or `list` containing positional arguments, or any
        other type of `object` if there is only one positional argument,
        to be used during instantiation of ``cls``

    cls_kwargs: `dict`, optional
        A `dict` containing the keyword arguments to be used during
        instantiation of ``cls``

    Examples
    --------
    >>> class SampleClass:
    ...     def __init__(self, arg1, kwarg1=None):
    ...         pass
    ...     @property
    ...     def attribute(self):
    ...         return 42
    >>> cls_args = (1,)
    >>> cls_kwargs = {'kwarg1': 2}
    >>> attribute_call_string(SampleClass, 'attribute', cls_args, cls_kwargs)
    'SampleClass(1, kwarg1=2).attribute'
    """
    cls_args = tuple() if cls_args is None else cls_args
    cls_kwargs = dict() if cls_kwargs is None else cls_kwargs
    return f"{call_string(cls, cls_args, cls_kwargs)}.{attr}"


def method_call_string(
    cls,
    method: str,
    cls_args: Any = tuple(),
    cls_kwargs: Dict = {},
    method_args: Any = tuple(),
    method_kwargs: Dict = {},
) -> str:
    """
    Return a string to represent calling a class method.

    Parameters
    ----------
    cls : `class`
        The class to be used in the string representation

    method: `str`
        The name of the method in class ``cls``

    cls_args : `tuple`, `list`, or any; optional
        A `tuple` or `list` containing positional arguments, or any
        other type of `object` if there is only one positional argument,
        to be used during instantiation of ``cls``


    cls_kwargs: `dict`, optional
        A `dict` containing the keyword arguments to be used during
        instantiation of ``cls``

    method_args : `tuple`, `list`, or any; optional
        A `tuple` or `list` containing the positional arguments to be
        used in the method call, or any other `object` if there is only
        one positional argument

    method_kwargs: `dict`, optional
        A `dict` containing the keyword arguments to be used during
        the method call

    Examples
    --------
    >>> class SampleClass:
    ...     def __init__(self, cls_arg, cls_kwarg=None):
    ...         pass
    ...     def method(self, method_arg, method_kwarg=None):
    ...         return 42
    >>> c_args = (1,)
    >>> c_kwargs = {'cls_kwarg': 2}
    >>> m_args = 3
    >>> m_kwargs = {'method_kwarg': 4}
    >>> method_call_string(SampleClass, 'method', c_args, c_kwargs, m_args, m_kwargs)
    'SampleClass(1, cls_kwarg=2).method(3, method_kwarg=4)'

    See Also
    --------
    call_string
    attribute_call_string
    """
    class_call_string = f"{call_string(cls, cls_args, cls_kwargs)}"
    method_args_and_kwargs = _code_repr_of_args_and_kwargs(method_args, method_kwargs)
    return f"{class_call_string}.{method}({method_args_and_kwargs})"
