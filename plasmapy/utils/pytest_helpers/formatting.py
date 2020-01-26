import inspect
from astropy import units as u
from typing import Callable, Any, Dict, Optional, Union, Tuple

__all__ = [
    "call_string",
    "class_attribute_call_string",
    "class_method_call_string",
]

# TODO: Choose a more appropriate filename...perhaps errmsg_formatting.py?


def _format_quantity(arg) -> str:
    """
    Transform a `~astropy.units.Quantity` into a format as would appear
    in a function call.
    """

    formatted = f"{arg.value}"
    for base, power in zip(arg.unit.bases, arg.unit.powers):
        if power == -1:
            formatted += f"/u.{base}"
        elif power == 1:
            formatted += f"*u.{base}"
        else:
            formatted += f"*u.{base}**{power}"
    return formatted


def _format_arg(arg) -> str:
    """Transform an argument into a format as would appear in a function call."""

    if hasattr(arg, "__name__"):
        return arg.__name__
    elif isinstance(arg, u.Quantity):
        return _format_quantity(arg)
    else:
        return repr(arg)


def _format_kw(keyword) -> str:
    """Transform a keyword into a format as would appear in a function call."""

    if isinstance(keyword, str):
        return str(keyword)
    elif hasattr(keyword, "__name__"):
        return keyword.__name__
    else:
        return repr(keyword)


def _format_args_and_kwargs(args: Any = tuple(), kwargs: Dict = {}) -> str:
    """
    Take positional and keyword arguments, and format them into a
    string as they would appear in a function call.
    """

    args = args if isinstance(args, tuple) else (args,)

    args_and_kwargs = ""

    for arg in args:
        args_and_kwargs += f"{_format_arg(arg)}, "

    for kwarg in kwargs:
        args_and_kwargs += f"{_format_kw(kwarg)}={_format_arg(kwargs[kwarg])}, "

    if args_and_kwargs[-2:] == ", ":
        args_and_kwargs = args_and_kwargs[:-2]

    return args_and_kwargs


def _exc_str(ex: Exception) -> str:
    """
    Return a string with an indefinite article and the name of
    exception ``ex``.
    """

    exception_name = ex.__name__
    use_an = exception_name[0] in "aeiouAEIOU" and exception_name[0:4] != "User"
    indefinite_article = "an" if use_an else "a"
    return f"{indefinite_article} {exception_name}"


def _get_object_name(obj: Any, showmodule=False) -> str:
    """
    Return the name of an `object`.  If the `object` has a "__name__"
    attribute and ``showmodule`` is `True`, then prepend the module
    name if not in `builtins`.
    """

    obj_name = obj.__name__ if hasattr(obj, "__name__") else repr(obj)

    if hasattr(obj, "__name__") and showmodule is True:
        module_name = inspect.getmodule(obj).__name__
        if module_name != "builtins":
            obj_name = f"{module_name}.{obj_name}"

    return obj_name


def _string_together_warnings_for_printing(warning_types, warning_messages):
    """
    Take a list of warning types with a list of corresponding warning
    messages, and create a string that prints out each warning type
    followed by the corresponding message, separated by a full line.
    """

    warnings_with_messages = [
        _get_object_name(warning, showmodule=False) + ": " + message
        for warning, message in zip(warning_types, warning_messages)
    ]

    return "\n\n".join(warnings_with_messages)

    # TODO: Add tests!


def call_string(f: Callable, args: Any = tuple(), kwargs: Dict = {}) -> str:
    """
    Return a string with the equivalent call of a function or class.

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
    """

    args_and_kwargs = _format_args_and_kwargs(args, kwargs)
    return f"{f.__name__}({args_and_kwargs})"


def class_attribute_call_string(
    cls,
    attr: str,
    cls_args: Optional[Union[Tuple, Any]] = None,
    cls_kwargs: Optional[Dict[str, Any]] = None,
) -> str:
    """
    Return a string to represent accessing a class attribute.

    Parameters
    ----------
    cls : `class`
        The class to be used in the string representation

    attr: `str`
        The name of the class attribute

    cls_args : `tuple`, `list`, or any; optional
        A `tuple` or `list` containing the positional arguments to be
        used during instantiation of `cls`, or any other `object` if
        there is only one positional argument

    cls_kwargs: `dict`, optional
        A `dict` containing the keyword arguments to be used during
        instantiation of `cls`

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
    >>> class_attribute_call_string(SampleClass, 'attribute', cls_args, cls_kwargs)
    'SampleClass(1, kwarg1=2).attribute'
    """

    cls_args = tuple() if cls_args is None else cls_args
    cls_kwargs = dict() if cls_kwargs is None else cls_kwargs
    return f"{call_string(cls, cls_args, cls_kwargs)}.{attr}"


def class_method_call_string(
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
        The name of the class method

    cls_args : `tuple`, `list`, or any; optional
        A `tuple` or `list` containing the positional arguments to be
        used during instantiation of `cls`, or any other `object` if
        there is only one positional argument

    cls_kwargs: `dict`, optional
        A `dict` containing the keyword arguments to be used during
        instantiation of `cls`

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
    >>> class_method_call_string(SampleClass, 'method', c_args, c_kwargs, m_args, m_kwargs)
    'SampleClass(1, cls_kwarg=2).method(3, method_kwarg=4)'
    """

    class_call_string = f"{call_string(cls, cls_args, cls_kwargs)}"
    method_args_and_kwargs = _format_args_and_kwargs(method_args, method_kwargs)
    return f"{class_call_string}.{method}({method_args_and_kwargs})"
