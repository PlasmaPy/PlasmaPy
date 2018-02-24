""""""
import pytest
from typing import Callable, Dict, Tuple, Any
import inspect
import astropy.units as u


def call_string(f: Callable, args: Any=tuple(), kwargs: Dict={}) -> str:
    """Return a string with the equivalent call of a function."""

    def format_arg(arg):
        return arg.__name__ if hasattr(arg, '__name__') else repr(arg)

    def format_kw(keyword):
        if isinstance(keyword, str):
            return str(keyword)
        elif hasattr(keyword, '__name__'):
            return keyword.__name__
        else:
            return repr(keyword)

    if not isinstance(args, tuple):
        args = (args,)

    args_and_kwargs = ""

    for arg in args:
        args_and_kwargs += f"{format_arg(arg)}, "

    for kwarg in kwargs:
        args_and_kwargs += f"{format_kw(kwarg)}={format_arg(kwargs[kwarg])}, "

    if args_and_kwargs[-2:] == ", ":
        args_and_kwargs = args_and_kwargs[:-2]

    return f"{f.__name__}({args_and_kwargs})"


def run_test(f: Callable,
             args: Any = (),
             kwargs: Dict = {},
             expected: Any = None):
    """
    Test that a function or class returns the expected result, raises
    the expected exception, or issues an expected warning for the
    supplied positional and keyword arguments.

    Parameters
    ----------
    f: `function` or `class`
        The callable to be tested.

    args: `tuple` or `object`
        The positional arguments to `f`.

    kwargs: `dict`
        The keyword arguments to `f`.

    expected: `object`
        The expected result, exception, or warning from `f(*args, **kwargs)`.
        This may also be a tuple of length two

    """

    if not isinstance(args, tuple):
        args = (args,)

    call_string = _function_call_string(f, args, kwargs)

    if inspect.isclass(expected) and issubclass(expected, Exception):

        article = 'a'

        missing_exception_errmsg = (
            f"When testing {f.__module__}, the command\n\n"
            f"    {call_string}\n\n"
            f"did not raise {article} {expected.__name__} as expected. "
        )

        with pytest.raises(expected, message=(missing_exception_errmsg)):
            f(*args, **kwargs)
