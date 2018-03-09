"""Utilities for """
import pytest
from typing import Callable, Dict, Tuple, Any, Union, Optional
import inspect
from inspect import isclass
import astropy.units as u
from collections import defaultdict


class TestError(Exception):
    pass


class MissingExceptionError(TestError):
    pass


class UnexpectedExceptionError(TestError):
    pass


class MissingWarningError(TestError):
    pass


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
             expected_outcome: Any = None):
    """
    Test that a function or class returns the expected result, raises
    the expected exception, or issues an expected warning for the
    supplied positional and keyword arguments.

    Parameters
    ----------
    f: callable
        The callable to be tested.

    args: tuple or object
        The positional arguments to `f`.

    kwargs: dict
        The keyword arguments to `f`.

    expected_outcome: object
        The expected result, exception, or warning from
        `f(*args, **kwargs)`. This may also be a tuple of length two
        that contains the expected result as the first item and the
        expected warning as the second item.

    """

    def exc_str(ex: Exception) -> str:
        """
        Return a string with an indefinite article and the name of
        exception `ex`.
        """
        article = 'an' if ex.__name__[0] in 'aeiouAEIOU' else 'a'
        return f"{article} {ex.__name__[0]}"

    if not isinstance(args, tuple):
        args = (args,)

    # We will want to include a string that can reproduce the test for
    # all error messages.

    call_str = call_string(f, args, kwargs)

    # There are many possibilities for expected outcomes that we must
    # keep track of, including exceptions being raised and warnings
    # being issued.

    expected = defaultdict(lambda: None)

    if inspect.isclass(expected):
        if issubclass(expected, Exception):
            expected['exception'] = expected_outcome
        elif issubclass(expected, Warning):
            expected['warning'] = expected_outcome
        elif issubclass(expected, u.Unit):
            expected['unit'] = expected_outcome

    # If a warning is issued, then there may also be an expected result.

    if isinstance(expected_outcome, tuple):
        tuple_contains_warning = inspect.isclass(expected_outcome[1])
        expected['warning'] = expected_outcome[1] if tuple_contains_warning else None
        expected['result'] = expected_outcome[0]

    if not expected['exception'] and not expected['result']:
        expected['result'] = expected_outcome

    expected['type'] = type(expected['result'])

    # First we go through all of the possibilities for when an exception
    # is expected to be raised.  If no exception is raised, then we want
    # an error message that includes the result.  If the wrong exception
    # is raised, then we want an error message that includes that
    # exception.  An alternative would be to use `with pytest.raises()`
    # but this makes it easier to break down what the error messages
    # should be.

    if expected['exception']:
        exception = expected['exception']

        try:
            result = f(*args, **kwargs)
        except exception:
            pass
        except Exception as exc:
            unexpected_exception = exc.__reduce__()[0]
            raise UnexpectedExceptionError(
                f"Running the command:\n\n"
                f"  {call_str}\n\n"
                f"did not raise {exc_str(exception)} as expected, "
                f"but instead raised {exc_str(unexpected_exception)}.\n"
            ) from exc
        else:
            raise MissingExceptionError(
                f"Running the command:\n\n"
                f"  {call_str}\n\n"
                f"did not raise a {exc_str(exception)} as "
                f"expected, but instead returned the value:\n\n"
                f"  {result}\n")
