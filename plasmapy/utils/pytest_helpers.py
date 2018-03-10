"""Utilities for """
import pytest
from typing import Callable, Dict, Tuple, Any, Union, Optional
import inspect
import warnings
from inspect import isclass
import astropy.units as u
from collections import defaultdict


class RunTestError(Exception):
    pass


class MissingExceptionError(RunTestError):
    pass


class UnexpectedExceptionError(RunTestError):
    pass


class MissingWarningError(RunTestError):
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


def run_test(func: Callable, args: Any = (), kwargs: Dict = {}, expected_outcome: Any = None):
    """
    Test that a function or class returns the expected result, raises
    the expected exception, or issues an expected warning for the
    supplied positional and keyword arguments.

    Parameters
    ----------
    func: callable
        The callable to be tested.

    args: tuple or object
        The positional arguments to `func`.

    kwargs: dict
        The keyword arguments to `func`.

    expected_outcome: object
        The expected result, exception, or warning from
        `func(*args, **kwargs)`. This may also be a tuple of length two
        that contains the expected result as the first item and the
        expected warning as the second item.

    Examples
    --------
    >>> from warnings import warn
    >>> run_test(lambda: 42, tuple(), dict(), 42)
    >>> run_test(lambda: raise ValueError, tuple(), dict(), ValueError)
    >>> run_test(lambda: warn("", UserWarning), tuple(), dict(), UserWarning)

    """

    def exc_str(ex: Exception) -> str:
        """
        Return a string with an indefinite article and the name of
        exception `ex`.
        """
        article = 'an' if ex.__name__[0] in 'aeiouAEIOU' else 'a'
        return f"{article} {ex.__name__}"

    if not isinstance(args, tuple):
        args = (args,)

    # We will want to include a string that can reproduce the test for
    # all error messages.

    call_str = call_string(func, args, kwargs)

    # print(f"call_str = {call_str}")

    # There are many possibilities for expected outcomes that we must
    # keep track of, including exceptions being raised and warnings
    # being issued.

    expected = defaultdict(lambda: None)

    if inspect.isclass(expected_outcome):
        subclass_of_Exception = issubclass(expected_outcome, Exception)
        subclass_of_Warning = issubclass(expected_outcome, Warning)

        if subclass_of_Warning:
            expected['warning'] = expected_outcome
        if subclass_of_Exception and not subclass_of_Warning:
            expected['exception'] = expected_outcome

    # If a warning is issued, then there may also be an expected result.

    if isinstance(expected_outcome, tuple):
        if len(expected_outcome) > 2 or not inspect.isclass(expected_outcome[1]):
            raise ValueError
        if issubclass(expected_outcome[1], Warning):
            expected['result'] = expected_outcome[0]
            expected['warning'] = expected_outcome[1]

        else:
            raise ValueError

    if expected['exception'] is None and expected['warning'] is None:
        expected['result'] = expected_outcome

    expected['type'] = type(expected['result']) if expected['result'] is not None else None

    for key in expected.keys():
        print(f"expected[{key}] = {expected[key]}")

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
            result = func(*args, **kwargs)
        except exception:
            pass
        except Exception as exc:
            unexpected_exception = exc.__reduce__()[0]
            raise UnexpectedExceptionError(f"\nRunning the command:\n\n"
                                           f"  {call_str}\n\n"
                                           f"did not raise {exc_str(exception)} as expected, "
                                           f"but instead raised {exc_str(unexpected_exception)}.\n") from exc
        else:
            raise MissingExceptionError(f"\nRunning the command:\n\n"
                                        f"  {call_str}\n\n"
                                        f"did not raise {exc_str(exception)} as "
                                        f"expected, but instead returned the value:\n\n"
                                        f"  {result}\n")
    else:

        try:
            with pytest.warns(expected['warning']):
                result = func(*args, **kwargs)
        except pytest.raises.Exception as exc:
            raise MissingWarningError(f"\nRunning the command:\n"
                                      f"  {call_str}\n"
                                      f"should issue {exc_str(expected['warning'])}, "
                                      f"but instead returned the value:\n"
                                      f"  {result}\n")
        except Exception as exc:

            raise UnexpectedExceptionError(f"Running the command {call_str} "
                                           f"unexpectedly raised {exc_str(exc.__reduce__()[0])} "
                                           f"instead of returning the expected value of {expected['result']}.")

        if expected['result'] is not None:
            assert result == expected['result'], (
                f"Running the command {call_str} yielded a result of "
                f"{result} instead of returning the expected value of "
                f"{expected['result']}.")
