"""Test helper utilities."""

import pytest
from typing import Callable, Dict, Tuple, Any, Union, Optional
import inspect
import warnings
from inspect import isclass
import astropy.units as u
import astropy.constants as const
from collections import defaultdict
from textwrap import fill
import numpy as np

class RunTestError(Exception):
    """Base exception for test failures. Derived from `Exception`."""

    pass


class MissingExceptionError(RunTestError):
    """
    Exception for when an expected exception is not raised.  Derived
    from `~plasmapy.utils.RunTestError`.
    """

    pass


class UnexpectedExceptionError(RunTestError):
    """
    Exception for when an exception is expected, but a different
    exception is raised instead.  Derived from
    `~plasmapy.utils.RunTestError`.
    """

    pass


class MissingWarningError(RunTestError):
    """
    Exception for when a warning is expected to be issued, but isn't.
    Derived from `~plasmapy.utils.RunTestError`.
    """

    pass


class IncorrectResultError(RunTestError):
    """
    Exception for when the result is
    """

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


def run_test(
        func: Callable,
        args: Any = (),
        kwargs: Dict = {},
        expected_outcome: Any = None,
        rtol: float = 0.0,
        atol: float = 0.0,
        ):
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

    rtol : float

    atol : float

    Examples
    --------
    >>> from warnings import warn
    >>> run_test(lambda: 42, tuple(), dict(), 42)
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

    # First we go through all of the possibilities for when an exception
    # is expected to be raised.  If no exception is raised, then we want
    # an error message that includes the result.  If the wrong exception
    # is raised, then we want an error message that includes that
    # exception.  An alternative would be to use `with pytest.raises()`
    # but this makes it easier to break down what the error messages
    # should be.

    if expected['exception']:

        expected_exception = expected['exception']

        try:
            result = func(*args, **kwargs)
        except expected_exception:
            return None
        except Exception as exc_unexpected_exception:
            unexpected_exception = exc_unexpected_exception.__reduce__()[0]
            raise UnexpectedExceptionError(
                f"Running the command:\n\n"
                f"  {call_str}\n\n"
                f"did not raise {exc_str(expected_exception)} as expected, "
                f"but instead raised {exc_str(unexpected_exception)}."
            ) from exc_unexpected_exception
        else:
            raise MissingExceptionError(
                f"Running the command:\n\n"
                f"  {call_str}\n\n"
                f"did not raise {exc_str(expected_exception)} as "
                f"expected, but instead returned the value:\n\n"
                f"  {result}\n"
            )

    try:
        with pytest.warns(expected['warning']):
            result = func(*args, **kwargs)
    except pytest.raises.Exception as missing_warning:

        raise MissingWarningError(
            f"Running the command:\n\n"
            f"  {call_str}\n\n"
            f"should issue {exc_str(expected['warning'])}, "
            f"but instead returned:\n\n"
            f"  {result}\n"
        ) from missing_warning

    except Exception as exception_no_warning:

        raise UnexpectedExceptionError(
            f"Running the command:\n\n"
            f"  {call_str}\n\n"
            f"unexpectedly raised {exc_str(exception_no_warning.__reduce__()[0])} "
            f"instead of returning the expected value of:\n\n"
            f"  {expected['result']}\n") from exception_no_warning

    if expected['result'] is None:
        return None

    try:
        if result == expected['result']:
            return None
    except NotImplemented as exc_equality:
        raise RunTestError(
            f"The equality of {result} and {expected['result']} "
            f"cannot be evaluated.") from exc_equality

    if isinstance(expected['result'], u.UnitBase):

        if isinstance(result, u.UnitBase):
            if result != expected['result']:
                raise u.UnitsError(
                    f"Running the command:\n\n"
                    f"  {call_str}\n\n"
                    f"returned {repr(result)} instead of the expected "
                    f"value of {repr(expected['result'])}.")
            return None

        if not isinstance(result, (u.Quantity, const.Constant, const.EMConstant)):
            raise u.UnitsError(
                f"Running the command:\n\n"
                f"  {call_str}\n\n"
                f"returned a value of:\n\n"
                f"  {result}\n\n"
                f"instead of a quantity or constant with units of "
                f"{expected['result']}.")

        if result.unit != expected['result']:
            raise u.UnitsError(
                f"Running the command:\n\n"
                f"  {call_str}\n\n"
                f"returned the result:\n\n"
                f"  {result}\n\n"
                f"which has units of {result.unit} instead of the"
                f"expected units of {expected['result']}.")

        return None

    if isinstance(expected['result'], (u.Quantity, const.Constant, const.EMConstant)):
        if not result.unit == expected['result'].unit:
            raise u.UnitsError(
                f"Running the command:\n\n"
                f"  {call_str}\n\n"
                f"returned the result:\n\n"
                f"  {result}\n\n"
                f"which has different units than the expected result of:\n\n"
                f"  {expected['result']}\n")

        if np.allclose(result.value, expected['result'].value):
            return None

    if np.allclose(result, expected['result'], rtol=rtol, atol=atol):
        return None

    raise RunTestError(
        f"Running the command:\n\n"
        f"  {call_str}\n\n"
        f"returned the result:\n\n"
        f"  {result}\n\n"
        f"instead of the expected value of:\n\n"
        f"  {expected['result']}\n")
