"""Test helper utilities."""
import functools
import pytest
import inspect
from collections import defaultdict
from typing import Callable, Dict, Any
import numpy as np
import astropy.units as u
import astropy.constants as const
import colorama

# These colors/styles are used to highlight certain parts of the error
# messages in consistent ways.

_bold = colorama.Style.BRIGHT
_magenta = colorama.Fore.MAGENTA
_blue = colorama.Fore.BLUE
_cyan = colorama.Fore.CYAN
_red = colorama.Fore.RED

_exception_color = f"{_magenta}{_bold}"
_type_color = f"{_magenta}{_bold}"
_func_color = f"{_cyan}{_bold}"
_result_color = f"{_blue}{_bold}"
_message_color = f"{_red}{_bold}"


class RunTestError(Exception):
    """Base exception for test failures. Derived from `Exception`."""


class UnexpectedResultError(RunTestError):
    """
    Exception for when the actual result differs from the expected
    result.  Derived from `~plasmapy.utils.RunTestError`.
    """


class InconsistentTypeError(RunTestError):
    """
    Exception for when the type of the actual result differs from the
    type of the expected result.  Derived from
    `~plasmapy.utils.RunTestError`.
    """


class MissingExceptionError(RunTestError):
    """
    Exception for when an expected exception is not raised.  Derived
    from `~plasmapy.utils.RunTestError`.
    """


class UnexpectedExceptionError(RunTestError):
    """
    Exception for when an exception is expected, but a different
    exception is raised instead.  Derived from
    `~plasmapy.utils.RunTestError`.
    """


class MissingWarningError(RunTestError):
    """
    Exception for when a warning is expected to be issued, but isn't.
    Derived from `~plasmapy.utils.RunTestError`.
    """


class IncorrectResultError(RunTestError):
    """
    Exception for when the actual result differs from the expected
    result by more than the allowed tolerance.  Derived from
    `~plasmapy.utils.RunTestError`.
    """


def call_string(f: Callable,
                args: Any=tuple(),
                kwargs: Dict={},
                color="",
                return_color="",
                ) -> str:
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

    if color and not return_color:
        return_color = _message_color

    if not isinstance(args, tuple):
        args = (args,)

    args_and_kwargs = ""

    for arg in args:
        args_and_kwargs += f"{format_arg(arg)}, "

    for kwarg in kwargs:
        args_and_kwargs += f"{format_kw(kwarg)}={format_arg(kwargs[kwarg])}, "

    if args_and_kwargs[-2:] == ", ":
        args_and_kwargs = args_and_kwargs[:-2]

    return f"{color}{f.__name__}({args_and_kwargs}){return_color}"


def _exc_str(ex: Exception, color=_exception_color) -> str:
    """
    Return a string with an indefinite article and the name of
    exception `ex`.
    """
    if color is None:
        color = ""
        return_color = ""
    else:
        return_color = _message_color
    exception_name = ex.__name__
    use_an = exception_name[0] in 'aeiouAEIOU' and exception_name[0:4] != "User"
    article = 'an' if use_an else 'a'
    return f"{article} {color}{exception_name}{return_color}"


def _represent_result(result: Any, color=_result_color) -> str:
    if color is None:
        color = ""
        return_color = ""
    else:
        return_color = _message_color

    if hasattr(result, '__name__'):
        return f"{color}{result.__name__}{return_color}"
    else:
        return f"{color}{repr(result)}{return_color}"


def _process_input(wrapped_function: Callable):
    """
    Allow run_test to take a single positional argument that is a `list`
    or `tuple` in lieu of using all positional and keyword arguments as
    usual.  If `len` of this solitary argument is `3`, then it assumes
    that `kwargs` is an empty `dict` and that the expected
    result/outcome is the last item.
    """
    def decorator(wrapped_function: Callable):
        """"""
        wrapped_signature = inspect.signature(wrapped_function)

        @functools.wraps(wrapped_function)
        def wrapper(*args, **kwargs):

            arguments = wrapped_signature.bind(*args, **kwargs).arguments

            if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], (list, tuple)):
                inputs = args[0]
                if len(inputs) not in (3, 4):
                    raise RuntimeError(f"{args} is an invalid input to run_test.")
                new_kwargs = {'func': inputs[0], 'args': inputs[1]}
                new_kwargs['kwargs'] = inputs[2] if len(inputs) == 4 else {}
                new_kwargs['expected'] = inputs[3] if len(inputs) == 4 else inputs[2]
            else:
                new_kwargs = {argname: argval for argname, argval in arguments.items()}

            return wrapped_function(**new_kwargs)

        return wrapper

    # Allow the decorator to be used either with or without arguments

    if wrapped_function is not None:
        return decorator(wrapped_function)
    else:
        return decorator


@_process_input
def run_test(
        func,
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
    func
        The `callable` to be tested.  The first (and sole) argument to
        `~plasmapy.utils.run_test` may alternatively be a list or tuple
        containing these arguments (optionally omitting `kwargs` if the
        `len` returns 3).

    args: tuple or object
        The positional arguments to `func`.

    kwargs: dict
        The keyword arguments to `func`.

    expected_outcome: object
        The expected result, exception, or warning from
        `func(*args, **kwargs)`. This may also be a `tuple` of length
        two that contains the expected result as the first item and the
        expected warning as the second item.

    rtol : float
        The relative tolerance to be used by `~numpy.allclose` in an
        element-wise comparison, defaulting to `0`.

    atol : float
        The absolute tolerance to be used by `~numpy.allclose` in an
        element-wise comparison, defaulting to `0`.

    Returns
    -------
    `None`

    Raises
    ------
    ~plasmapy.utils.UnexpectedResultError
        If the test returns a result that is different from the expected
        result.

    ~plasmapy.utils.InconsistentTypeError
        If the actual result is of a different type than the expected
        result.

    ~plasmapy.utils.UnexpectedExceptionError
        If an exception occurs when no exception or a different
        exception is expected.

    ~plasmapy.utils.MissingExceptionError
        If no exception is raised when an exception is expected.

    ~plasmapy.utils.MissingWarningError
        An expected warning is not issued.

    ~astropy.units.UnitsError
        If the result has different units than expected.

    TypeError
        If the equality of the actual result and expected result cannot
        be determined (e.g., for a class lacking an `__eq__` method.

    Examples
    --------
    The simplest way to use `~plasmapy.utils.run_test` is with inputs
    for the function to be tests, the positional arguments in a `tuple`
    or `list`, the keyword arguments in a `dict`, and then finally the
    expected result or outcome.

    >>> args = tuple()
    >>> kwargs = dict()
    >>> run_test(lambda: 0, args, kwargs, 0)

    If `expected` is a an exception or warning, then
    `~plasmapy.atomic.run_test` will raise an exception if the expected
    exception is not raised or the expected warning is not issued.

    >>> from warnings import warn

    >>> issue_warning = lambda: warn("Electrons are weird!", UserWarning)
    >>> run_test(issue_warning, args, kwargs, UserWarning)

    >>> def raise_exception(): raise RuntimeError
    >>> run_test(raise_exception, args, kwargs, RuntimeError)

    This function is also flexible enough that it can accept a `tuple`
    or `list` as its sole argument, with the arguments in the same
    order as in the function signature.

    >>> return_arg = lambda x: x
    >>> inputs = (return_arg, 42, {}, 42)
    >>> run_test(inputs)

    If the `tuple` or `list` has a length of `3`, then
    `~plasmapy.utils.run_test` assumes that `kwargs` is missing.

    >>> inputs_without_kwargs = [return_arg, 42, 42]
    >>> run_test(inputs_without_kwargs)

    """

    if not isinstance(args, tuple):
        args = (args,)

    # By including the function call that is run during a test in error
    # messages, we can make it easier to reproduce the error in an
    # interactive session.

    call_str = call_string(func, args, kwargs, color=_func_color, return_color=_message_color)

    # There are many possibilities for expected outcomes that we must
    # keep track of, including exceptions being raised and warnings
    # being issued.

    expected = defaultdict(lambda: None)

    if inspect.isclass(expected_outcome):
        subclass_of_Exception = issubclass(expected_outcome, Exception)
        subclass_of_Warning = issubclass(expected_outcome, Warning)
        if subclass_of_Warning:
            expected['warning'] = expected_outcome
        elif subclass_of_Exception and not subclass_of_Warning:
            expected['exception'] = expected_outcome

    # If a warning is issued, then there may also be an expected result.

    if isinstance(expected_outcome, tuple):
        length_not_two = len(expected_outcome) != 2
        is_not_class = not inspect.isclass(expected_outcome[1])
        is_not_warning = True if is_not_class else not issubclass(expected_outcome[1], Warning)
        if length_not_two or is_not_warning:
            raise ValueError("Invalid expected outcome in run_test.")
        expected['result'] = expected_outcome[0]
        expected['warning'] = expected_outcome[1]

    if expected['exception'] is None and expected['warning'] is None:
        expected['result'] = expected_outcome

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
        except expected_exception as unexpected_exception:
            if unexpected_exception.__reduce__()[0].__name__ == expected_exception.__name__:
                return None
            else:
                raise UnexpectedExceptionError(
                    f"The command {call_str} did not specifically raise "
                    f"{_exc_str(expected_exception)} as expected, but "
                    f"instead raised {_exc_str(unexpected_exception)} "
                    f"which is a subclass of the expected exception.")
        except Exception as exc_unexpected_exception:
            unexpected_exception = exc_unexpected_exception.__reduce__()[0]
            raise UnexpectedExceptionError(
                f"The command {call_str} did not raise "
                f"{_exc_str(expected_exception)} as expected, "
                f"but instead raised {_exc_str(unexpected_exception)}."
            ) from exc_unexpected_exception
        else:
            raise MissingExceptionError(
                f"The command {call_str} did not raise "
                f"{_exc_str(expected_exception)} as expected, but instead "
                f"returned {_represent_result(result)}.")

    try:
        with pytest.warns(expected['warning']):
            result = func(*args, **kwargs)
    except pytest.raises.Exception as missing_warning:
        raise MissingWarningError(
            f"The command {call_str} should issue "
            f"{_exc_str(expected['warning'])}, but instead returned "
            f"{_represent_result(result)}."
        ) from missing_warning
    except Exception as exception_no_warning:
        raise UnexpectedExceptionError(
            f"The command {call_str} unexpectedly raised "
            f"{_exc_str(exception_no_warning.__reduce__()[0])} "
            f"instead of returning the expected value of "
            f"{_represent_result(expected['result'])}."
        ) from exception_no_warning

    if isinstance(expected['result'], u.UnitBase):

        if isinstance(result, u.UnitBase):
            if result != expected['result']:
                raise u.UnitsError(
                    f"The command {call_str} returned "
                    f"{_represent_result(result)} instead of the expected "
                    f"value of {_represent_result(expected['result'])}.")
            return None

        if not isinstance(result, (u.Quantity, const.Constant, const.EMConstant)):
            raise u.UnitsError(
                f"The command {call_str} returned "
                f"{_represent_result(result)} instead of a quantity or "
                f"constant with units of "
                f"{_represent_result(expected['result'])}.")

        if result.unit != expected['result']:
            raise u.UnitsError(
                f"The command {call_str} returned "
                f"{_represent_result(result)}, which has units of "
                f"{result.unit} instead of the expected units of "
                f"{_represent_result(expected['result'])}.")

        return None

    if isinstance(expected['result'], (u.Quantity, const.Constant, const.EMConstant)):
        if not result.unit == expected['result'].unit:
            raise u.UnitsError(
                f"The command {call_str} returned "
                f"{_represent_result(result)} which has different units "
                f"than the expected result of "
                f"{_represent_result(expected['result'])}.")

        if np.allclose(result.value, expected['result'].value):
            return None

    if expected['result'] is None:
        return None

    if type(result) != type(expected['result']):
        raise InconsistentTypeError(
            f"The command {call_str} returned "
            f"{_represent_result(result)} which has type "
            f"{_represent_result(type(result), color=_type_color)}, "
            f"instead of the expected value of "
            f"{_represent_result(expected['result'])} which has type "
            f"{_represent_result(type(expected['result']), color=_type_color)}."
        )

    try:
        if result == expected['result']:
            return None
    except Exception as exc_equality:  # coveralls: ignore
        raise TypeError(
            f"The equality of {_represent_result(result)} and "
            f"{_represent_result(expected['result'])} "
            f"cannot be evaluated.") from exc_equality

    try:
        different_length = len(result) != len(expected['result'])
    except Exception:
        different_length = False

    try:
        all_close = np.allclose(expected['result'], result, rtol=rtol, atol=atol)
        if all_close and not different_length:
            return None
    except Exception:
        pass

    errmsg = (
        f"The command {call_str} returned "
        f"{_represent_result(result)} instead of the expected "
        f"value of {_represent_result(expected['result'])}."
    )

    if atol or rtol:
        errmsg += " with "
        if atol:
            errmsg += f"atol = {atol}"
        if atol and rtol:
            errmsg += " and "
        if rtol:
            errmsg += f"rtol = {rtol}"
    errmsg += "."

    raise UnexpectedResultError(errmsg)
