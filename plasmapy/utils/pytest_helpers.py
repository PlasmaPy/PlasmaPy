"""Test helper utilities."""
import functools
import pytest
import inspect
import collections
import typing
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


class InvalidTestError(RunTestError):
    """
    Exception for when the inputs to a test are not valid.
    """

def call_string(f: typing.Callable,
                args: typing.Any=tuple(),
                kwargs: typing.Dict={},
                color="",
                return_color="",
                ) -> str:
    """Return a string with the equivalent call of a function."""

    def format_quantity(arg):
        formatted = f'{arg.value}'
        for base, power in zip(arg.unit.bases, arg.unit.powers):
            if power == -1:
                formatted += f"/u.{base}"
            elif power == 1:
                formatted += f"*u.{base}"
            else:
                formatted += f"*u.{base}**{power}"
        return formatted

    def format_arg(arg):
        if hasattr(arg, '__name__'):
            return arg.__name__
        elif isinstance(arg, u.quantity.Quantity):
            return format_quantity(arg)
        else:
            return repr(arg)

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


def _represent_result(result: typing.Any, color=_result_color) -> str:
    if color is None:
        color = ""
        return_color = ""
    else:
        return_color = _message_color

    if hasattr(result, '__name__'):
        return f"{color}{result.__name__}{return_color}"
    else:
        return f"{color}{repr(result)}{return_color}"


def _process_input(wrapped_function: typing.Callable):
    """
    Allow `run_test` to take a single positional argument that is a
    `list` or `tuple` in lieu of using multiple positional/keyword
    arguments as usual.  If `len` of this argument returns `3`, then
    it assumes that `kwargs` is an empty `dict` and that the expected
    result/outcome is the last item.
    """
    def decorator(wrapped_function: typing.Callable):
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
                new_kwargs['expected_outcome'] = inputs[3] if len(inputs) == 4 else inputs[2]
            else:
                new_kwargs = {argname: argval for argname, argval in arguments.items()}
            return wrapped_function(**new_kwargs)

        return wrapper

    return decorator(wrapped_function)


@_process_input
def run_test(
        func,
        args: typing.Any = (),
        kwargs: typing.Dict = {},
        expected_outcome: typing.Any = None,
        rtol: float = 0.0,
        atol: float = 0.0,
        ):
    """
    Test that a function or class returns the expected result, raises
    the expected exception, or issues an expected warning for the
    supplied positional and keyword arguments.

    Parameters
    ----------
    func: callable, list, or tuple
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
    UnexpectedResultError
        If the test returns a result that is different from the expected
        result.

    InconsistentTypeError
        If the actual result is of a different type than the expected
        result.

    UnexpectedExceptionError
        If an exception occurs when no exception or a different
        exception is expected.

    MissingExceptionError
        If no exception is raised when an exception is expected.

    MissingWarningError
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

    For warnings, `~plasmapy.utils.run_test` can accept a `tuple` of two
    items where the first item is the expected result and the second
    item is the expected warning.

    .. code-block:: python

        def return_arg_and_warn(x):
            warn("", UserWarning)
            return x

        run_test(return_arg_and_warn, 1, {}, (1, UserWarning))

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

    .. code-block:: python

        import pytest

        def func(x, raise_exception=False, issue_warning=False):
            if raise_exception:
                raise ValueError("I'm sorry, Dave. I'm afraid I can't do that.")
            elif issue_warning:
                warn("Open the pod bay doors, HAL.", UserWarning)
            return x

        inputs_table = [
            (func, 1, 1),
            (func, (2,), {}, 2),
            (func, 3, {'raise_exception': True}, ValueError),
            (func, 4, {'issue_warning': True}, UserWarning),
            (func, 5, {'issue_warning': True}, (5, UserWarning)),
        ]

        @pytest.mark.parametrize('inputs', inputs_table)
        def test_func(inputs):
            run_test(inputs)

    """

    if not isinstance(args, tuple):
        args = (args,)

    if not callable(func):
        raise InvalidTestError(f"The argument func = {func} to run_test must be callable.")

    # By including the function call that is run during a test in error
    # messages, we can make it easier to reproduce the error in an
    # interactive session.

    call_str = call_string(func, args, kwargs, color=_func_color, return_color=_message_color)

    # There are many possibilities for expected outcomes that we must
    # keep track of, including exceptions being raised and warnings
    # being issued.

    expected = collections.defaultdict(lambda: None)

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
            raise InvalidTestError("Invalid expected outcome in run_test.")
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
        except expected_exception as exc_result:
            resulting_exception = exc_result.__reduce__()[0]
            if resulting_exception.__name__ == expected_exception.__name__:
                return None
            else:
                raise UnexpectedExceptionError(
                    f"The command {call_str} did not specifically raise "
                    f"{_exc_str(expected_exception)} as expected, but "
                    f"instead raised {_exc_str(resulting_exception)} "
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


def run_test_equivalent_calls(*test_inputs, require_same_type: bool = True):
    """
    Test that different functions/inputs return equivalent results.

    Parameters
    ----------
    test_inputs
        The functions and inputs to the tests in an allowed format, as
        described below.

    require_same_type: bool
        If `True` (the default), then all of the results are required to
        be of the same type.  If `False`, results do not need to be of
        the same type (e.g., cases like `1.0 == 1` will not raise an
        exception).

    Raises
    ------
    ~plasmapy.utils.UnexpectedResultError
        If not all of the results are equivalent, or not all of the
        results are of the same type and `require_same_type` evaluates
        to `True`.

    ~plasmapy.utils.UnexpectedExceptionError
        If an exception is raised whilst attempting to run one of the
        test cases.

    ~plasmapy.utils.InvalidTestError
        If there is an error associated with the inputs or the test is
        set up incorrectly.

    Examples
    --------
    There are several possible formats that can be accepted by this
    `~plasmapy.utils.run_test_equivalent_calls` to test that different
    combinations of functions (or other `callable` objects), positional
    arguments, and keyword arguments return equivalent results.

    To test a single function that takes a single positional argument,
    then `test_inputs` may be the function followed by an arbitrary
    number of positional arguments to be included into the function.

    >>> def f(x): return x ** 2
    >>> run_test_equivalent_calls(f, -1, 1)

    To test a single function with an arbitrary number of positional and
    keyword arguments, the first argument should be the function,
    followed by an arbitrary number of `tuple` or `list` objects that
    contain a `tuple` or `list` containing the positional arguments, and
    a `dict` containing the keyword arguments.

    >>> def g(x, y, z): return x + y + z
    >>> run_test_equivalent_calls(g, ((1, 2, 3), {}), ((3, 2), {'z': 1}))

    If there is only one positional argument, then it is not necessary
    to include it in a `tuple` or `list`.

    >>> run_test_equivalent_calls(f, ([1], {}), ([1], {}))
    >>> run_test_equivalent_calls(f, (1, {}), (1, {}))

    To test multiple functions with an arbitrary number of positional
    and keyword arguments, use a series of `tuple` or `list` objects
    that contain the function for each test, a `tuple` or `list` with
    the positional arguments, and a `dict` with the keyword arguments.

    >>> def p(x, y=None): return x + y if y else x
    >>> def q(x, y=None): return x + 1 if y else x

    >>> run_test_equivalent_calls([p, (1,), {'y': 1}], [q, (2,), {'y': False}])

    The inputs may also be passed in as a whole as a `tuple` or `list`.

    >>> run_test_equivalent_calls(f, -1, 1)
    >>> run_test_equivalent_calls([f, -1, 1])

    If `require_same_type` is `False`, then an exception will not be
    raised if the results are of different types.

    >>> run_test_equivalent_calls(f, -1, 1.0, require_same_type=False)

    """

    if len(test_inputs) == 1:
        test_inputs = test_inputs[0]

    if not isinstance(test_inputs, (tuple, list)):
        raise InvalidTestError(
            f"The argument to run_test_equivalent_calls must be a tuple "
            f"or list.  The provided inputs are: {test_inputs}"
        )

    if callable(test_inputs[0]):
        func = test_inputs[0]
        test_inputs = test_inputs[1:]
    else:
        func = None

    # Make sure everything is a list to allow f(*args)

    test_inputs = [test_input if isinstance(test_input, (list, tuple)) else [test_input]
                  for test_input in test_inputs]

    # Construct a list of dicts, of which each dict contains the
    # function, positional arguments, and keyword arguments for each
    # test case.

    test_cases = []

    for inputs in test_inputs:
        test_case = {}

        test_case['function'] = func if func else inputs[0]
        test_case['args'] = inputs[0] if func else inputs[1]

        if not isinstance(test_case['args'], (list, tuple)):
            test_case['args'] = [test_case['args']]

        if func:
            test_case['kwargs'] = inputs[1] if len(inputs) == 2 else {}
        else:
            test_case['kwargs'] = inputs[2] if len(inputs) == 3 else {}

        try:
            test_case['call string'] = call_string(
                test_case['function'], test_case['args'], test_case['kwargs'])
        except Exception:
            test_case['call string'] = (
                f"function = {test_case['function']}, "
                f"args = {test_case['args']}, and "
                f"kwargs = {test_case['kwargs']}")

        test_cases.append(test_case)

    if len(test_cases) < 2:
        raise InvalidTestError("At least two tests are needed for run_test_equivalent_calls")

    # Check to make sure that each function is callable, each set of
    # args is a list or tuple, and each set of kwargs is a dict.  Make
    # sure that the error message contains all of the problems.

    bad_inputs_errmsg = ""

    for test_case in test_cases:
        if not callable(test_case['function']):
            bad_inputs_errmsg += f"\n{test_case['function']} is not callable "
        if not isinstance(test_case['args'], (tuple, list)):
            bad_inputs_errmsg += f"\n{test_case['args']} is not a list or tuple "
        if not isinstance(test_case['kwargs'], dict):
            bad_inputs_errmsg += f"\n{test_case['kwargs']} is not a dict "

    if bad_inputs_errmsg:
        raise InvalidTestError(bad_inputs_errmsg)

    # Now we can get the results for each test case.

    for test_case in test_cases:
        try:
            f, args, kwargs = test_case['function'], test_case['args'], test_case['kwargs']
            test_case['result'] = f(*args, **kwargs)
            test_case['type'] = type(test_case['result'])
        except Exception as exc:
            raise UnexpectedExceptionError(
                f"Unable to evaluate {test_case['call string']}.")

    # Make sure that all of the results evaluate as equal to the first
    # result.

    results = [test_case['result'] for test_case in test_cases]
    types = [test_case['type'] for test_case in test_cases]

    try:
        equals_first_result = [result == results[0] for result in results]
    except Exception as exc:  # coveralls: ignore
        raise UnexpectedExceptionError(
            f"Unable to determine equality properties of results."
        ) from exc

    equals_first_type = [result_type == types[0] for result_type in types]

    not_all_equal = not all(equals_first_result)
    not_all_same_type = not all(equals_first_type)

    if not_all_equal:
        errmsg = f"The following tests did not all produce identical results:"
    elif not_all_same_type and require_same_type:
        errmsg = f"The following tests did not all produce results of the same type:"

    if not_all_equal or (not_all_same_type and require_same_type):

        for test_case in test_cases:
            errmsg += (
                f"\n  {test_case['call string']} yielded {test_case['result']} "
                f"of type {test_case['type']}")

        raise UnexpectedResultError(errmsg)
