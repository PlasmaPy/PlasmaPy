"""Utilities to help with testing."""

__all__ = [
    "assert_can_handle_nparray",
    "run_test",
    "run_test_equivalent_calls",
]

import astropy.constants as const
import astropy.tests.helper as astrohelper
import astropy.units as u
import collections
import functools
import inspect
import numpy as np
import pytest
import warnings

from typing import Any, Callable, Dict

from plasmapy.tests.helpers.exceptions import (
    InvalidTestError,
    MissingExceptionFail,
    MissingWarningFail,
    TypeMismatchFail,
    UnexpectedExceptionFail,
    UnexpectedResultFail,
)
from plasmapy.utils.code_repr import _name_with_article, _object_name, call_string
from plasmapy.utils.exceptions import PlasmaPyWarning


def _process_input(wrapped_function: Callable):  # coverage: ignore
    """
    Allow `run_test` to take a single positional argument that is a
    `list` or `tuple` in lieu of using multiple positional/keyword
    arguments as usual.  If `len` of this argument returns `3`, then
    it assumes that `kwargs` is an empty `dict` and that the expected
    result/outcome is the last item.
    """

    def decorator(wrapped_function: Callable):
        wrapped_signature = inspect.signature(wrapped_function)

        @functools.wraps(wrapped_function)
        def wrapper(*args, **kwargs):
            arguments = wrapped_signature.bind(*args, **kwargs).arguments
            if (
                len(args) == 1
                and len(kwargs) == 0
                and isinstance(args[0], (list, tuple))
            ):
                inputs = args[0]
                if len(inputs) not in (3, 4):
                    raise RuntimeError(f"{args} is an invalid input to run_test.")
                new_kwargs = {"func": inputs[0], "args": inputs[1]}
                new_kwargs["kwargs"] = inputs[2] if len(inputs) == 4 else {}
                new_kwargs["expected_outcome"] = (
                    inputs[3] if len(inputs) == 4 else inputs[2]
                )
            else:
                new_kwargs = {argname: argval for argname, argval in arguments.items()}
            return wrapped_function(**new_kwargs)

        return wrapper

    return decorator(wrapped_function)


@_process_input
def run_test(
    func,
    args: Any = (),
    kwargs: Dict = None,
    expected_outcome: Any = None,
    rtol: float = 0.0,
    atol: float = 0.0,
):  # coverage: ignore
    """
    Test that a function or class returns the expected result, raises
    the expected exception, or issues an expected warning for the
    supplied positional and keyword arguments.

    Parameters
    ----------
    func: callable, list, or tuple
        The callable to be tested.  The first (and sole) argument to
        `~plasmapy.utils.pytest_helpers.pytest_helpers.run_test`
        may alternatively be a `list` or `tuple` containing these
        arguments (optionally omitting ``kwargs`` if the `len` returns
        3).

    args: tuple or object
        The positional arguments to ``func``.

    kwargs: dict
        The keyword arguments to ``func``.

    expected_outcome: object
        The expected result, exception, or warning from
        ``func(*args, **kwargs)``. This may also be a `tuple` of length
        two that contains the expected result as the first item and the
        expected warning as the second item.

    rtol : float
        The relative tolerance to be used by `~numpy.allclose` in an
        element-wise comparison, defaulting to ``0``.

    atol : float
        The absolute tolerance to be used by `~numpy.allclose` in an
        element-wise comparison, defaulting to ``0``.

    Returns
    -------
    `None`

    Raises
    ------
    ~plasmapy.tests.helpers.exceptions.UnexpectedResultFail
        If the test returns a result that is different from the expected
        result.

    ~plasmapy.tests.helpers.exceptions.TypeMismatchFail
        If the actual result is of a different type than the expected
        result.

    ~plasmapy.tests.helpers.exceptions.UnexpectedExceptionFail
        If an exception occurs when no exception or a different
        exception is expected.

    ~plasmapy.tests.helpers.exceptions.MissingExceptionFail
        If no exception is raised when an exception is expected.

    ~plasmapy.tests.helpers.exceptions.MissingWarningFail
        An expected warning is not issued.

    ~astropy.units.UnitsError
        If the result has different units than expected.

    TypeError
        If the equality of the actual result and expected result cannot
        be determined (e.g., for a class lacking an ``__eq__`` method.

    Examples
    --------
    The simplest way to use `~plasmapy.utils.pytest_helpers.pytest_helpers.run_test`
    is with inputs for the function to be tests, the positional arguments
    in a `tuple` or `list`, the keyword arguments in a `dict`, and then
    finally the expected result or outcome.

    >>> args = tuple()
    >>> kwargs = dict()
    >>> run_test(lambda: 0, args, kwargs, 0)

    If ``expected`` is an exception or warning, then
    `~plasmapy.utils.pytest_helpers.pytest_helpers.run_test` will raise
    an exception if the expected exception is not raised or the expected
    warning is not issued.

    >>> from warnings import warn

    >>> issue_warning = lambda: warn("Electrons are weird!", UserWarning)
    >>> run_test(issue_warning, args, kwargs, UserWarning)

    >>> def raise_exception(): raise RuntimeError
    >>> run_test(raise_exception, args, kwargs, RuntimeError)

    For warnings, `~plasmapy.utils.pytest_helpers.pytest_helpers.run_test`
    can accept a `tuple` of two items where the first item is the
    expected result and the second item is the expected warning.

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

    If the `tuple` or `list` has a length of ``3``, then
    `~plasmapy.utils.pytest_helpers.pytest_helpers.run_test` assumes
    that ``kwargs`` is missing.

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

    if kwargs is None:
        kwargs = {}

    if not type(args) in [tuple, list]:
        args = (args,)

    if not callable(func):
        raise InvalidTestError(
            f"The argument func = {func} to run_test must be callable."
        )

    # By including the function call that is run during a test in error
    # messages, we can make it easier to reproduce the error in an
    # interactive session.

    call_str = call_string(func, args, kwargs)

    # There are many possibilities for expected outcomes that we must
    # keep track of, including exceptions being raised and warnings
    # being issued.

    expected = collections.defaultdict(lambda: None)

    if inspect.isclass(expected_outcome):
        subclass_of_Exception = issubclass(expected_outcome, Exception)
        subclass_of_Warning = issubclass(expected_outcome, Warning)
        if subclass_of_Warning:
            expected["warning"] = expected_outcome
        elif subclass_of_Exception and not subclass_of_Warning:
            expected["exception"] = expected_outcome

    # If a warning is issued, then there may also be an expected result.

    if isinstance(expected_outcome, tuple):
        length_not_two = len(expected_outcome) != 2
        is_not_class = not inspect.isclass(expected_outcome[1])
        is_not_warning = (
            True if is_not_class else not issubclass(expected_outcome[1], Warning)
        )
        if length_not_two or is_not_warning:
            raise InvalidTestError("Invalid expected outcome in run_test.")
        expected["result"] = expected_outcome[0]
        expected["warning"] = expected_outcome[1]

    if expected["exception"] is None and expected["warning"] is None:
        expected["result"] = expected_outcome

    # First we go through all of the possibilities for when an exception
    # is expected to be raised.  If no exception is raised, then we want
    # an error message that includes the result.  If the wrong exception
    # is raised, then we want an error message that includes that
    # exception.  An alternative would be to use `with pytest.raises()`
    # but this makes it easier to break down what the error messages
    # should be.

    if expected["exception"]:

        expected_exception = expected["exception"]

        try:
            result = func(*args, **kwargs)
        except expected_exception as exc_result:
            resulting_exception = exc_result.__reduce__()[0]
            if resulting_exception.__name__ == expected_exception.__name__:
                return None
            else:
                raise UnexpectedExceptionFail(
                    f"The command {call_str} did not specifically raise "
                    f"{_name_with_article(expected_exception)} as expected, but "
                    f"instead raised {_name_with_article(resulting_exception)} "
                    f"which is a subclass of the expected exception."
                ) from exc_result
        except Exception as exc_unexpected_exception:
            unexpected_exception = exc_unexpected_exception.__reduce__()[0]
            raise UnexpectedExceptionFail(
                f"The command {call_str} did not raise "
                f"{_name_with_article(expected_exception)} as expected, "
                f"but instead raised {_name_with_article(unexpected_exception)}."
            ) from exc_unexpected_exception
        else:
            raise MissingExceptionFail(
                f"The command {call_str} did not raise "
                f"{_name_with_article(expected_exception)} as expected, but instead "
                f"returned {_object_name(result)}."
            )

    try:
        with pytest.warns(expected["warning"]):
            result = func(*args, **kwargs)
    except pytest.raises.Exception as missing_warning:
        raise MissingWarningFail(
            f"The command {call_str} should issue "
            f"{_name_with_article(expected['warning'])}, but instead returned "
            f"{_object_name(result)}."
        ) from missing_warning
    except Exception as exception_no_warning:
        raise UnexpectedExceptionFail(
            f"The command {call_str} unexpectedly raised "
            f"{_name_with_article(exception_no_warning.__reduce__()[0])} "
            f"instead of returning the expected value of "
            f"{_object_name(expected['result'])}."
        ) from exception_no_warning

    if isinstance(expected["result"], u.UnitBase):

        if isinstance(result, u.UnitBase):
            if result != expected["result"]:
                raise u.UnitsError(
                    f"The command {call_str} returned "
                    f"{_object_name(result)} instead of the expected "
                    f"value of {_object_name(expected['result'])}."
                )
            return None

        if not isinstance(result, (u.Quantity, const.Constant, const.EMConstant)):
            raise u.UnitsError(
                f"The command {call_str} returned "
                f"{_object_name(result)} instead of a quantity or "
                f"constant with units of "
                f"{_object_name(expected['result'])}."
            )

        if result.unit != expected["result"]:
            raise u.UnitsError(
                f"The command {call_str} returned "
                f"{_object_name(result)}, which has units of "
                f"{result.unit} instead of the expected units of "
                f"{_object_name(expected['result'])}."
            )

        return None

    if isinstance(expected["result"], (u.Quantity, const.Constant, const.EMConstant)):
        if not result.unit == expected["result"].unit:
            raise u.UnitsError(
                f"The command {call_str} returned "
                f"{_object_name(result)} which has different units "
                f"than the expected result of "
                f"{_object_name(expected['result'])}."
            )

        if np.allclose(result.value, expected["result"].value):
            return None

    if expected["result"] is None:
        return None

    if type(result) != type(expected["result"]):
        raise TypeMismatchFail(
            f"The command {call_str} returned "
            f"{_object_name(result)} which has type "
            f"{_object_name(type(result))}, "
            f"instead of the expected value of "
            f"{_object_name(expected['result'])} which has type "
            f"{_object_name(type(expected['result']))}."
        )

    try:
        if result == expected["result"]:
            return None
    except Exception as exc_equality:  # coverage: ignore
        raise TypeError(
            f"The equality of {_object_name(result)} and "
            f"{_object_name(expected['result'])} "
            f"cannot be evaluated."
        ) from exc_equality

    try:
        different_length = len(result) != len(expected["result"])
    except Exception:
        different_length = False

    try:
        all_close = np.allclose(expected["result"], result, rtol=rtol, atol=atol)
        if all_close and not different_length:
            return None
    except Exception:
        pass

    errmsg = (
        f"The command {call_str} returned "
        f"{_object_name(result)} instead of the expected "
        f"value of {_object_name(expected['result'])}."
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

    raise UnexpectedResultFail(errmsg)


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
        the same type (e.g., cases like ``1.0 == 1`` will not raise an
        exception).

    Raises
    ------
    ~plasmapy.tests.helpers.exceptions.UnexpectedResultFail
        If not all of the results are equivalent, or not all of the
        results are of the same type and ``require_same_type`` evaluates
        to `True`.

    ~plasmapy.tests.helpers.exceptions.UnexpectedExceptionFail
        If an exception is raised whilst attempting to run one of the
        test cases.

    ~plasmapy.tests.helpers.exceptions.InvalidTestError
        If there is an error associated with the inputs or the test is
        set up incorrectly.

    Examples
    --------
    There are several possible formats that can be accepted by this
    `~plasmapy.utils.pytest_helpers.pytest_helpers.run_test_equivalent_calls`
    to test that different combinations of functions (or other callable
    objects), positional arguments, and keyword arguments return
    equivalent results.

    To test a single function that takes a single positional argument,
    then ``test_inputs`` may be the function followed by an arbitrary
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

    If ``require_same_type``  is `False`, then an exception will not be
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

    test_inputs = [
        test_input if isinstance(test_input, (list, tuple)) else [test_input]
        for test_input in test_inputs
    ]

    # Construct a list of dicts, of which each dict contains the
    # function, positional arguments, and keyword arguments for each
    # test case.

    test_cases = []

    for inputs in test_inputs:
        test_case = {
            "function": func or inputs[0],
            "args": inputs[0] if func else inputs[1],
        }

        if not isinstance(test_case["args"], (list, tuple)):
            test_case["args"] = [test_case["args"]]

        if func:
            test_case["kwargs"] = inputs[1] if len(inputs) == 2 else {}
        else:
            test_case["kwargs"] = inputs[2] if len(inputs) == 3 else {}

        try:
            test_case["call string"] = call_string(
                test_case["function"], test_case["args"], test_case["kwargs"]
            )
        except Exception:
            test_case["call string"] = (
                f"function = {test_case['function']}, "
                f"args = {test_case['args']}, and "
                f"kwargs = {test_case['kwargs']}"
            )

        test_cases.append(test_case)

    if len(test_cases) < 2:
        raise InvalidTestError(
            "At least two tests are needed for run_test_equivalent_calls"
        )

    # Check to make sure that each function is callable, each set of
    # args is a list or tuple, and each set of kwargs is a dict.  Make
    # sure that the error message contains all of the problems.

    bad_inputs_errmsg = ""

    for test_case in test_cases:
        if not callable(test_case["function"]):
            bad_inputs_errmsg += f"\n{test_case['function']} is not callable "
        if not isinstance(test_case["args"], (tuple, list)):
            bad_inputs_errmsg += f"\n{test_case['args']} is not a list or tuple "
        if not isinstance(test_case["kwargs"], dict):
            bad_inputs_errmsg += f"\n{test_case['kwargs']} is not a dict "

    if bad_inputs_errmsg:
        raise InvalidTestError(bad_inputs_errmsg)

    # Now we can get the results for each test case.

    for test_case in test_cases:
        try:
            f, args, kwargs = (
                test_case["function"],
                test_case["args"],
                test_case["kwargs"],
            )
            test_case["result"] = f(*args, **kwargs)
            test_case["type"] = type(test_case["result"])
        except Exception as exc:
            raise UnexpectedExceptionFail(
                f"Unable to evaluate {test_case['call string']}."
            ) from exc

    # Make sure that all of the results evaluate as equal to the first
    # result.

    results = [test_case["result"] for test_case in test_cases]
    types = [test_case["type"] for test_case in test_cases]

    try:
        equals_first_result = [result == results[0] for result in results]
    except Exception as exc:  # coverage: ignore
        raise UnexpectedExceptionFail(
            "Unable to determine equality properties of results."
        ) from exc

    equals_first_type = [result_type == types[0] for result_type in types]

    not_all_equal = not all(equals_first_result)
    not_all_same_type = not all(equals_first_type)

    if not_all_equal:
        errmsg = "The following tests did not all produce identical results:"
    elif not_all_same_type and require_same_type:
        errmsg = "The following tests did not all produce results of the same type:"

    if not_all_equal or (not_all_same_type and require_same_type):

        for test_case in test_cases:
            errmsg += (
                f"\n  {test_case['call string']} yielded {test_case['result']} "
                f"of type {test_case['type']}"
            )

        raise UnexpectedResultFail(errmsg)


def assert_can_handle_nparray(
    function_to_test,
    insert_some_nans=None,
    insert_all_nans=None,
    kwargs=None,
):
    """
    Test for ability to handle numpy array quantities.

    Parameters
    ----------
    function_to_test
        The function to be tested for ability to handle numpy array quantities.
        Arguments are automatically given a vector input based on their
        variable name. Current args that are interpreted as vectors are:
        ``["T", "T_i", "T_e", "temperature"]``,
        ``["n", "n_i", "n_e", "density"]``,
        ``["B"]``,
        ``["V", "Vperp"]``,
        ``["coulomb_log"]``,
        ``["characteristic_length"]``.

    insert_some_nans: `list`
        List of argument names in which to insert some `~numpy.nan`
        values. These must be arguments that will be tested as vectors
        as listed above.

    insert_all_nans: `list`
        List of argument names to fill entirely with `~numpy.nan` values.

    kwargs: `dict`
        Arguments to pass directly to the function in under test, in the
        normal kwargs python dictionary format.

    Raises
    ------
    ValueError
        If this function cannot interpret a parameter of ``function_to_test``.

    Examples
    --------
    >>> from plasmapy.formulary import Alfven_speed, gyrofrequency
    >>> assert_can_handle_nparray(Alfven_speed)
    >>> assert_can_handle_nparray(gyrofrequency, kwargs={"signed": True})
    >>> assert_can_handle_nparray(gyrofrequency, kwargs={"signed": False})
    """

    if insert_some_nans is None:
        insert_some_nans = []

    if insert_all_nans is None:
        insert_all_nans = []

    if kwargs is None:
        kwargs = {}

    def _prepare_input(
        param_name, param_default, insert_some_nans, insert_all_nans, kwargs
    ):
        """
        Parse parameter names and set up values to input for 0d, 1d, and 2d array tests.
        """
        # first things first: let any passed in kwarg right through (VIP access)
        if param_name in kwargs.keys():
            return (kwargs[param_name],) * 4

        # else, if it's a recognized variable name, give it a reasonable unit and magnitude
        elif param_name in ["particle", "ion_particle", "ion"]:
            if not (param_default is inspect._empty or param_default is None):
                return (param_default,) * 4
            else:
                return ("p",) * 4
        elif param_name == "particles" or param_name == "species":
            if not (param_default is inspect._empty):
                return (param_default,) * 4
            else:
                return (("e", "p"),) * 4
        elif param_name in ["T", "T_i", "T_e", "temperature"]:
            unit = u.eV
            magnitude = 1.0
        elif param_name in ["n", "n_i", "n_e", "density"]:
            unit = u.m**-3
            magnitude = 1e20
        elif param_name == "B":
            unit = u.G
            magnitude = 1e3
        elif param_name in ["V", "Vperp"]:
            unit = u.m / u.s
            magnitude = 1e5
        elif param_name == "coulomb_log":
            unit = 1.0
            magnitude = 1e1
        elif param_name == "characteristic_length":
            unit = u.m
            magnitude = 1.0
        elif param_name == "k":
            unit = u.m**-1
            magnitude = 1.0

        # else, last resort, if it has a default argument, go with that:
        elif not (param_default is inspect._empty):
            return (param_default,) * 4

        else:
            raise ValueError(f"Unrecognized function input: {param_name}")

        # now knowing unit and magnitude, set up the 0d, 1d, 2d, and 3d arrays:
        input_data_3d = np.reshape(np.arange(1.0, 9.0, 1.0), (2, 2, 2))
        input_data_2d = np.reshape(np.arange(1.0, 5.0, 1.0), (2, 2))
        input_data_1d = np.arange(1.0, 5.0, 1.0)
        if param_name in insert_some_nans:
            input_data_3d[0, 0, 1] = np.nan
            input_data_3d[0, 1, 0] = np.nan
            input_data_2d[0, 1] = np.nan
            input_data_2d[1, 0] = np.nan
            input_data_1d[1] = np.nan
        elif param_name in insert_all_nans:
            input_data_3d = np.ones((2, 2, 2)) * np.nan
            input_data_2d = np.ones((2, 2)) * np.nan
            input_data_1d = np.ones(4) * np.nan
        input_data_3d *= magnitude
        input_data_3d *= unit
        input_data_2d *= magnitude
        input_data_2d *= unit
        input_data_1d *= magnitude
        input_data_1d *= unit
        input_data_0d = input_data_1d[3]
        return input_data_0d, input_data_1d, input_data_2d, input_data_3d

    # call _prepare_input to prepare 0d, 1d, and 2d sets of arguments for the function:
    function_sig = inspect.signature(function_to_test)
    function_params = function_sig.parameters
    args_0d = dict()
    args_1d = dict()
    args_2d = dict()
    args_3d = dict()
    param_names = [elm for elm in function_params.keys()]
    for idx, key in enumerate(function_params):
        args_0d[key], args_1d[key], args_2d[key], args_3d[key] = _prepare_input(
            param_names[idx],
            function_params[key].default,
            insert_some_nans,
            insert_all_nans,
            kwargs,
        )

    # call the function with the prepared argument sets:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=PlasmaPyWarning)
        result_0d = function_to_test(**args_0d)
        result_1d = function_to_test(**args_1d)
        result_2d = function_to_test(**args_2d)
        result_3d = function_to_test(**args_3d)

    # assert that the 1d, 2d, 3d versions get the same result (elementwise) as the 0d version:
    # (if the function returns multiple values, loop through and test each)
    try:
        scalar_testable = result_0d.value
    except AttributeError:
        scalar_testable = result_0d
    if np.isscalar(scalar_testable):
        astrohelper.assert_quantity_allclose(result_0d, result_1d[3])
        astrohelper.assert_quantity_allclose(result_0d, result_2d[1, 1])
        astrohelper.assert_quantity_allclose(result_0d, result_3d[0, 1, 1])
    else:
        for idx, res_0d in enumerate(result_0d):
            astrohelper.assert_quantity_allclose(res_0d, result_1d[idx][3])
            astrohelper.assert_quantity_allclose(res_0d, result_2d[idx][1, 1])
            astrohelper.assert_quantity_allclose(res_0d, result_3d[idx][0, 1, 1])
