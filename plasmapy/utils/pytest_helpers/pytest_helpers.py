"""Utilities to help with testing."""

__all__ = [
    "assert_can_handle_nparray",
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
    ~plasmapy.tests.helpers.exceptions.UnexpectedResultFail
        If not all of the results are equivalent, or not all of the
        results are of the same type and `require_same_type` evaluates
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

    test_inputs = [
        test_input if isinstance(test_input, (list, tuple)) else [test_input]
        for test_input in test_inputs
    ]

    # Construct a list of dicts, of which each dict contains the
    # function, positional arguments, and keyword arguments for each
    # test case.

    test_cases = []

    for inputs in test_inputs:
        test_case = {}

        test_case["function"] = func if func else inputs[0]
        test_case["args"] = inputs[0] if func else inputs[1]

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
            )

    # Make sure that all of the results evaluate as equal to the first
    # result.

    results = [test_case["result"] for test_case in test_cases]
    types = [test_case["type"] for test_case in test_cases]

    try:
        equals_first_result = [result == results[0] for result in results]
    except Exception as exc:  # coverage: ignore
        raise UnexpectedExceptionFail(
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
                f"of type {test_case['type']}"
            )

        raise UnexpectedResultFail(errmsg)


def assert_can_handle_nparray(
    function_to_test, insert_some_nans=None, insert_all_nans=None, kwargs=None,
):
    """
    Test for ability to handle numpy array quantities.

    Parameters
    ----------
    function_to_test
        The function to be tested for ability to handle numpy array quantities.
        Arguments are automatically given a vector input based on their
        variable name. Current args that are interpreted as vectors are:
        `["T", "T_i", "T_e", "temperature"]`
        `["n", "n_i", "n_e", "density"]`
        `["B"]`
        `["V", "Vperp"]`
        `["coulomb_log"]`
        `["characteristic_length"]`

    insert_some_nans: `list`
        List of argument names in which to insert some np.nan values.
        These must be arguments that will be tested as vectors as listed
        above.

    insert_all_nans: `list`
        List of argument names to fill entirely with np.nan values.

    kwargs: `dict`
        Arguments to pass directly to the function in under test, in the
        normal kwargs python dictionary format.

    Raises
    ------
    ValueError
        If this function cannot interpret a parameter of function_to_test.

    Examples
    --------
    >>> from plasmapy.formulary.parameters import Alfven_speed, gyrofrequency
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
            unit = u.m ** -3
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
            unit = u.m ** -1
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
