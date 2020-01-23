import warnings
import pytest
import numpy as np
from astropy import units as u

from plasmapy.utils.pytest_helpers.actual import ActualTestOutcome
from plasmapy.utils.pytest_helpers.expected import ExpectedTestOutcome
from plasmapy.utils.pytest_helpers.inputs import FunctionTestInputs
from plasmapy.utils.pytest_helpers.comparators import CompareActualExpected, CompareValues
from plasmapy.utils.pytest_helpers.exceptions import InvalidTestError


class SampleException(Exception):
    pass


class SampleExceptionSubclass(SampleException):
    pass


class SampleWarning(Warning):
    pass


class SampleWarningSubclass(SampleWarning):
    pass


def return_42() -> int:
    return 42


def return_42_meters() -> u.Quantity:
    return 42.0 * u.m


def return_np_array(*args) -> np.array:
    return np.array(args)


def issue_sample_warning_and_return_42() -> int:
    warnings.warn("warning message", SampleWarning)
    return 42


def raise_sample_exception():
    raise SampleException("exception message")


def return_sum_of_two_args_and_kwargs(arg1, arg2, *, kw1=None, kw2=None):
    return arg1 + arg2 + kw1 + kw2


def return_none():
    return None


class SampleClass:
    def __init__(self, *args, **kwargs):
        pass

    def arg_plus_kwarg(self, arg, *, kwarg=None):
        return arg + kwarg

    @property
    def forty(self) -> int:
        return 40


noargs = ()
nokwargs = {}


@pytest.mark.parametrize(
    "function, args, kwargs, expected, test_should_pass",
    [
        (raise_sample_exception, noargs, nokwargs, SampleException, True),
        (raise_sample_exception, noargs, nokwargs, Exception, False),
        (raise_sample_exception, noargs, nokwargs, SampleExceptionSubclass, False),
        (raise_sample_exception, noargs, nokwargs, 42, False),
        (raise_sample_exception, noargs, nokwargs, SampleWarning, False),
        (return_42, noargs, nokwargs, 42, True),
        (return_42, noargs, nokwargs, 6 * 9, False),
        (return_42, noargs, nokwargs, 42.0, False),
        (return_42, noargs, nokwargs, SampleException, False),
        (return_42, noargs, nokwargs, (SampleWarning, 42), False),
        (return_42, noargs, nokwargs, (42, SampleWarning), False),
        (return_42_meters, noargs, nokwargs, 42.0 * u.m, True),
        (return_42_meters, noargs, nokwargs, 42 * (1 + 1e-15) * u.m, True),
        (return_42_meters, noargs, nokwargs, 4200.0 * u.cm, False),
        (return_42_meters, noargs, nokwargs, 42 * u.kg, False),
        (return_42_meters, noargs, nokwargs, 6 * 9 * u.m, False),
        (return_42_meters, noargs, nokwargs, u.m, True),
        (return_42_meters, noargs, nokwargs, u.cm, False),
        (return_none, noargs, nokwargs, None, True),
        (return_none, noargs, nokwargs, False, False),
        (issue_sample_warning_and_return_42, noargs, nokwargs, SampleWarning, True),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            SampleWarningSubclass,
            False,
        ),
        (issue_sample_warning_and_return_42, noargs, nokwargs, Warning, False),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (SampleWarning, 42),
            True,
        ),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (42, SampleWarning),
            True,
        ),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (SampleWarning, 6 * 9),
            False,
        ),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (6 * 9, SampleWarning),
            False,
        ),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (SampleWarningSubclass, 42),
            False,
        ),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (42, SampleWarningSubclass),
            False,
        ),
        (issue_sample_warning_and_return_42, noargs, nokwargs, (Warning, 42), False),
        (issue_sample_warning_and_return_42, noargs, nokwargs, (42, Warning), False),
        # if unexpected warnings should cause a failure
        (issue_sample_warning_and_return_42, noargs, nokwargs, 42, False),
        (return_sum_of_two_args_and_kwargs, (2, 3), {"kw1": 5, "kw2": 7}, 17, True),
        (return_sum_of_two_args_and_kwargs, (2, 3), {"kw1": 5, "kw2": 7}, 42, False),
        (
            return_sum_of_two_args_and_kwargs,
            (2, 3),
            {"kw1": 5, "kw2": 7},
            Warning,
            False,
        ),
        (
            return_sum_of_two_args_and_kwargs,
            (2, 3),
            {"kw1": 5, "kw2": 7},
            Exception,
            False,
        ),
        (
            return_sum_of_two_args_and_kwargs,
            (2, 3),
            {"kw1": 5, "kw2": 7},
            (Warning, 17),
            False,
        ),
        (
            return_sum_of_two_args_and_kwargs,
            (2, 3),
            {"kw1": 5, "kw2": 7},
            (17, Warning),
            False,
        ),
        (return_sum_of_two_args_and_kwargs, (2, 3), {"kw1": 5, "kw2": 7}, "17", False),
    ],
)
def test_comparator_actual_expected(function, args, kwargs, expected, test_should_pass):
    """
    Test that `CompareActualExpected` correctly determines whether or
    not a test passes.

    Parameters
    ----------
    function
        The sample function to be tested

    args
        A `tuple` or `list` of positional arguments

    kwargs
        A dictionary of keyword arguments

    expected
        The expected outcome of the test to be provided to
        `CompareActualExpected`.  The value of this ``expected`` will
        sometimes be correct and sometimes  be incorrect in order to
        perform test this class.

    test_should_pass
        `True` if ``function(*args, **kwargs)`` should result in
        ``expected``, and `False` otherwise.  Equivalently, this
        argument should be `True` if ``expected`` is correct and
        `False` if ``expected`` is incorrect.
    """

    test_information = (
        "\n\n"
        f"function = {function.__name__}\n"
        f"args = {repr(args)}\n"
        f"kwargs = {repr(kwargs)}\n"
        f"expected = {repr(expected)}\n"
        f"test_should_pass = {test_should_pass}"
    )

    try:
        inputs = FunctionTestInputs(function, args, kwargs)
        actual = ActualTestOutcome(inputs)
        expected = ExpectedTestOutcome(expected)
    except Exception:
        pytest.fail(
            "Problem instantiating classes needed to instantiate "
            "CompareActualExpected" + test_information
        )

    try:
        comparison = CompareActualExpected(actual, expected, atol=None, rtol=1e-7)
    except Exception:
        raise InvalidTestError(
            f"Unable to instantiate CompareActualExpected." f"{test_information}"
        )

    test_information += f"\ncomparison.test_passed = {comparison.test_passed}"

    if comparison.error_message:
        test_information += (
            f"\n\n" f"Presumably incorrect error message: " + comparison.error_message
        )

    if test_should_pass and not comparison.test_passed:
        pytest.fail("Test should have passed, but did not." + test_information)

    if not test_should_pass and comparison.test_passed:
        pytest.fail("Test should not have passed, but did." + test_information)




@pytest.mark.parametrize(
    "this, that, attribute, expected",
    [
        (1, 1, "values", (1, 1)),
        ("1", "1", "types", (str, str)),
        (1, 1, "units", (None, None)),
        (1, 1, "units_are_identical", True),
        (1, 1, "units_are_compatible", True),
        (None, None, "are_identical", True),
        (1, 1, "are_equal", True),
        (0, 0.0, "are_equal", True),
        (1, 1, "are_quantities", False),
        (1, 1, "have_same_types", True),
        (5 * u.m, u.m, "units", (u.m, u.m)),
        (5 * u.m, u.cm, "are_quantity_and_unit", True),
        (5 * u.m, 6 * u.m, "units_are_identical", True),
        (5 * u.m, 6 * u.cm, "units_are_identical", False),
        (5 * u.m, u.cm, "units_are_compatible", True),
        (5 * u.m, 5.000001 * u.m, "are_allclose", True),
        (5 * u.m, 5.01 * u.m, "are_allclose", False),
        ([5.0, 6.0] * u.m, [5.0, 6.0] * u.m, "are_allclose", True),
        ([5.0, 6.0] * u.m, [5.0, 6.00000001] * u.m, "are_allclose", True),
        ([5.0, 6.0] * u.m, [5.0, 6.01] * u.m, "are_allclose", False),
        (5 * u.m, 5.000001 * u.m, "are_equal", False),
        (5 * u.m, 5.01 * u.m, "are_equal", False),
        ([5.0, 6.0] * u.m, [5.0, 6.0] * u.m, "are_equal", True),
        ([5.0, 6.0] * u.m, [5.0, 6.00000001] * u.m, "are_equal", False),
        ([5.0, 6.0] * u.m, [5.0, 6.01] * u.m, "are_equal", False),
    ]
)
def test_compare_values_attributes(this, that, attribute, expected):
    """
    Test that the attributes of CompareValues return the correct results.
    """

    try:
        comparison = CompareValues(this, that, rtol=1e-4)
    except Exception:
        pytest.fail(f"Unable to instantiate CompareValues for {this} and {that}.")

    try:
        value_of_attribute = getattr(comparison, attribute)
    except Exception:
        pytest.fail(
            f"Unable to access attribute {attribute} for the CompareValues "
            f"instance for {this} and {that}.")

    if value_of_attribute != expected and value_of_attribute is not expected:
        pytest.fail(
            f"For the CompareValues instance for {this} and {that}, "
            f"the {repr(attribute)} attribute returns a value of "
            f"{value_of_attribute}, which differs from the expected "
            f"value of {expected}."
        )


@pytest.mark.parametrize(
    "this, that, when_made_boolean",
    [
        (1, 1, True),
        (1, 2, False),
        ("1", "1", True),
        (5 * u.m, 6 * u.m, False),
        (None, None, True),
        (0, 0.0, False),
        (5 * u.m, u.m, True),
        (5 * u.m, u.cm, False),
        (5 * u.m, 6 * u.m, False),
        (5 * u.m, 6 * u.cm, False),
        (5 * u.m, 500 * u.cm, False),
        (5 * u.m, 5.000001 * u.m, True),
        (5 * u.m, 5.01 * u.m, False),
        ([5.0, 6.0] * u.m, [5.0, 6.0] * u.m, True),
        ([5.0, 6.0] * u.m, [5.0, 6.00000001] * u.m, True),
        ([5.0, 6.0] * u.m, [5.0, 6.01] * u.m, False),
        (np.nan, np.nan, True),
        (np.inf, np.inf, True),
        (np.inf, np.nan, False),
        (1, np.int32(1), False),
    ]
)
def test_compare_values_bool(this, that, when_made_boolean):
    """
    Test that CompareValues.__bool__ returns the expected results.
    """

    try:
        comparison = CompareValues(this, that, rtol=1e-4)
    except Exception:
        pytest.fail(f"Unable to instantiate CompareValues for {this} and {that}.")

    try:
        made_boolean = bool(comparison)
    except Exception:
        pytest.fail(
            f"The CompareValues instance for {this} and {that} cannot "
            f"be made boolean."
        )

    if made_boolean is not when_made_boolean:
        pytest.fail(
            f"The CompareValues instance for {this} and {that} corresponds "
            f"to a boolean value of {made_boolean}, when it should "
            f"actually correspond to {when_made_boolean}."
        )
