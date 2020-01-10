import warnings
import pytest
import numpy as np
from astropy import units as u

from plasmapy.utils.pytest_helpers.actual import ActualTestOutcome
from plasmapy.utils.pytest_helpers.expected import ExpectedTestOutcome
from plasmapy.utils.pytest_helpers.inputs import FunctionTestInputs
from plasmapy.utils.pytest_helpers.comparators import CompareActualExpected
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
        (return_42_meters, noargs, nokwargs, 42.00000000000001 * u.m, True),
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
        comparison = CompareActualExpected(actual, expected)
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


def test_comparator_atol():
    raise NotImplementedError


def test_comparator_rtol():
    raise NotImplementedError
