import warnings
import pytest

from plasmapy.utils.pytest_helpers.actual import ActualTestOutcome
from plasmapy.utils.pytest_helpers.expected import ExpectedTestOutcome
from plasmapy.utils.pytest_helpers.inputs import FunctionTestInputs
from plasmapy.utils.pytest_helpers.comparators import CompareActualExpected


class SampleException(Exception):
    pass


class SampleExceptionSubclass(SampleException):
    pass


class SampleWarning(Warning):
    pass


class SampleWarningSubclass(SampleWarning):
    pass


def return_42():
    return 42


def issue_sample_warning_and_return_42():
    warnings.warn('warning message', SampleWarning)
    return 42


def raise_sample_exception():
    raise SampleException("exception message")


def return_sum_of_two_args_and_kwargs(arg1, arg2, *, kw1=None, kw2=None):
    return arg1 + arg2 + kw1 + kw2


class SampleClass:

    def __init__(self, *args, **kwargs):
        pass

    def arg_plus_kwarg(self, arg, *, kwarg=None):
        return arg + kwarg

    @property
    def forty(self):
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
        (return_42, noargs, nokwargs, SampleException, False),
        (return_42, noargs, nokwargs, (SampleWarning, 42), False),
        (return_42, noargs, nokwargs, (42, SampleWarning), False),

        (issue_sample_warning_and_return_42, noargs, nokwargs, SampleWarning, True),
        (issue_sample_warning_and_return_42, noargs, nokwargs, SampleWarningSubclass, False),
        (issue_sample_warning_and_return_42, noargs, nokwargs, Warning, False),

        (issue_sample_warning_and_return_42, noargs, nokwargs, (SampleWarning, 42), True),
        (issue_sample_warning_and_return_42, noargs, nokwargs, (42, SampleWarning), True),
        (issue_sample_warning_and_return_42, noargs, nokwargs, (SampleWarning, 6 * 9), False),
        (issue_sample_warning_and_return_42, noargs, nokwargs, (6 * 9, SampleWarning), False),

        (issue_sample_warning_and_return_42, noargs, nokwargs, (SampleWarningSubclass, 42), False),
        (issue_sample_warning_and_return_42, noargs, nokwargs, (42, SampleWarningSubclass), False),
        (issue_sample_warning_and_return_42, noargs, nokwargs, (Warning, 42), False),
        (issue_sample_warning_and_return_42, noargs, nokwargs, (42, Warning), False),

        # if unexpected warnings should cause a failure
        (issue_sample_warning_and_return_42, noargs, nokwargs, 42, False),

        (return_sum_of_two_args_and_kwargs, (2, 3), {'kw1': 5, 'kw2': 7}, 17, True),
        (return_sum_of_two_args_and_kwargs, (2, 3), {'kw1': 5, 'kw2': 7}, 42, False),
        (return_sum_of_two_args_and_kwargs, (2, 3), {'kw1': 5, 'kw2': 7}, Warning, False),
        (return_sum_of_two_args_and_kwargs, (2, 3), {'kw1': 5, 'kw2': 7}, Exception, False),
        (return_sum_of_two_args_and_kwargs, (2, 3), {'kw1': 5, 'kw2': 7}, (Warning, 17), False),
        (return_sum_of_two_args_and_kwargs, (2, 3), {'kw1': 5, 'kw2': 7}, (17, Warning), False),
        (return_sum_of_two_args_and_kwargs, (2, 3), {'kw1': 5, 'kw2': 7}, '17', False),

    ])
def test_comparator_actual_expected(function, args, kwargs, expected, test_should_pass):

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
        pytest.fail(
            f"Unable to instantiate CompareActualExpected."
            f"{test_information}")

    test_information += f"\ncomparison.test_passed = {comparison.test_passed}"

    if comparison.error_message:
        test_information += (
            f"\n\n"
            f"Presumably incorrect error message: "
            + comparison.error_message
        )

    if test_should_pass and not comparison.test_passed:
        pytest.fail("Test should have passed, but did not." + test_information)
    if not test_should_pass and comparison.test_passed:
        pytest.fail("Test should not have passed, but did." + test_information)
