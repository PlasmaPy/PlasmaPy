import pytest
import warnings
from plasmapy.utils.pytest_helpers.actual import ActualTestOutcome
from plasmapy.utils.pytest_helpers.inputs import FunctionTestInputs


def return_something():
    return 42


def issue_warning():
    warnings.warn(UserWarning)
    return 24


def raise_exception():
    raise RuntimeError


@pytest.mark.parametrize("function, attribute, expected", [
    (return_something, 'warning_was_issued', False),
    (return_something, 'exception_was_raised', False),
    (return_something, 'result', 42),
    (issue_warning, 'warning_was_issued', True),
    (issue_warning, 'exception_was_raised', False),
    (issue_warning, 'result', 24),
    (raise_exception, 'warning_was_issued', False),
    (raise_exception, 'exception_was_raised', True),
])
def test_actual_test_outcome(function, attribute, expected):
    inputs = FunctionTestInputs(function)
    actual = ActualTestOutcome(inputs)
    if getattr(actual, attribute) != expected:
        pytest.fail(
            f"The attribute {repr(attribute)} of ActualTestOutcome for"
            f" function {function.__name__} was {getattr(actual, attribute)}, "
            f"instead of the expected value of {expected}.",
        )
