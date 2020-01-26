import pytest
import warnings
from plasmapy.utils.pytest_helpers.actual import ActualTestOutcome
from plasmapy.utils.pytest_helpers.inputs import FunctionTestInputs


def return_something():
    return 42


def issue_warning():
    warnings.warn("warning message", UserWarning)
    return 24


def raise_exception():
    raise SyntaxError("exception message")


@pytest.mark.parametrize(
    "function, attribute, expected",
    [
        (return_something, "warning_was_issued", False),
        (return_something, "exception_was_raised", False),
        (return_something, "value", 42),
        (return_something, "value_was_returned", True),
        (return_something, "call_string", "return_something()"),
        (issue_warning, "warning_was_issued", True),
        (issue_warning, "exception_was_raised", False),
        (issue_warning, "value", 24),
        (issue_warning, "value_was_returned", True),
        (issue_warning, "warning_messages", ["warning message"]),
        (issue_warning, "warning_types", [UserWarning]),
        (raise_exception, "warning_was_issued", False),
        (raise_exception, "exception_was_raised", True),
        (raise_exception, "value_was_returned", False),
        (raise_exception, "exception_message", "exception message"),
        (raise_exception, "exception_type", SyntaxError),
    ],
)
def test_actual_test_outcome(function, attribute, expected):
    inputs = FunctionTestInputs(function)
    actual = ActualTestOutcome(inputs)
    if getattr(actual, attribute) != expected:
        pytest.fail(
            f"The attribute {repr(attribute)} of ActualTestOutcome for"
            f" function {function.__name__} was {getattr(actual, attribute)}, "
            f"instead of the expected value of {expected}."
        )


expected_exception = RuntimeError


@pytest.mark.parametrize(
    "function, attribute, expected_exception",
    [
        (return_something, "exception_info", expected_exception),
        (return_something, "warnings_record", expected_exception),
        (return_something, "exception_type", expected_exception),
        (issue_warning, "exception_info", expected_exception),
        (issue_warning, "exception_type", expected_exception),
        (issue_warning, "exception_message", expected_exception),
        (raise_exception, "warnings_record", expected_exception),
        (raise_exception, "value", expected_exception),
        (raise_exception, "warning_messages", expected_exception),
        (raise_exception, "warning_types", expected_exception),
    ],
)
def test_actual_test_outcome_errors(function, attribute, expected_exception):
    inputs = FunctionTestInputs(function)
    actual = ActualTestOutcome(inputs)
    with pytest.raises(expected_exception):
        getattr(actual, attribute)
