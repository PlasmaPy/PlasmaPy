import pytest
from plasmapy.tests.helper.actual import ActualTestOutcome
from plasmapy.tests.helper.inputs import FunctionTestInputs

from plasmapy.tests.helper.tests.sample_functions import (
    return_42,
    issue_warning,
    issue_warning_return_42,
    raise_exception,
    SampleWarning,
    SampleException,
)


@pytest.mark.parametrize(
    "function, attribute, expected",
    [
        (return_42, "warning_was_issued", False),
        (return_42, "exception_was_raised", False),
        (return_42, "value", 42),
        (return_42, "value_was_returned", True),
        (return_42, "call_string", f"{return_42.__name__}()"),
        (issue_warning, "warning_was_issued", True),
        (issue_warning, "exception_was_raised", False),
        (issue_warning, "warning_messages", ["warning message"]),
        (issue_warning, "warning_types", [SampleWarning]),
        (issue_warning_return_42, "value", 42),
        (issue_warning_return_42, "value_was_returned", True),
        (raise_exception, "warning_was_issued", False),
        (raise_exception, "exception_was_raised", True),
        (raise_exception, "value_was_returned", False),
        (raise_exception, "exception_message", "exception message"),
        (raise_exception, "exception_type", SampleException),
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
        (return_42, "exception_info", expected_exception),
        (return_42, "warnings_record", expected_exception),
        (return_42, "exception_type", expected_exception),
        (issue_warning, "exception_info", expected_exception),
        (issue_warning, "exception_type", expected_exception),
        (issue_warning, "exception_message", expected_exception),
        (raise_exception, "warnings_record", expected_exception),
        (raise_exception, "value", expected_exception),
        (raise_exception, "warning_messages", expected_exception),
        (raise_exception, "warning_types", expected_exception),
    ],
)
def test_actual_test_outcome_errors(function, attribute, expected_exception: Exception):
    inputs = FunctionTestInputs(function)
    actual = ActualTestOutcome(inputs)
    with pytest.raises(expected_exception):
        getattr(actual, attribute)
