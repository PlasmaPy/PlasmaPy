"""Test the classes that represent the actual result of a test."""

from typing import Callable

import pytest
from plasmapy.tests.helpers.actual import ActualTestOutcome
from plasmapy.tests.helpers.exceptions import InvalidTestError
from plasmapy.tests.helpers.inputs import FunctionTestInputs
from plasmapy.tests.helpers.tests.sample_functions import (
    SampleException,
    SampleWarning,
    issue_warning,
    issue_warning_return_42,
    raise_exception,
    return_42,
)

func_attr_expected = [
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
]


@pytest.mark.parametrize("function, attr, expected", func_attr_expected)
def test_actual_test_outcome(function: Callable, attr: str, expected):
    """Test that `ActualTestOutcome` attributes return the expected values."""

    inputs = FunctionTestInputs(function)
    actual = ActualTestOutcome(inputs)

    if getattr(actual, attr) != expected:
        pytest.fail(
            f"The attribute {repr(attr)} of ActualTestOutcome for"
            f" function {function.__name__} was {getattr(actual, attr)}, "
            f"instead of the expected value of {expected}."
        )


function_and_exception_raising_attribute = [
    (return_42, "exception_info"),
    (return_42, "warnings_record"),
    (return_42, "exception_type"),
    (issue_warning, "exception_info"),
    (issue_warning, "exception_type"),
    (issue_warning, "exception_message"),
    (raise_exception, "warnings_record"),
    (raise_exception, "value"),
    (raise_exception, "warning_messages"),
    (raise_exception, "warning_types"),
]


@pytest.mark.parametrize(
    "function, attribute", function_and_exception_raising_attribute
)
def test_actual_test_outcome_errors(function: Callable, attribute: str):
    """Test that `ActualTestOutcome` attributes raise exceptions as appropriate."""

    inputs = FunctionTestInputs(function)
    actual = ActualTestOutcome(inputs)

    with pytest.raises(InvalidTestError):
        getattr(actual, attribute)
