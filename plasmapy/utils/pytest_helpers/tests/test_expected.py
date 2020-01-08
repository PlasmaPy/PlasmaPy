import collections
import inspect

import pytest

from plasmapy.utils.pytest_helpers.expected import (
    ExpectedTestOutcome,
    _is_warning,
    _is_exception,
    _is_warning_and_value,
)


@pytest.mark.parametrize(
    "argument, expected",
    [(Warning, True), (UserWarning, True), (Exception, False), ("", False),],
)
def test__is_warning(argument, expected):
    """
    Test that `~plasmapy.utils.pytest_helpers.expected._is_warning`
    returns `True` for warnings and `False` for other objects.
    """
    assert _is_warning(argument) is expected


@pytest.mark.parametrize(
    "argument, expected",
    [(Warning, False), (UserWarning, False), (Exception, True), ("", False),],
)
def test__is_exception(argument, expected):
    """
    Test that `~plasmapy.utils.pytest_helpers.expected._is_exception`
    returns `True` for exceptions and `False` for other objects.
    """
    assert _is_exception(argument) is expected


@pytest.mark.parametrize(
    "argument, expected",
    [
        ((Warning, ""), True),
        (["", UserWarning], True),
        ((Warning, UserWarning), False),
        (Warning, False),
        (UserWarning, False),
        (Exception, False),
        ("", False),
    ],
)
def test__is_warning_and_value(argument, expected):
    """
    Test that `~plasmapy.utils.pytest_helpers.expected._is_warning_and_value`
    returns `True` for a `tuple` or `list` containing a warning and an
    object that is not a warning, and `False` for other objects.
    """
    assert _is_warning_and_value(argument) is expected


expected_exception = KeyError
expected_warning = UserWarning
expected_value = 42

Case = collections.namedtuple("Case", ["argument", "attribute", "correct_outcome"])

expected_value_and_warning = (expected_value, expected_warning)
expected_warning_and_value = (expected_warning, expected_value)

cases = [
    Case(expected_exception, "expected_exception", expected_exception),
    Case(expected_exception, "expected_outcome", expected_exception),
    Case(expected_exception, "expecting_an_exception", True),
    Case(expected_exception, "expecting_a_warning", False),
    Case(expected_exception, "expecting_a_value", False),
    Case(expected_warning, "expected_warning", expected_warning),
    Case(expected_warning, "expected_outcome", expected_warning),
    Case(expected_warning, "expecting_an_exception", False),
    Case(expected_warning, "expecting_a_warning", True),
    Case(expected_warning, "expecting_a_value", False),
    Case(expected_value, "expected_value", expected_value),
    Case(expected_value, "expected_outcome", expected_value),
    Case(expected_value, "expecting_an_exception", False),
    Case(expected_value, "expecting_a_warning", False),
    Case(expected_value, "expecting_a_value", True),
    Case(expected_value_and_warning, "expected_warning", expected_warning),
    Case(expected_value_and_warning, "expected_value", expected_value),
    Case(expected_value_and_warning, "expecting_an_exception", False),
    Case(expected_value_and_warning, "expecting_a_warning", True),
    Case(expected_value_and_warning, "expecting_a_value", True),
    Case(expected_value_and_warning, "expected_outcome", expected_warning_and_value),
    Case(expected_warning_and_value, "expected_warning", expected_warning),
    Case(expected_warning_and_value, "expected_value", expected_value),
    Case(expected_warning_and_value, "expecting_an_exception", False),
    Case(expected_warning_and_value, "expecting_a_warning", True),
    Case(expected_warning_and_value, "expecting_a_value", True),
    Case(expected_warning_and_value, "expected_outcome", expected_warning_and_value),
]


@pytest.mark.parametrize("case", cases)
def test_expected_outcome(case):
    expected_outcome = ExpectedTestOutcome(case.argument)
    result = expected_outcome.__getattribute__(case.attribute)
    if result is not case.correct_outcome and result != case.correct_outcome:
        pytest.fail(
            f"ExpectedTestOutcome({case.argument}) results in {repr(result)} "
            f"but should result in {case.correct_outcome}."
        )


exception_to_be_raised = RuntimeError

exception_raising_cases = [
    Case(expected_exception, "expected_warning", exception_to_be_raised),
    Case(expected_exception, "expected_value", exception_to_be_raised),
    Case(expected_warning, "expected_exception", exception_to_be_raised),
    Case(expected_warning, "expected_value", exception_to_be_raised),
    Case(expected_value, "expected_exception", exception_to_be_raised),
    Case(expected_value, "expected_warning", exception_to_be_raised),
    Case(expected_value_and_warning, "expected_exception", exception_to_be_raised),
    Case(expected_warning_and_value, "expected_exception", exception_to_be_raised),
]


@pytest.mark.parametrize("case", exception_raising_cases)
def test_attribute_exceptions(case):
    if not issubclass(case.correct_outcome, Exception):
        raise TypeError(
            "Incorrect test setup: the expected outcome must be an exception."
        )
    with pytest.raises(case.correct_outcome):
        expected_outcome = ExpectedTestOutcome(case.argument)
        result = expected_outcome.__getattribute__(case.attribute)
        print(result)
