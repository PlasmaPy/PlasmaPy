"""Test the helper functions and classes used to run tests."""

from abc import ABC, abstractmethod
from typing import NoReturn, Optional

import astropy.units as u
import numpy as np
import pytest
from plasmapy.tests.helpers.cases import FunctionTestCase, MethodTestCase, AttrTestCase
from plasmapy.tests.helpers.exceptions import (
    ExceptionMismatchError,
    Failed,
    InconsistentTypeError,
    InvalidTestError,
    MissingExceptionError,
    MissingWarningError,
    UnexpectedExceptionError,
    UnexpectedResultError,
    UnexpectedWarningError,
    WarningMismatchError,
)
from plasmapy.tests.helpers.runners import test_runner
from plasmapy.tests.helpers.tests.sample_functions import (
    SampleClass1,
    SampleClass2,
    SampleException,
    SampleExceptionSubclass,
    SampleWarning,
    SampleWarningSubclass,
    issue_warning_return_42,
    raise_exception,
    return_42,
    return_42_meters,
)
from plasmapy.utils.formatting.formatting import _object_name, call_string, class_attribute_call_string, class_method_call_string


cases_and_expected_exceptions = [
    (FunctionTestCase(expected=42, function=return_42), None),
    (
        FunctionTestCase(expected=SampleException, function=return_42),
        MissingExceptionError,
    ),
    (FunctionTestCase(expected=SampleWarning, function=return_42), MissingWarningError),
    (FunctionTestCase(expected=42.0 * u.m, function=return_42_meters), None),
    (FunctionTestCase(expected=u.m, function=return_42_meters), None),
    (FunctionTestCase(expected=42.0 * u.cm, function=return_42_meters), u.UnitsError),
    (FunctionTestCase(expected=u.cm, function=return_42_meters), u.UnitsError),
    (FunctionTestCase(expected=u.kg, function=return_42_meters), u.UnitsError),
    (FunctionTestCase(expected=SampleException, function=raise_exception), None),
    (FunctionTestCase(expected=SampleExceptionSubclass, function=raise_exception),
        ExceptionMismatchError,
    ),
    (
        FunctionTestCase(expected=BaseException, function=raise_exception),
        ExceptionMismatchError,
    ),
    (FunctionTestCase(expected=43, function=return_42), UnexpectedResultError),
    (
        FunctionTestCase(expected=np.int32(42), function=return_42),
        InconsistentTypeError,
    ),
    (FunctionTestCase(expected=42, function=raise_exception), UnexpectedExceptionError),
    (
        FunctionTestCase(expected=42, function=issue_warning_return_42),
        UnexpectedWarningError,
    ),
    (FunctionTestCase(expected=SampleWarning, function=issue_warning_return_42), None),
    (
        FunctionTestCase(expected=Warning, function=issue_warning_return_42),
        WarningMismatchError,
    ),
    (
        FunctionTestCase(
            expected=SampleWarningSubclass, function=issue_warning_return_42
        ),
        WarningMismatchError,
    ),
    (FunctionTestCase(expected=42, function=return_42, kwargs=42), InvalidTestError),
    (FunctionTestCase(expected="..", function=lambda x: 2 * x, args=(".")), None),
    (
        FunctionTestCase(expected="", function=lambda x: 2 * x, args="."),
        UnexpectedResultError,
    ),
    (AttrTestCase(expected=40, cls=SampleClass1, attribute="forty"), None),
    (
        AttrTestCase(expected=41, cls=SampleClass1, attribute="forty"),
        UnexpectedResultError,
    ),
    (
        MethodTestCase(
            expected=5,
            cls=SampleClass1,
            method="arg_plus_kwarg",
            method_args=1,
            method_kwargs={"kwarg": 4},
        ),
        None,
    ),
    (
        MethodTestCase(
            expected=5,
            cls=SampleClass1,
            method="arg_plus_kwarg",
            method_args=(1,),
            method_kwargs={"kwarg": 4},
        ),
        None,
    ),
    (
        MethodTestCase(
            expected=np.int64(36),
            cls=SampleClass2,
            method="method",
            cls_args=(1, 2),
            cls_kwargs={"cls_kwarg1": 3, "cls_kwarg2": 4},
            method_args=(5, 6),
            method_kwargs={"method_kwarg1": 7, "method_kwarg2": 8},
        ),
        None,
    ),
]

# TODO: Add more attribute and method test cases to make sure everything is working


def _get_call_string(case):
    """Interface between ``case`` and the functions that generate call strings."""
    if isinstance(case, FunctionTestCase):
        return call_string(case.function, case.args, case.kwargs)
    elif isinstance(case, MethodTestCase):
        return class_method_call_string(case.cls, case.method, case.cls_args, case.cls_kwargs, case.method_args, case.method_kwargs)
    elif isinstance(case, AttrTestCase):
        return class_attribute_call_string(case.cls, case.attribute, case.cls_args, case.cls_kwargs)


@pytest.mark.parametrize("case, exception", cases_and_expected_exceptions)
def test_the_test_runners(case, exception: Optional[Exception]):
    """
    Test that the test runners pass and fail as they are supposed to.
    """

    test_should_pass = exception is None

    if test_should_pass:

        try:

            test_runner(case)

        except BaseException as exc:

            exception_name = exc.__class__.__name__

            unexpected_exception_errmsg = (
                f"For {_get_call_string(case)} with an expected outcome of "
                f"{_object_name(case.expected)}, the test should have "
                f"passed but instead raised {exception_name}."
            )

            raise Failed(unexpected_exception_errmsg) from exc

    else:

        with pytest.raises(BaseException) as exc_info:

            test_runner(case)

            should_have_failed_errmsg = (
                f"For {_get_call_string(case)} with an expected outcome of "
                f"{case.expected}, the test should have failed but "
                f"instead passed."
            )

            raise Failed(should_have_failed_errmsg)

        if exc_info.type is not exception:

            wrong_exception_errmsg = (
                f"For {_get_call_string(case)} with an expected outcome of "
                f"{case.expected}, the test raised "
                f"{exc_info.typename} but was expected to raise a "
                f"{exception.__name__}."
            )

            raise Failed(wrong_exception_errmsg)
