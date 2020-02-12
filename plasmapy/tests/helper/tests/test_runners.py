import pytest
from typing import Optional, NoReturn
from abc import ABC, abstractmethod
import numpy as np
import astropy.units as u

from plasmapy.tests.helper.runners import (
    function_test_runner,
    attr_test_runner,
    method_test_runner,
)

from plasmapy.tests.helper.exceptions import (
    Failed,
    UnexpectedResultError,
    InconsistentTypeError,
    UnexpectedExceptionError,
    MissingExceptionError,
    UnexpectedWarningError,
    MissingWarningError,
    InvalidTestError,
    ExceptionMismatchError,
    WarningMismatchError,
)

from plasmapy.tests.helper.tests.sample_functions import (
    return_42,
    return_42_meters,
    issue_warning_return_42,
    raise_exception,
    SampleWarningSubclass,
    SampleWarning,
    SampleExceptionSubclass,
    SampleException,
    SampleClass1,
    SampleClass2,
)

from plasmapy.utils.formatting.formatting import (
    call_string,
    class_attribute_call_string,
    class_method_call_string,
    _name_with_article,
    _object_name,
)


class BaseTestCase(ABC):
    """Define the interface to be used when testing the test cases."""

    @abstractmethod
    def run_test(self) -> NoReturn:

        pass

    @property
    @abstractmethod
    def call_string(self) -> str:

        pass


class FunctionTestCase(BaseTestCase):
    """Contains information for a test case for `function_test_runner`."""

    def __init__(
        self,
        expected,
        function,
        args=None,
        kwargs=None,
        exception_upon_failure: Optional[Exception] = None,
        rtol=1e-8,
        atol=None,
    ):

        self.expected = expected
        self.function = function
        self.args = args if args is not None else ()
        self.kwargs = kwargs if kwargs is not None else {}
        self.exception_upon_failure = exception_upon_failure
        self.test_should_pass = exception_upon_failure is None
        self.rtol = rtol
        self.atol = atol

    def run_test(self) -> NoReturn:
        """Perform the test using `function_test_runner`."""

        function_test_runner(
            expected=self.expected,
            function=self.function,
            args=self.args,
            kwargs=self.kwargs,
            rtol=self.rtol,
            atol=self.atol,
        )

    @property
    def call_string(self) -> str:
        """Return a string to replicate the function call."""

        return call_string(self.function, self.args, self.kwargs)


class AttributeTestCase(BaseTestCase):
    """Contains information for a test case for `attr_test_runner`."""

    def __init__(
        self,
        expected,
        cls,
        attribute: str,
        cls_args=None,
        cls_kwargs=None,
        exception_upon_failure: Optional[Exception] = None,
        rtol=1e-8,
        atol=None,
    ):

        self.expected = expected
        self.cls = cls
        self.attribute = attribute
        self.cls_args = cls_args if cls_args is not None else ()
        self.cls_kwargs = cls_kwargs if cls_kwargs is not None else {}
        self.test_should_pass = exception_upon_failure is None
        self.exception_upon_failure = exception_upon_failure
        self.rtol = rtol
        self.atol = atol

    def run_test(self) -> NoReturn:
        """Perform the test using `attr_test_runner`."""

        attr_test_runner(
            expected=self.expected,
            cls=self.cls,
            attribute=self.attribute,
            cls_args=self.cls_args,
            cls_kwargs=self.cls_kwargs,
            atol=self.atol,
            rtol=self.rtol,
        )

    @property
    def call_string(self) -> str:
        """Return a string to replicate accessing the attribute."""

        return class_attribute_call_string(
            cls=self.cls,
            attr=self.attribute,
            cls_args=self.cls_args,
            cls_kwargs=self.cls_kwargs,
        )


class MethodTestCase(BaseTestCase):
    """Contains information for a test case for `method_test_runner`."""

    def __init__(
        self,
        expected,
        cls,
        method: str,
        cls_args=None,
        cls_kwargs=None,
        method_args=None,
        method_kwargs=None,
        exception_upon_failure: Optional[Exception] = None,
        rtol=1e-8,
        atol=None,
    ):

        self.expected = expected
        self.cls = cls
        self.method = method
        self.cls_args = cls_args if cls_args is not None else ()
        self.cls_kwargs = cls_kwargs if cls_kwargs is not None else {}
        self.method_args = method_args if method_args is not None else ()
        self.method_kwargs = method_kwargs if method_kwargs is not None else {}
        self.test_should_pass = exception_upon_failure is None
        self.exception_upon_failure = exception_upon_failure
        self.rtol = rtol
        self.atol = atol

    def run_test(self) -> NoReturn:
        """Perform the test using `method_test_runner`."""

        __tracebackhide__ = True

        method_test_runner(
            expected=self.expected,
            cls=self.cls,
            method=self.method,
            cls_args=self.cls_args,
            cls_kwargs=self.cls_kwargs,
            method_args=self.method_args,
            method_kwargs=self.method_kwargs,
            rtol=1e-8,
            atol=None,
        )

    @property
    def call_string(self) -> str:
        """Return a string to replicate calling the method."""

        return class_method_call_string(
            cls=self.cls,
            method=self.method,
            cls_args=self.cls_args,
            cls_kwargs=self.cls_kwargs,
            method_args=self.method_args,
            method_kwargs=self.method_kwargs,
        )


cases = [
    FunctionTestCase(
        expected=42,
        function=return_42,
        exception_upon_failure=None,
    ),
    FunctionTestCase(
        expected=SampleException,
        function=return_42,
        exception_upon_failure=MissingExceptionError,
    ),
    FunctionTestCase(
        expected=SampleWarning,
        function=return_42,
        exception_upon_failure=MissingWarningError,
    ),
    FunctionTestCase(
        expected=42.0 * u.m,
        function=return_42_meters,
        exception_upon_failure=None,
    ),
    FunctionTestCase(
        expected=u.m,
        function=return_42_meters,
        exception_upon_failure=None,
    ),
    FunctionTestCase(
        expected=42.0 * u.cm,
        function=return_42_meters,
        exception_upon_failure=u.UnitsError,
    ),
    FunctionTestCase(
        expected=u.cm,
        function=return_42_meters,
        exception_upon_failure=u.UnitsError,
    ),
    FunctionTestCase(
        expected=u.kg,
        function=return_42_meters,
        exception_upon_failure=u.UnitsError,
    ),
    FunctionTestCase(
        expected=SampleException,
        function=raise_exception,
        exception_upon_failure=None,
    ),
    FunctionTestCase(
        expected=SampleExceptionSubclass,
        function=raise_exception,
        exception_upon_failure=ExceptionMismatchError,
    ),
    FunctionTestCase(
        expected=BaseException,
        function=raise_exception,
        exception_upon_failure=ExceptionMismatchError,
    ),
    FunctionTestCase(
        expected=43,
        function=return_42,
        exception_upon_failure=UnexpectedResultError,
    ),
    FunctionTestCase(
        expected=np.int32(42),
        function=return_42,
        exception_upon_failure=InconsistentTypeError,
    ),
    FunctionTestCase(
        expected=42,
        function=raise_exception,
        exception_upon_failure=UnexpectedExceptionError,
    ),
    FunctionTestCase(
        expected=42,
        function=issue_warning_return_42,
        exception_upon_failure=UnexpectedWarningError,
    ),
    FunctionTestCase(
        expected=SampleWarning,
        function=issue_warning_return_42,
        exception_upon_failure=None,
    ),
    FunctionTestCase(
        expected=Warning,
        function=issue_warning_return_42,
        exception_upon_failure=WarningMismatchError,
    ),
    FunctionTestCase(
        expected=SampleWarningSubclass,
        function=issue_warning_return_42,
        exception_upon_failure=WarningMismatchError,
    ),
    FunctionTestCase(
        expected=42,
        function=return_42,
        kwargs=42,
        exception_upon_failure=InvalidTestError,
    ),
    FunctionTestCase(
        expected="..",
        function=lambda x: 2 * x,
        args=("."),
        exception_upon_failure=None,
    ),
    FunctionTestCase(
        expected="",
        function=lambda x: 2 * x,
        args=".",
        exception_upon_failure=UnexpectedResultError,
    ),
    AttributeTestCase(
        expected=40,
        cls=SampleClass1,
        attribute="forty",
        exception_upon_failure=None,
    ),
    AttributeTestCase(
        expected=41,
        cls=SampleClass1,
        attribute="forty",
        exception_upon_failure=UnexpectedResultError,
    ),
    MethodTestCase(
        expected=5,
        cls=SampleClass1,
        method="arg_plus_kwarg",
        method_args=1,
        method_kwargs={"kwarg": 4},
        exception_upon_failure=None,
    ),
    MethodTestCase(
        expected=5,
        cls=SampleClass1,
        method="arg_plus_kwarg",
        method_args=(1,),
        method_kwargs={"kwarg": 4},
        exception_upon_failure=None,
    ),
    MethodTestCase(
        expected=np.int64(36),
        cls=SampleClass2,
        method="method",
        cls_args=(1, 2),
        cls_kwargs={"cls_kwarg1": 3, "cls_kwarg2": 4},
        method_args=(5, 6),
        method_kwargs={"method_kwarg1": 7, "method_kwarg2": 8},
        exception_upon_failure=None,
    ),
]

# TODO: Add more attribute and method test cases to make sure everything is working


@pytest.mark.parametrize("case", cases)
def test_the_test_runners(case: BaseTestCase):
    """
    Test that the test runners pass and fail as they are supposed to.
    """

    if case.test_should_pass:

        try:

            case.run_test()

        except BaseException as exc:

            exception_name = exc.__class__.__name__

            unexpected_exception_errmsg = (
                f"For {case.call_string} with an expected outcome of "
                f"{_object_name(case.expected)}, the test should have "
                f"passed but instead raised {exception_name}."
            )

            raise Failed(unexpected_exception_errmsg) from exc

    else:

        with pytest.raises(BaseException) as exc_info:

            case.run_test()

            should_have_failed_errmsg = (
                f"For {case.call_string} with an expected outcome of "
                f"{case.expected}, the test should have failed but "
                f"instead passed."
            )

            raise Failed(should_have_failed_errmsg)

        if exc_info.type is not case.exception_upon_failure:

            wrong_exception_errmsg = (
                f"For {case.call_string} with an expected outcome of "
                f"{case.expected}, the test raised "
                f"{exc_info.typename} but was expected to raise a "
                f"{case.exception_upon_failure.__name__}."
            )

            raise Failed(wrong_exception_errmsg)
