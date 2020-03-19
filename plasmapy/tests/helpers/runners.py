"""Functions and classes that run tests."""

__all__ = ["test_runner"]

from typing import NoReturn, Union

from plasmapy.tests.helpers.actual import ActualTestOutcome
from plasmapy.tests.helpers.cases import FunctionTestCase, MethodTestCase, AttrTestCase
from plasmapy.tests.helpers.comparators import CompareActualExpected
from plasmapy.tests.helpers.exceptions import InvalidTestError
from plasmapy.tests.helpers.expected import ExpectedTestOutcome
from plasmapy.tests.helpers.inputs import (
    AbstractTestInputs,
    ClassAttributeTestInputs,
    ClassMethodTestInputs,
    FunctionTestInputs,
)


def _get_test_inputs(
    test_case: Union[FunctionTestCase, AttrTestCase, MethodTestCase]
) -> AbstractTestInputs:
    """
    Select and instantiate the appropriate subclass of
    ~plasmapy.tests.helper.inputs.AbstractTestInputs`.

    Parameters
    ----------
    test_case : `FunctionTestCase`, `AttrTestCase`, or `MethodTestCase`
    """
    if isinstance(test_case, FunctionTestCase):

        inputs = FunctionTestInputs(
            function=test_case.function, args=test_case.args, kwargs=test_case.kwargs
        )

    elif isinstance(MethodTestCase):

        inputs = ClassMethodTestInputs(
            cls=test_case.cls,
            method=test_case.method,
            cls_args=test_case.cls_args,
            cls_kwargs=test_case.cls_kwargs,
            method_args=test_case.method_args,
            method_kwargs=test_case.method_kwargs,
        )

    elif isinstance(AttrTestCase):

        inputs = ClassAttributeTestInputs(
            cls=test_case.cls,
            attribute=test_case.attribute,
            cls_args=test_case.cls_args,
            cls_kwargs=test_case.cls_kwargs,
        )

    else:

        raise TypeError(
            "A test case must be an instance of FunctionTestCase, "
            "MethodTestCase, or AttrTestCase."
        )

    return inputs


def test_runner(
    test_case: Union[FunctionTestCase, MethodTestCase, AttrTestCase]
) -> NoReturn:
    """
    Test that calling a function or method with particular arguments or
    accessing a class attribute results in the expected outcome.

    This test runner works well in combination with `~pytest.mark.parametrize`.

    Parameters
    ----------
    test_case : `FunctionTestCase`, `AttrTestCase`, or `MethodTestCase`

    Raises
    ------
    ~plasmapy.tests.helpers.exceptions.InvalidTestError
        If the test is not set up correctly.

    ~plasmapy.tests.helpers.exceptions.UnexpectedResultError
        If the value returned by the test does not match the expected value.

    ~plasmapy.tests.helpers.exceptions.InconsistentTypeError
        If the type of the value returned by the test is not the same
        type as the expected value.

    ~plasmapy.tests.helpers.exceptions.UnexpectedExceptionError
        If an exception was raised unexpectedly.

    ~plasmapy.tests.helpers.exceptions.MissingExceptionError
        If an exception was expected to be raised, but was not.

    ~plasmapy.tests.helpers.exceptions.ExceptionMismatchError
        If an exception was raised that was different from the expected exception.

    ~plasmapy.tests.helpers.exceptions.UnexpectedWarningError
        If a warning was issued unexpectedly.

    ~plasmapy.tests.helpers.exceptions.MissingWarningError
        If a warning was expected to be issued, but was not.

    ~plasmapy.tests.helpers.exceptions.WarningMismatchError
        If the expected warning was not issued, but one or more different warnings were.

    ~plasmapy.tests.helpers.exceptions.Failed
        If the test fails in a different way or more than one way.

    See Also
    --------
    FunctionTestCase
    MethodTestCase
    AttrTestCase

    Examples
    --------
    This test runner checks that calling a function supplied with certain
    positional and keyword arguments results in the expected outcome. The
    expected outcome may be a value, a warning, an exception, or a
    combination of a warning and a value.

    Suppose we create a function that takes no arguments and returns ``42``.

    >>> return_42 = lambda: 42

    The first step is to specify a test case.  Because we are testing a
    function, we use `FunctionTestCase`.

    >>> test_case = FunctionTestCase(expected=42, function=return_42)

    We then supply ``test_case`` to `test_runner` in order to perform
    the test.

    >>> test_runner(test_case)

    No exception is raised when running the test, which means that this
    test passes.

    Next we wish to test a function that adds the values of a positional
    argument and a keyword argument.

    >>> def add(arg, kwarg=None):
    ...     return arg + kwarg

    Positional arguments are generally given as a `tuple` or `list`
    while keyword arguments are provided as a `dict`.  If there is only
    one positional argument, then it may be supplied as itself (if it is
    not a `tuple` or `list`).  This test also passes.

    >>> args = (5,)
    >>> kwargs = {"kwarg": 4}
    >>> test_add = FunctionTestCase(expected=9, function=add, args=args, kwargs=kwargs)
    >>> test_runner(test_add)

    We may use `test_runner` to check that a function raises an
    appropriate exception.

    >>> def raise_exception():
    ...     raise Exception
    >>> test_exception = FunctionTestCase(expected=Exception, function=raise_exception)
    >>> test_runner(test_exception)

    Similarly, we may test that a function issues a particular warning.

    >>> import warnings
    >>> def issue_warning():
    ...     warnings.warn("warning message", Warning)
    ...     return 42
    >>> test_warning = FunctionTestCase(expected=Warning, function=issue_warning)
    >>> test_runner(test_warning)

    Alternatively, ``expected`` may be a `tuple` containing both the
    expected warning and the expected value (in either order).

    >>> test_warning_and_value = FunctionTestCase(expected=(Warning, 42), function=issue_warning)
    >>> test_runner(test_warning_and_value)

    We may use ``test_runner`` to test that methods and attributes in
    class instances are behaving correctly. Suppose we have a class with
    an attribute named ``attribute`` and a method named ``method``.

    >>> class SampleClass:
    ...     def __init__(self, cls_arg, *, cls_kwarg=None):
    ...         self.cls_arg = cls_arg
    ...         self.cls_kwarg = cls_kwarg
    ...     @property
    ...     def attribute(self):
    ...         return self.cls_arg + self.cls_kwarg
    ...     def method(self, method_arg, *, method_kwarg=None):
    ...         return self.cls_arg + self.cls_kwarg + self.method_arg + self.method_kwarg

    We may set up a test case for an attribute using
    `~plasmapy.tests.helper.cases.AttrTestCase`.

    >>> test_attribute = AttrTestCase(
    ...     cls=SampleClass,
    ...     attribute="attribute",  # name of instance attribute
    ...     cls_args=(1,),  # positional arguments provided to class
    ...     cls_kwargs={"cls_kwarg": 2},  # keyword arguments provided to class
    ...     expected=3,
    ... )
    >>> test_runner(test_attribute)

    We may similarly set up a test case for a method using
    `~plasmapy.tests.helper.cases.MethodTestCase`.

    >>> test_method = MethodTestCase(
    ...     expected=10,
    ...     cls=SampleClass,
    ...     method="method",  # name of instance method
    ...     cls_args=(1,),  # positional arguments provided to class
    ...     cls_kwargs={"cls_kwarg": 2},  # keyword arguments provided to class
    ...     method_args=(3,),  # positional arguments provided to method
    ...     method_kwargs={"method_kwarg": 4},  # keyword arguments provided to method
    ... )
    >>> test_runner(test_method)
    """

    __tracebackhide__ = True

    try:
        test_inputs = _get_test_inputs(test_case)
        actual_outcome = ActualTestOutcome(test_inputs)
        expected_outcome = ExpectedTestOutcome(test_case.expected)
        comparison = CompareActualExpected(
            actual_outcome, expected_outcome, rtol=test_case.rtol, atol=test_case.atol
        )
    except Exception as exc:
        raise InvalidTestError("Unable to run test.") from exc

    if not comparison.test_passed:
        raise comparison.exception(comparison.error_message)
