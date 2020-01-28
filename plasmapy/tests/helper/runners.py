import pytest

from typing import Callable, Union, Optional, Any, AnyStr, Dict
from numbers import Number

from astropy import units as u

from plasmapy.tests.helper.inputs import (
    AbstractTestInputs,
    FunctionTestInputs,
    ClassAttributeTestInputs,
    ClassMethodTestInputs,
)

from plasmapy.tests.helper.actual import ActualTestOutcome
from plasmapy.tests.helper.expected import ExpectedTestOutcome

from plasmapy.tests.helper.comparators import CompareActualExpected

from plasmapy.tests.helper.exceptions import InvalidTestError


__all__ = ["function_test_runner", "method_test_runner", "attr_test_runner"]


def _test_runner(inputs: AbstractTestInputs, expected, *, rtol=1e-8, atol=None):
    """
    Perform the parts of the test that are common among the different
    test runners for the result test pattern.
    """

    __tracebackhide__ = True

    try:
        actual_outcome = ActualTestOutcome(inputs)
        expected_outcome = ExpectedTestOutcome(expected)
        comparison = CompareActualExpected(actual_outcome, expected_outcome, rtol=rtol, atol=atol)
    except Exception as exc:
        raise InvalidTestError("Unable to run test.") from exc

    if not comparison.test_passed:
        raise comparison.exception(comparison.error_message)


def function_test_runner(
    expected,
    function: Callable,
    args=None,
    kwargs: Optional[Dict[AnyStr, Any]] = None,
    *,
    rtol: Union[Number, u.Quantity] = 1e-8,
    atol: Optional[Union[Number, u.Quantity]] = None,
):
    """
    Test that calling a function with particular arguments results in
    the expected outcome.

    Parameters
    ----------
    expected
        The expected outcome of the test, which can be an exception,
        a warning, the resulting object, or a tuple that contains
        a warning and the resulting object.

    function
        The function to be tested.

    args : optional
        The positional arguments to be provided to the function that is
        being tested.  If this is a `tuple` or `list`, then the
        arguments will be each of the items in the collection. If there
        is only one positional argument, then it may be inputted as is
        without being put into a `tuple` or `list`.

    kwargs : `dict`, optional
        The keyword arguments to be provided to the function that is
        being tested. If provided, the keys of this `dict` must be
        strings.

    rtol : dimensionless number, optional, keyword-only
        The relative tolerance to be used in `astropy.units.allclose`.
        Defaults to ``1e-8``.  If ``rtol`` is a `~astropy.units.Quantity`
        instance, then it must be dimensionless.

    atol : number or `~astropy.units.Quantity`, optional, keyword-only
        The absolute tolerance to be supplied to `astropy.units.isclose`
        or `~astropy.units.allclose`.  If ``atol`` is a
        `~astropy.units.Quantity`, then it must have the same units as
        ``this`` and ``that``.  Defaults to zero in the appropriate units.

    Raises
    ------
    InvalidTestError
        If the test is not set up correctly.
    """

    __tracebackhide__ = True
    inputs = FunctionTestInputs(function, args, kwargs)
    _test_runner(inputs, expected, rtol=rtol, atol=atol)


def method_test_runner(
    expected,
    cls,
    method: AnyStr,
    *,
    cls_args=None,
    cls_kwargs: Optional[Dict[AnyStr, Any]] = None,
    method_args=None,
    method_kwargs: Optional[Dict[str, Any]] = None,
    rtol: Union[Number, u.Quantity] = 1e-8,
    atol: Optional[Union[Number, u.Quantity]] = None,
):
    """
    Test that calling a class method results in the expected outcome.

    Parameters
    ----------
    expected
        The expected outcome of the test, which can be an exception,
        a warning, the resulting object, or a tuple that contains
        a warning and the resulting object.

    cls
        The class containing the method to be tested.

    cls_args : optional
        The positional arguments to be provided to the class that is
        being tested during instantiation. If this is a `tuple` or
        `list`, then the arguments will be each of the items in the
        collection. If there is only one positional argument, then it
        may be inputted as is without being put into a `tuple` or `list`.

    cls_kwargs : `dict`, optional
        The keyword arguments to be provided to the class that is
        being tested during instantiation. If provided, the keys of
        this `dict` must be strings.

    method : str
        The name of the method contained in ``cls`` to be tested.

    method_args : optional
        The positional arguments to be provided to the method that is
        being tested.  If this is a `tuple` or `list`, then the
        arguments will be each of the items in the collection. If there
        is only one positional argument, then it may be inputted as is
        without being put into a `tuple` or `list`.

    method_kwargs : `dict`, optional
        The keyword arguments to be provided to the method that is being
        tested. If provided, the keys of this `dict` must be strings.

    rtol : dimensionless number, optional, keyword-only
        The relative tolerance to be used in `astropy.units.allclose`.
        Defaults to ``1e-8``.  If ``rtol`` is a `~astropy.units.Quantity`
        instance, then it must be dimensionless.

    atol : number or `~astropy.units.Quantity`, optional, keyword-only
        The absolute tolerance to be supplied to `astropy.units.isclose`
        or `~astropy.units.allclose`.  If ``atol`` is a
        `~astropy.units.Quantity`, then it must have the same units as
        ``this`` and ``that``.  Defaults to zero in the appropriate units.

    Raises
    ------
    InvalidTestError
        If the test is not set up correctly.
    """

    __tracebackhide__ = True
    inputs = ClassMethodTestInputs(cls, method, cls_args, cls_kwargs, method_args, method_kwargs)
    _test_runner(inputs, expected, rtol=rtol, atol=atol)


def attr_test_runner(
    expected,
    cls,
    attribute: AnyStr,
    cls_args=None,
    cls_kwargs: Optional[Dict[AnyStr, Any]] = None,
    *,
    rtol: Union[Number, u.Quantity] = 1e-8,
    atol: Optional[Union[Number, u.Quantity]] = None,
):
    """
    Test that accessing a class attribute results in the expected outcome.

    Parameters
    ----------
    expected
        The expected outcome of the test, which can be an exception,
        a warning, the resulting object, or a tuple that contains
        a warning and the resulting object.

    cls
        The class containing the method to be tested.

    attribute : str
        The name of the attribute contained in ``cls`` to be tested.

    cls_args : optional
        The positional arguments to be provided to the class that is
        being tested during instantiation. If this is a `tuple` or
        `list`, then the arguments will be each of the items in the
        collection. If there is only one positional argument, then it
        may be inputted as is without being put into a `tuple` or `list`.

    cls_kwargs : `dict`, optional
        The keyword arguments to be provided to the class that is
        being tested during instantiation. If provided, the keys of
        this `dict` must be strings.

    rtol : dimensionless number, optional, keyword-only
        The relative tolerance to be used in `astropy.units.allclose`.
        Defaults to ``1e-8``.  If ``rtol`` is a `~astropy.units.Quantity`
        instance, then it must be dimensionless.

    atol : number or `~astropy.units.Quantity`, optional, keyword-only
        The absolute tolerance to be supplied to `astropy.units.isclose`
        or `~astropy.units.allclose`.  If ``atol`` is a
        `~astropy.units.Quantity`, then it must have the same units as
        ``this`` and ``that``.  Defaults to zero in the appropriate units.
    """

    __tracebackhide__ = True
    inputs = ClassAttributeTestInputs(cls, attribute, cls_args, cls_kwargs)
    _test_runner(inputs, expected, rtol=rtol, atol=atol)
