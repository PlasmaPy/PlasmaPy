"""Classes to contain all inputs for a test case and the expected outcome."""

__all__ = ["FunctionTestCase", "MethodTestCase", "AttrTestCase"]

from abc import ABC, abstractmethod
from numbers import Real
from typing import Any, Callable, Dict, Optional


class _AbstractTestCase(ABC):
    """An interface for storing the inputs and expected result of a test case."""

    @property
    @abstractmethod
    def expected(self):
        pass

    @property
    @abstractmethod
    def rtol(self):
        pass

    @property
    @abstractmethod
    def atol(self):
        pass


class FunctionTestCase(_AbstractTestCase):
    """
    Stores objects that define a test case for a function.

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

    Examples
    --------
    `FunctionTestCase` stores the objects needed by
    `~plasmapy.tests.helpers.test_runner` to test that calling
    a function with certain positional and/or keyword arguments results
    in the expected outcome (which could be a value, an exception, a
    warning, or a unit).

    This example demonstrates how to set up and run a test case to check
    that a function that takes no arguments returns a particular value.

    >>> from plasmapy.tests.helpers import FunctionTestCase, test_runner
    >>> def return_42():
    ...     return 42
    >>> test_case1 = FunctionTestCase(function=return_42, expected=42)
    >>> test_runner(test_case1)

    We may similarly specify a test case for a function that takes
    positional and keyword arguments.

    >>> def return_args_plus_kwarg(arg1, arg2, *, kwarg=0):
    ...     return arg1 + arg2 + kwarg
    >>> test_case2 = FunctionTestCase(
    ...     function=return_args_plus_kwarg,
    ...     args=(1, 2),
    ...     kwargs={"kwarg": 3},
    ...     expected=6,
    ... )
    >>> test_runner(test_case2)

    For tests with floating point operations, we may specify a relative
    tolerance with the ``rtol`` keyword and an absolute tolerance with
    the ``atol`` keyword.

    >>> from astropy import units as u
    >>> def double(x):
    ...     return 2.0 * x
    >>> test_case3 = FunctionTestCase(
    ...     expected=0.7360212041 * u.kg,
    ...     function=double,
    ...     args=0.36801 * u.kg,
    ...     rtol=0.0001,  # relative tolerance; must be dimensionless
    ...     atol=0.000001 * u.kg,  # absolute tolerance; needs compatible units
    ... )
    >>> test_runner(test_case3)

    Next suppose we have a function that issues a `Warning` and returns
    a value.

    >>> import warnings
    >>> def issue_warning():
    ...     warnings.warn("warning message", Warning)
    ...     return 42

    We may set up a test case to check that ``issue_warning`` issues a
    `Warning` by setting ``expected=Warning``.

    >>> test_case4 = FunctionTestCase(function=issue_warning, expected=Warning)
    >>> test_runner(test_case4)

    We may set up a test case that checks both that a `Warning` is
    issued and the expected value is returned by having ``expected`` be
    a `tuple` containing a warning and the expected value (in either
    order).

    >>> test_case5 = FunctionTestCase(function=issue_warning, expected=(Warning, 42))
    >>> test_runner(test_case5)

    Similarly, we may test that a function raises an exception by
    setting ``expected`` to that exception.

    >>> def raise_exception():
    ...     raise Exception
    >>> test_case6 = FunctionTestCase(function=raise_exception, expected=Exception)
    >>> test_runner(test_case6)

    See Also
    --------
    ~plasmapy.tests.helpers.test_runner
    MethodTestCase
    AttrTestCase

    """

    def __init__(
        self,
        expected: Any,
        function: Callable,
        args=(),
        kwargs: Optional[Dict[str, Any]] = None,
        *,
        rtol: Real = 1e-8,
        atol=None,
    ):

        self._expected = expected
        self._function = function
        self._args = args
        self._kwargs = {} if kwargs is None else kwargs
        self._atol = atol
        self._rtol = rtol

    @property
    def expected(self):
        """The expected outcome of the test."""
        return self._expected

    @property
    def function(self) -> Callable:
        """The function to be tested."""
        return self._function

    @property
    def args(self):
        """The positional arguments to be provided to the function during the test."""
        return self._args

    @property
    def kwargs(self) -> Dict[str, Any]:
        """The keyword arguments to be provided to the function during the test."""
        return self._kwargs

    @property
    def rtol(self):
        """The relative tolerance for numerical comparisons."""
        return self._rtol

    @property
    def atol(self):
        """The absolute tolerance for numerical comparisons."""
        return self._atol


class MethodTestCase:
    """
    Stores the objects needed to test a method of a class instance.

    Parameters
    ----------
    expected
        The expected outcome of the test, which can be an exception,
        a warning, the resulting object, or a tuple that contains
        a warning and the resulting object.

    cls
        The class containing the method to be tested.

    method : str
        The name of the method to be tested.

    cls_args : optional, keyword-only
        The positional arguments to be provided to ``cls`` during
        instantiation.  If this is a `tuple` or `list`, then the
        arguments will be each of the items in the collection. If there
        is only one positional argument, then it may be inputted as is
        without being put into a `tuple` or `list`.

    cls_kwargs : `dict`, optional, keyword-only
        The keyword arguments to be provided to ``cls`` during
        instantiation. If provided, the keys of this `dict` must be
        strings.

    method_args : optional, keyword-only
        The positional arguments to be provided to the method during the
        test.  If this is a `tuple` or `list`, then the arguments will
        be each of the items in the collection. If there is only one
        positional argument, then it may be inputted as is without being
        put into a `tuple` or `list`.

    method_kwargs : `dict`, optional, keyword-only
        The keyword arguments to be provided to the method during the
        test. If provided, the keys of this `dict` must be
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

    Examples
    --------
    `MethodTestCase` stores the objects needed by
    `~plasmapy.tests.helpers.test_runner` to test that calling a method
    of a class instance results in the expected outcome (which could be
    a value, an exception, a warning, or a unit).

    Suppose we have a simple class with a simple method that returns a
    value.  Neither the class nor the method takes any arguments.

    >>> class SimpleClass:
    ...     def simple_method(self):
    ...         return 6

    We may set up a test case with `MethodTestCase` and run that case
    with `~plasmapy.tests.helpers.test_runner`.

    >>> from plasmapy.tests.helpers import MethodTestCase, test_runner
    >>> test_case = MethodTestCase(cls=SimpleClass, method="simple_method", expected=6)
    >>> test_runner(test_case)

    A more complex class may accept positional and/or keyword arguments
    upon instantiation, and methods in that class may also accept
    positional and/or keyword arguments.

    >>> class ComplexClass:
    ...     def __init__(self, cls_arg1, cls_arg2, *, cls_kwarg1=0, cls_kwarg2=0):
    ...         self.everything = [cls_arg1, cls_arg2, cls_kwarg1, cls_kwarg2]
    ...     def complex_method(self, method_arg, *, method_kwarg=0):
    ...         return sum(self.everything) + method_arg + method_kwarg

    We may use the ``cls_args`` and ``cls_kwargs`` arguments to specify
    the positional and keyword arguments to be supplied to the class
    during instantiation within the test, and ``method_args`` and
    ``method_kwargs`` to specify the arguments to be passed to the
    method.

    >>> more_complex_test_case = MethodTestCase(
    ...     expected=21,
    ...     cls=ComplexClass,
    ...     method="complex_method",
    ...     cls_args=(1, 2),
    ...     cls_kwargs={"cls_kwarg1": 3, "cls_kwarg2": 4},
    ...     method_args=5,
    ...     method_kwargs={"method_kwarg": 6},
    ... )
    >>> test_runner(more_complex_test_case)

    If either of ``cls_args`` or ``method_args`` is a `tuple` or `list`,
    then it will be treated as a collection of positional arguments.
    If either of these arguments is any other type of `object`, then it
    will be treated as the sole positional argument.

    See Also
    --------
    ~plasmapy.tests.helpers.test_runner
    FunctionTestCase
    AttrTestCase
    """

    def __init__(
        self,
        expected: Any,
        cls,
        method: str,
        *,
        cls_args=(),
        cls_kwargs: Optional[Dict[str, Any]] = None,
        method_args=(),
        method_kwargs: Optional[Dict[str, Any]] = None,
        atol=None,
        rtol=1e-8,
        purpose: Optional[str] = None,
    ):

        self._expected = expected
        self._cls = cls
        self._method = method
        self._cls_args = cls_args
        self._cls_kwargs = {} if cls_kwargs is None else cls_kwargs
        self._method_args = method_args
        self._method_kwargs = {} if method_kwargs is None else method_kwargs
        self._purpose = str(purpose)
        self._atol = atol
        self._rtol = rtol

    @property
    def expected(self):
        """The expected outcome of the test."""
        return self._expected

    @property
    def cls(self):
        """The class containing the method to be tested."""
        return self._cls

    @property
    def method(self) -> str:
        """The name of the method to be tested."""
        return self._method

    @property
    def cls_args(self):
        """
        The positional arguments to be provided to ``cls`` during
        instantiation.
        """
        return self._cls_args

    @property
    def cls_kwargs(self) -> Dict[str, Any]:
        """
        The keyword arguments to be provided to ``cls`` during
        instantiation.
        """
        return self._cls_kwargs

    @property
    def method_args(self):
        """
        The positional arguments to be provided to the method during the
        test.
        """
        return self._method_args

    @property
    def method_kwargs(self) -> Dict[str, Any]:
        """
        The keyword arguments to be provided to ``cls`` during
        instantiation.
        """
        return self._method_kwargs

    @property
    def rtol(self):
        """The relative tolerance for numerical comparisons."""
        return self._rtol

    @property
    def atol(self):
        """The absolute tolerance for numerical comparisons."""
        return self._atol


class AttrTestCase:
    """
    Stores the objects needed to test an attribute of a class instance.

    Parameters
    ----------
    expected
        The expected outcome of the test, which can be an exception,
        a warning, the resulting object, or a tuple that contains
        a warning and the resulting object.

    cls
        The class containing the attribute to be tested.

    attribute : str
        The name of the attribute to be tested.

    cls_args : optional, keyword-only
        The positional arguments to be provided to ``cls`` during
        instantiation.  If this is a `tuple` or `list`, then the
        arguments will be each of the items in the collection. If there
        is only one positional argument, then it may be inputted as is
        without being put into a `tuple` or `list`.

    cls_kwargs : `dict`, optional, keyword-only
        The keyword arguments to be provided to ``cls`` during
        instantiation. If provided, the keys of this `dict` must be
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


    See Also
    --------
    ~plasmapy.tests.helpers.test_runner
    FunctionTestCase
    MethodTestCase
    """

    def __init__(
        self,
        expected: Any,
        cls,
        attribute: str,
        *,
        cls_args=(),
        cls_kwargs: Optional[Dict[str, Any]] = None,
        atol=None,
        rtol=1e-8,
    ):

        self._expected = expected
        self._cls = cls
        self._attribute = attribute
        self._cls_args = cls_args
        self._cls_kwargs = {} if cls_kwargs is None else cls_kwargs
        self._atol = atol
        self._rtol = rtol

    @property
    def expected(self):
        """The expected outcome of the test."""
        return self._expected

    @property
    def cls(self):
        """The class containing the method to be tested."""
        return self._cls

    @property
    def attribute(self) -> str:
        """The name of the attribute to be tested."""
        return self._attribute

    @property
    def cls_args(self):
        """
        The positional arguments to be provided to ``cls`` during
        instantiation.
        """
        return self._cls_args

    @property
    def cls_kwargs(self) -> Dict[str, Any]:
        """
        The keyword arguments to be provided to ``cls`` during
        instantiation.
        """
        return self._cls_kwargs

    @property
    def rtol(self):
        """The relative tolerance for numerical comparisons."""
        return self._rtol

    @property
    def atol(self):
        """The absolute tolerance for numerical comparisons."""
        return self._atol
