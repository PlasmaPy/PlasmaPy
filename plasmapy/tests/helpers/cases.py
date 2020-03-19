"""Classes to contain all inputs for a test case and the expected outcome."""

__all__ = ["FunctionTestCase", "MethodTestCase", "AttrTestCase"]

from abc import ABC, abstractmethod
from numbers import Real
from typing import Any, Optional, Dict, Callable


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
    Stores the objects needed to test a function.

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
    >>> return_arg = lambda arg: arg
    >>> test_case = FunctionTestCase(function=return_arg, args=42, expected=42)
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
