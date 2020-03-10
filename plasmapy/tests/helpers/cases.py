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
    Contains the function to be tested, the arguments to be provided to
    the function during the test, and the test's expected outcome.

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


class MethodTestCase(_AbstractTestCase):
    """..."""

    def __init__(
        self,
        expected: Any,
        cls,
        method: str,
        *,
        cls_args=(),
        cls_kwargs: Optional[Dict[str, Any]],
        method_args=(),
        method_kwargs: Optional[Dict[str, Any]],
        atol=None,
        rtol=1e-8,
        purpose: Optional[str] = None,
    ):

        self.expected = expected
        self.cls = cls
        self.method = method
        self.cls_args = cls_args
        self.cls_kwargs = {} if cls_kwargs is None else cls_kwargs
        self.method_args = method_args
        self.method_kwargs = method_kwargs
        self.purpose = str(purpose)
        self.atol = atol
        self.rtol = rtol


class AttrTestCase(_AbstractTestCase):
    """..."""

    def __init__(
        self,
        expected: Any,
        cls,
        attribute: str,
        *,
        cls_args=(),
        cls_kwargs: Optional[Dict[str, Any]] = None,
        purpose: Optional[str] = None,
        atol=None,
        rtol=1e-8,
    ):

        self.expected = expected
        self.cls = cls
        self.attribute = attribute
        self.cls_args = cls_args
        self.cls_kwargs = {} if cls_kwargs is None else cls_kwargs
        self.purpose = purpose
        self.atol = atol
        self.rtol = rtol
