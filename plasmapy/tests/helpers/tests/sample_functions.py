"""Sample functions and classes to be used in tests."""

import warnings
from typing import NoReturn

import astropy.units as u
import numpy as np


class SampleException(Exception):
    """A sample exception to be used for testing purposes."""

    pass


class SampleExceptionSubclass(SampleException):
    """
    The subclass of the sample exception to used for testing purposes.

    If the `~pytest.raises` context manager expects a certain exception
    to be raised, then the test will pass if a subclass of that exception
    is raised.  This subclass is used to test that PlasmaPy's test
    helpers functionality catches situations like that.
    """

    pass


class SampleWarning(Warning):
    """A sample warning to be used for testing purposes."""

    pass


class SampleWarningSubclass(SampleWarning):
    """
    The subclass of the sample warning to be used for testing purposes.

    If the `~pytest.warns` context manager expects a certain warning to
    be issued, then the test will pass if a subclass of that warning is
    issued.  This subclass is used to test that PlasmaPy's test helpers
    functionality catches situations like that.
    """

    pass


def return_42() -> int:
    """A sample function for testing purposes that returns ``42``."""

    return 42


def return_42_meters() -> u.Quantity:
    """
    A sample function for testing purposes that returns an
    `~astropy.units.Quantity` with a value of ``42.0 m``.
    """

    return 42.0 * u.m


def return_np_array(*args) -> np.array:
    """A function to be used when testing `~numpy.array` instances."""

    return np.array(args)


def issue_warning() -> NoReturn:
    """A sample function that issues a `SampleWarning`."""

    warnings.warn("warning message", SampleWarning)


def issue_warning_return_42() -> int:
    """
    A sample function for testing purposes that issues a `SampleWarning`
    and returns ``42``.
    """

    warnings.warn("warning message", SampleWarning)
    return 42


def raise_exception():
    """A sample function for testing purposes that raises a `SampleException`."""

    raise SampleException("exception message")


def sum_of_args_and_kwargs(arg1, arg2, *, kw1=None, kw2=None):
    """
    A sample function for testing purposes that returns the sum of
    two positional arguments and two keyword arguments.
    """

    return arg1 + arg2 + kw1 + kw2


def return_none() -> None:
    """A sample function for testing purposes that returns `None`."""

    return None


class SampleClass1:
    """A sample class to be used for testing purposes."""

    def __init__(self, *args, **kwargs):

        pass

    @classmethod
    def arg_plus_kwarg(self, arg, *, kwarg=None):
        """
        A sample method that returns the sum of a positional argument
        and a keyword argument.
        """

        return arg + kwarg

    @property
    def forty(self) -> int:
        """
        A sample attribute that returns ``40``.
        """

        return 40

    def raise_exception(self):
        """A sample method that raises a `SampleException`."""

        raise SampleException("error message")

    def issue_warning(self) -> NoReturn:
        """A sample method that issues a `SampleWarning`."""

        warnings.warn("warning message", SampleWarning)


class SampleClass2:
    """A sample class to be used for testing purposes."""

    def __init__(self, cls_arg1, cls_arg2, *, cls_kwarg1=None, cls_kwarg2=None):

        self.cls_arg1 = cls_arg1
        self.cls_arg2 = cls_arg2
        self.cls_kwarg1 = cls_kwarg1
        self.cls_kwarg2 = cls_kwarg2

    def method(
        self, method_arg1, method_arg2, *, method_kwarg1=None, method_kwarg2=None
    ):
        """
        Return the sum of the positional and keyword arguments supplied
        to the class upon instantiation plus the sum of the positional
        and keyword arguments supplied to the method when it is called.
        """

        return sum(
            [
                self.cls_arg1,
                self.cls_arg2,
                self.cls_kwarg1,
                self.cls_kwarg2,
                method_arg1,
                method_arg2,
                method_kwarg1,
                method_kwarg2,
            ]
        )

    @property
    def attr(self):
        """
        Return the sum of the positional and keyword arguments supplied
        to the class upon instantiation.
        """

        return self.cls_arg1 + self.cls_arg2 + self.cls_kwarg1 + self.cls_kwarg2
