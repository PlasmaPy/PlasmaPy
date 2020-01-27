"""Sample functions and classes to be used in tests."""

import warnings
import numpy as np
import astropy.units as u


class SampleException(Exception):
    pass


class SampleExceptionSubclass(SampleException):
    pass


class SampleWarning(Warning):
    pass


class SampleWarningSubclass(SampleWarning):
    pass


def return_42() -> int:
    return 42


def return_42_meters() -> u.Quantity:
    return 42.0 * u.m


def return_np_array(*args) -> np.array:
    return np.array(args)


def issue_warning_return_42() -> int:
    warnings.warn("warning message", SampleWarning)
    return 42


def raise_exception():
    raise SampleException("exception message")


def sum_of_args_and_kwargs(arg1, arg2, *, kw1=None, kw2=None):
    return arg1 + arg2 + kw1 + kw2


def return_none():
    return None


class SampleClass:
    def __init__(self, *args, **kwargs):
        pass

    def arg_plus_kwarg(self, arg, *, kwarg=None):
        return arg + kwarg

    @property
    def forty(self) -> int:
        return 40
