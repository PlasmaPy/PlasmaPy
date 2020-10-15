"""Exceptions that describe different types of test failures."""

__all__ = [
    "ExceptionMismatchFail",
    "InvalidTestError",
    "MissingExceptionFail",
    "MissingWarningFail",
    "TestFailure",
    "TypeMismatchFail",
    "UnexpectedExceptionFail",
    "UnexpectedResultFail",
    "UnexpectedWarningFail",
    "WarningMismatchFail",
]

import pytest


class TestFailure(pytest.fail.Exception):
    """
    Base exception for test failures.

    Notes
    -----
    This exception was derived from `~pytest.fail.Exception`, which in
    turn was derived from `BaseException` (not `Exception`).
    """

    pass


class TypeMismatchFail(TestFailure):
    """
    Exception for when the type of the actual result differs from the
    type of the expected result.
    """

    pass


class MissingExceptionFail(TestFailure):
    """
    Exception for when the expected exception is not raised.
    """

    pass


class UnexpectedExceptionFail(TestFailure):
    """
    Exception for when an exception is raised unexpectedly.
    """

    pass


class ExceptionMismatchFail(UnexpectedExceptionFail, MissingExceptionFail):
    """
    Exception for when an exception is expected, but a different
    exception is raised.
    """

    pass


class UnexpectedWarningFail(TestFailure):
    """
    Exception for when a warning is issued unexpectedly.
    """

    pass


class MissingWarningFail(TestFailure):
    """
    Exception for when an expected warning is not issued.
    """

    pass


class WarningMismatchFail(UnexpectedWarningFail, MissingWarningFail):
    """
    Exception for when a warning is expected, but one or more
    different warnings were issued instead.
    """

    pass


class UnexpectedResultFail(TestFailure):
    """
    Exception for when the actual result differs from the expected
    result by more than the allowed tolerance.
    """

    pass


class InvalidTestError(Exception):
    """
    Exception for when the inputs to a test are not valid.
    """

    pass
