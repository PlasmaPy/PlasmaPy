"""Exceptions that describe different types of test failures."""

__all__ = [
    "ExceptionMismatchFail",
    "InvalidTestError",
    "MissingExceptionFail",
    "MissingWarningFail",
    "TestFailed",
    "TypeMismatchFail",
    "UnexpectedExceptionFail",
    "UnexpectedResultFail",
    "UnexpectedWarningFail",
    "WarningMismatchFail",
]

from _pytest.outcomes import Failed


class TestFailed(Failed):
    """
    Base exception for test failures.

    Notes
    -----
    This exception was derived from ``pytest.fail.Exception``, which in
    turn was derived from `BaseException` (not `Exception`).
    """


class TypeMismatchFail(TestFailed):
    """
    Exception for when the type of the actual result differs from the
    type of the expected result.
    """


class MissingExceptionFail(TestFailed):
    """
    Exception for when the expected exception is not raised.
    """


class UnexpectedExceptionFail(TestFailed):
    """
    Exception for when an exception is raised unexpectedly.
    """


class ExceptionMismatchFail(UnexpectedExceptionFail, MissingExceptionFail):
    """
    Exception for when an exception is expected, but a different
    exception is raised.
    """


class UnexpectedWarningFail(TestFailed):
    """
    Exception for when a warning is issued unexpectedly.
    """


class MissingWarningFail(TestFailed):
    """
    Exception for when an expected warning is not issued.
    """


class WarningMismatchFail(UnexpectedWarningFail, MissingWarningFail):
    """
    Exception for when a warning is expected, but one or more
    different warnings were issued instead.
    """


class UnexpectedResultFail(TestFailed):
    """
    Exception for when the returned value does not match the value that
    was expected.
    """


class InvalidTestError(Exception):
    """
    Exception for when the inputs to a test are not valid.
    """
