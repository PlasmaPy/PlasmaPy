"""Exceptions that describe different types of test failures."""

__all__ = [
    "ExceptionMismatch",
    "TestFailure",
    "InconsistentType",
    "InvalidTest",
    "MissingException",
    "MissingWarning",
    "UnexpectedException",
    "UnexpectedResult",
    "UnexpectedWarning",
    "WarningMismatch",
]

import pytest

# TODO: get the error message while raising this to be simply "TestFailure"
#       rather than plasmapy.tests.helpers.exceptions.TestFailure


class TestFailure(pytest.fail.Exception):
    """
    Base exception for test failures.

    Notes
    -----
    This exception was derived from `~pytest.fail.Exception`, which in
    turn was derived from `BaseException` (not `Exception`).
    """

    pass


class TypeMismatch(TestFailure):
    """
    Exception for when the type of the actual result differs from the
    type of the expected result.
    """

    pass


class MissingException(TestFailure):
    """
    Exception for when the expected exception is not raised.
    """

    pass


class UnexpectedException(TestFailure):
    """
    Exception for when an exception is raised unexpectedly.
    """

    pass


class ExceptionMismatch(UnexpectedException, MissingException):
    """
    Exception for when an exception is expected, but a different
    exception is raised.
    """

    pass


class UnexpectedWarning(TestFailure):
    """
    Exception for when a warning is issued unexpectedly.
    """

    pass


class MissingWarning(TestFailure):
    """
    Exception for when an expected warning is not issued.
    """

    pass


class WarningMismatch(UnexpectedWarning, MissingWarning):
    """
    Exception for when a warning is expected, but one or more
    different warnings were issued instead.
    """

    pass


class UnexpectedResult(TestFailure):
    """
    Exception for when the actual result differs from the expected
    result by more than the allowed tolerance.
    """

    pass


class InvalidTest(Exception):
    """
    Exception for when the inputs to a test are not valid.
    """

    pass
