import pytest

__all__ = [
    "Failed",
    "UnexpectedResultError",
    "InconsistentTypeError",
    "MissingExceptionError",
    "UnexpectedExceptionError",
    "MissingExceptionError",
    "UnexpectedWarningError",
    "MissingWarningError",
    "InvalidTestError",
    "ExceptionMismatchError",
    "WarningMismatchError",
]


# TODO: get the error message while raising this to be simply "Failed"
#       rather than plasmapy.tests.helpers.exceptions.Failed


class Failed(pytest.fail.Exception):
    """
    Base exception for test failures.

    Notes
    -----
    This exception was derived from `~pytest.fail.Exception`, which in
    turn was derived from `BaseException` (not `Exception`).
    """

    pass


class InconsistentTypeError(Failed):
    """
    Exception for when the type of the actual result differs from the
    type of the expected result.
    """

    pass


class MissingExceptionError(Failed):
    """
    Exception for when the expected exception is not raised.
    """

    pass


class UnexpectedExceptionError(Failed):
    """
    Exception for when an exception is raised unexpectedly.
    """

    pass


class ExceptionMismatchError(UnexpectedExceptionError, MissingExceptionError):
    """
    Exception for when an exception is expected, but a different
    exception is raised.
    """

    pass


class UnexpectedWarningError(Failed):
    """
    Exception for when a warning is issued unexpectedly.
    """

    pass


class MissingWarningError(Failed):
    """
    Exception for when an expected warning is not issued.
    """

    pass


class WarningMismatchError(UnexpectedWarningError, MissingWarningError):
    """
    Exception for when a warning is expected, but one or more
    different warnings were issued instead.
    """

    pass


class UnexpectedResultError(Failed):
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
