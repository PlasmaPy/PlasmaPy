__all__ = [
    "RunTestError",
    "UnexpectedResultError",
    "InconsistentTypeError",
    "MissingExceptionError",
    "UnexpectedExceptionError",
    "UnexpectedResultError",
    "MissingExceptionError",
    "MissingWarningError",
    "IncorrectResultError",
    "InvalidTestError",
]


class RunTestError(Exception):
    """Base exception for test failures. Derived from `Exception`."""


class UnexpectedResultError(RunTestError):
    """
    Exception for when the actual result differs from the expected
    result.  Derived from `~plasmapy.pytest_helpers.utils.RunTestError`.
    """


class InconsistentTypeError(RunTestError):
    """
    Exception for when the type of the actual result differs from the
    type of the expected result.  Derived from
    `~plasmapy.utils.pytest_helpers.RunTestError`.
    """


class MissingExceptionError(RunTestError):
    """
    Exception for when an expected exception is not raised.  Derived
    from `~plasmapy.utils.pytest_helpers.RunTestError`.
    """


class UnexpectedExceptionError(RunTestError):
    """
    Exception for when an exception is expected, but a different
    exception is raised instead.  Derived from
    `~plasmapy.utils.pytest_helpers.RunTestError`.
    """


class MissingWarningError(RunTestError):
    """
    Exception for when a warning is expected to be issued, but isn't.
    Derived from `~plasmapy.utils.pytest_helpers.RunTestError`.
    """


class IncorrectResultError(RunTestError):
    """
    Exception for when the actual result differs from the expected
    result by more than the allowed tolerance.  Derived from
    `~plasmapy.utils.pytest_helpers.RunTestError`.
    """


class InvalidTestError(Exception):
    """
    Exception for when the inputs to a test are not valid.
    """
