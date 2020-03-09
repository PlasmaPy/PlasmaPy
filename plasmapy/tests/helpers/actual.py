"""Classes that represent the actual outcome of running a test."""

__all__ = ["ActualTestOutcome"]

import warnings
from typing import Any, List, Union

import pytest
from plasmapy.tests.helpers.exceptions import InvalidTestError
from plasmapy.tests.helpers.inputs import (
    AbstractTestInputs,
    ClassAttributeTestInputs,
    ClassMethodTestInputs,
    FunctionTestInputs,
)


class _ExitPytestRaises(Exception):
    """An `Exception` to be raised to exit `pytest.raises` context manager."""

    pass


class _ExitPytestWarns(Warning):
    """A `Warning` to be raised to exit `pytest.warns` context manager."""

    pass


class ActualTestOutcome:
    def __init__(
        self,
        inputs: Union[
            FunctionTestInputs, ClassMethodTestInputs, ClassAttributeTestInputs
        ],
    ):
        """
        A class to record the actual outcome of a test.

        Parameters
        ----------
        inputs : `plasmapy.utils.pytest_helpers.inputs.AbstractTestInputs`
            The instance of a subclass of
            `plasmapy.utils.pytest_helpers.inputs.AbstractTestInputs`
            that contains the function of class to be tested, along with
            corresponding positional and keyword arguments.
        """

        self._inputs = inputs

        if not isinstance(inputs, AbstractTestInputs):
            raise TypeError(
                "Expecting an instance of a subclass of AbstractTestInputs,"
                "such as FunctionTestInputs, ClassAttributeTestInputs, or "
                "ClassMethodTestInputs."
            )

        with pytest.warns(Warning) as warnings_record:
            warnings.warn(
                "So we can exit pytest.warns context manager", _ExitPytestWarns
            )
            with pytest.raises(BaseException) as exception_info:
                self._value = inputs.call()
                raise _ExitPytestRaises("So we can exit pytest.raises context manager")

        self._exception_info = exception_info

        # Remove _MockWarning from the list of warnings before storing it
        warnings_record.pop()
        self._warnings_record = warnings_record

    @property
    def value_was_returned(self) -> bool:
        """
        `True` if the `callable` being tested returned a value or
        `None`, and `False` otherwise.
        """

        return hasattr(self, "_value")

    @property
    def exception_was_raised(self) -> bool:
        """
        `True` if an exception was raised by the `callable` being
        tested, and `False` otherwise.
        """

        return self._exception_info.type is not _ExitPytestRaises

    @property
    def warning_was_issued(self) -> bool:
        """
        `True` if a warning was issued by the `callable` being tested,
        and `False` otherwise.
        """

        return bool(self._warnings_record.list)

    @property
    def exception_info(self):
        """
        The ``ExceptionInfo`` instance created while using the `pytest.raises`
        context manager.

        Raises
        ------
        InvalidTestError
            If no exception was raised during the test.
        """

        if self.exception_was_raised:
            return self._exception_info
        else:
            raise InvalidTestError("No exception was raised.")

    @property
    def exception_type(self):
        """
        Return the type of exception that was raised.

        Raises
        ------
        InvalidTestError
            If no exception was raised during the test.
        """

        return self.exception_info.type

    @property
    def exception_message(self) -> str:
        """
        The error message of the exception that was raised during the
        test.

        Raises
        ------
        InvalidTestError
            If no exception was raised during the test.
        """

        return str(self.exception_info.value)

    @property
    def warnings_record(self):
        """
        The ``WarningsRecorder`` instance created by the `pytest.warns`
        context manager.

        Raises
        ------
        InvalidTestError
            If no warning was issued during the test.
        """

        if self.warning_was_issued:
            return self._warnings_record
        else:
            raise InvalidTestError(
                f"The warnings_record attribute is not available because "
                f" no warnings were issued when running the command: "
                f"{self._inputs.call_string}"
            )

    @property
    def warning_types(self) -> List[Warning]:
        """
        A `tuple` containing the warnings that were issued during the
        test, corresponding to the ``warning_messages`` attribute.

        Raises
        ------
        InvalidTestError
            If no warning was issued during the test.
        """

        if self.warning_was_issued:
            return [warning.category for warning in self.warnings_record]
        else:
            raise InvalidTestError(
                f"The warning_types attribute is not available because "
                f"no warnings were issued when running the command: "
                f"{self._inputs.call_string}"
            )

    @property
    def warning_messages(self) -> List[str]:
        """
        A `tuple` containing the warning messages that were issued
        during the test, corresponding to the ``warning_types``
        attribute.

        Raises
        ------
        InvalidTestError
            If no warning was issued during the test.
        """

        if self.warning_was_issued:
            return [str(warning.message.args[0]) for warning in self.warnings_record]
        else:
            raise InvalidTestError(
                f"The warning_messages attribute is not available because "
                f"no warnings were issued when running the command: "
                f"{self._inputs.call_string}"
            )

    @property
    def value(self) -> Any:
        """
        The value returned by the function, class, or class method that
        is being tested.

        Raises
        ------
        InvalidTestError
            If no value was returned during the test, for example if
            the test raised an exception.
        """

        if self.value_was_returned:
            return self._value
        else:
            raise InvalidTestError("No value was returned.")

    @property
    def call_string(self) -> str:
        """A string that reproduces the call that is being tested."""

        return self._inputs.call_string

    def __len__(self):
        """
        Return the length of the actual value returned in the test, if a
        value was returned.  If ``__len__`` is undefined in the actual
        value, then return ``1``.
        """

        if self.value_was_returned:
            has_a_len = hasattr(self.value, "__len__")
            return len(self.value) if has_a_len else 1
        else:
            return NotImplemented
