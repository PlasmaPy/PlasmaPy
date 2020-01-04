import pytest
import warnings
from typing import List, Any
from plasmapy.utils.pytest_helpers.inputs import AbstractTestInputs

__all__ = ["ActualTestOutcome"]


class _MockException(Exception):
    pass


class _MockWarning(Warning):
    pass


class ActualTestOutcome:
    def __init__(self, inputs):
        """
        A class to record the actual outcome of a test.

        Parameters
        ----------
        inputs : instance of subclass of `plasmapy.utils.pytest_helpers.inputs.AbstractTestInputs`

        """

        self._inputs = inputs

        if not isinstance(inputs, AbstractTestInputs):
            raise TypeError(
                "Expecting an instance of a subclass of AbstractTestInputs,"
                "such as FunctionTestInputs, ClassAttributeTestInputs, or "
                "ClassMethodTestInputs.")

        with pytest.warns(Warning) as warnings_record:
            warnings.warn('So we can exit pytest.warns context manager', _MockWarning)
            with pytest.raises(Exception) as exception_info:
                self._value = inputs.call()
                raise _MockException('So we can exit pytest.raises context manager')

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
        return hasattr(self, '_value')

    @property
    def exception_was_raised(self) -> bool:
        """
        `True` if an exception was raised by the `callable` being
        tested, and `False` otherwise.
        """
        return self._exception_info.type is not _MockException

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
        RuntimeError
            If no exception was raised during the test.

        """
        if self.exception_was_raised:
            return self._exception_info
        else:
            raise RuntimeError("No exception was raised.")

    @property
    def exception_type(self):
        """
        Return the type of exception that was raised.

        Raises
        ------
        RuntimeError
            If no exception was raised during the test.

        """
        return self.exception_info.type

    @property
    def exception_message(self):
        """
        The error message of the exception that was raised during the
        test.

        Raises
        ------
        RuntimeError
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
        RuntimeError
            If no warning was issued during the test.

        """
        if self.warning_was_issued:
            return self._warnings_record
        else:
            raise RuntimeError("Warnings information is not available")

    @property
    def warning_types(self) -> List:
        """
        A `tuple` containing the warnings that were issued during the
        test, corresponding to the ``warning_messages`` attribute.

        Raises
        ------
        RuntimeError
            If no warning was issued during the test.
        """
        return [warning.category for warning in self.warnings_record]

    @property
    def warning_messages(self) -> List[str]:
        """
        A `tuple` containing the warning messages that were issued
        during the test, corresponding to the ``warning_types``
        attribute.

        Raises
        ------
        RuntimeError
            If no warning was issued during the test.

        """
        return [str(warning.message.args[0]) for warning in self.warnings_record]

    @property
    def value(self) -> Any:
        """
        The value returned by the function, class, or class method that
        is being tested.

        Raises
        ------
        RuntimeError
            If no value was returned during the test, for example if
            the test raised an exception.

        """
        if self.value_was_returned:
            return self._value
        else:
            raise RuntimeError("No value was returned.")

    @property
    def call_string(self) -> str:
        """A string that reproduces the call that is being tested."""
        return self._inputs.call_string
