"""Contains test comparators: tools that automate comparisons."""

from plasmapy.utils.pytest_helpers.expected import ExpectedTestOutcome
from plasmapy.utils.pytest_helpers.actual import ActualTestOutcome

from plasmapy.utils.pytest_helpers.error_messages import (
    _exc_str,
    _string_together_warnings_for_printing,
)

__all__ = ['CompareActualExpected']

# TODO: Test that error messages are working correctly!  This is going to
#       be hard to manually check.


class CompareActualExpected:

    def __init__(self, actual: ActualTestOutcome, expected: ExpectedTestOutcome):

        if not isinstance(actual, ActualTestOutcome):
            raise TypeError("Expecting an instance of ActualTestOutcome")
        if not isinstance(expected, ExpectedTestOutcome):
            raise TypeError("Expecting an instance of ExpectedTestOutcome")

        self._actual = actual
        self._expected = expected

        self._make_errmsg_if_necessary()

    @property
    def actual(self) -> ActualTestOutcome:
        """
        The instance of `~plasmapy.utils.pytest_helpers.actual.ActualTestOutcome`
        that is being compared.
        """

        return self._actual

    @property
    def expected(self) -> ExpectedTestOutcome:
        """
        The instance of `~plasmapy.utils.pytest_helpers.expected.ExpectedTestOutcome`
        that is being compared.
        """

        return self._expected

    @property
    def test_passed(self) -> bool:
        """
        Returns `True` if the actual outcome matches the expected outcome,
        and `False` otherwise.
        """

        return not bool(self.error_message)

    @property
    def error_message(self) -> str:
        """
        If the test failed, then return an appropriate error message.
        If the test passed, then return an empty string.
        """

        return ' '.join(self._error_messages_list)

    def _add_errmsg(self, errmsg):
        """Add an error message to the list of error messages."""
        self._error_messages_list.append(errmsg)

    def _make_unexpected_warnings_errmsg(self):
        """
        Compose an error message for tests where warnings were
        unexpectedly issued.
        """

        _string_together_warnings_for_printing(
            self.actual.warning_types,
            self.actual.warning_messages,
        )

        is_first_error = not self._error_messages_list
        subject = f"The command {self.actual.call_string}" if is_first_error else "This command"

        warnings_for_printing = _string_together_warnings_for_printing(
            self.actual.warning_types,
            self.actual.warning_messages,
        )

        number_of_warnings = len(self.actual.warning_types)

        errmsg = (
            f"{subject} unexpectedly issued the following warnings"
            f"{'s' if number_of_warnings > 1 else ''}:"
            f"\n\n"
            f"{warnings_for_printing}"
        )

        self._add_errmsg(errmsg)

    def _make_exception_mismatch_errmsg_if_necessary(self):
        """
        Compose an error message for tests where a certain type of
        exception was raised, but a different type of exception was
        expected.  If the expected exception was actually raised, do
        nothing.
        """

        expected_exception = self.expected.expected_exception
        actual_exception = self.actual.exception_type

        if actual_exception is not expected_exception:

            self._add_errmsg(
                f"The command {self.actual.call_string} raised "
                f"{_exc_str(self.actual.exception_type)}, instead of "
                f"{_exc_str(self.expected.expected_exception)} as expected."
            )

        # TODO: Create a different error message for cases where the
        #       expected error message is a subclass of the actual
        #       error message, or vice versa.

    def _make_missing_exception_errmsg(self):
        """
        Compose an error message for tests where an exception should
        have been raised, but was not.
        """

        self._add_errmsg(
            f"The command {self.actual.call_string} did not raise "
            f"{_exc_str(self.expected.expected_exception)} as expected. "
            f"Instead, this command returned the unexpected value of "
            f"{repr(self.actual.value)}."
        )

        # TODO: improve representation of the value (as repr doesn't
        #       always result in something particularly readable.)

    def _make_unexpected_exception_errmsg(self):
        """
        Compose an error message for tests where an exception was
        unexpectedly raised.
        """

        self._add_errmsg(
            f"The command {self.actual.call_string} unexpectedly raised "
            f"{_exc_str(self.actual.exception_type)}."
        )

    def _make_value_mismatch_errmsg_if_necessary(self):
        """
        Compose an error message for tests where the expected value does
        not match the value that was actually returned.  If the expected
        and actual values match, then do nothing.
        """

        # TODO: deal with all of the weird comparison cases (like making
        #       sure that a Quantity has a unit, a Quantity matches another
        #       Quantity, arrays match, etc.

        expected_value = self.expected.expected_value
        actual_value = self.actual.value

        if actual_value is not expected_value and actual_value != expected_value:
            self._add_errmsg('_make_value_mismatch_errmsg_if_necessary')

    def _make_missing_warning_errmsg(self):
        """
        Compose an error message for tests where a warning was expected
        to be issued, but was not.
        """

        is_first_error = not self._error_messages_list
        subject = f"The command {self.actual.call_string}" if is_first_error else "This command "

        self._add_errmsg(
            f"{subject} did not raise {_exc_str(self.expected.expected_warning)}"
            f"as expected."
        )

    def _make_warning_mismatch_errmsg_if_necessary(self):
        """
        Compose an error message for tests where the expected warnings
        were not issued but the
        """

        expected_warning = self.expected.expected_warning
        actual_warnings = self.actual.warning_types
        warning_messages = self.actual.warning_messages

        number_of_warnings = len(actual_warnings)

        is_first_error = bool(self._error_messages_list)

        if expected_warning in actual_warnings:
            return

        subject = f"The command {self.actual.call_string}" if is_first_error else "This command"

        errmsg = (
            f"{subject} was expected to issue a {_exc_str(expected_warning)}, "
            f"but instead issued the following warning"
            f"{'' if number_of_warnings == 1 else 's'}:"
        )
        errmsg += "\n\n"

        for warning, message in zip(actual_warnings, warning_messages):
            errmsg += warning.__name__ + message + "\n\n"

        self._add_errmsg(errmsg)

        # TODO: Figure out a way to deal to deal with deprecation warnings.
        #       We should not count those as test failures, but those should
        #       show up in the test report as an actual warning.

    def _make_errmsg_if_necessary(self):
        """
        Determine whether or not the test passes or fails.  If the test
        fails, create a descriptive test failure error message.  If the
        test passes, the error message that will be generated will be an
        empty string.
        """

        self._error_messages_list = []

        expecting_an_exception = self.expected.expecting_an_exception
        expecting_a_value = self.expected.expecting_a_value
        expecting_a_warning = self.expected.expecting_a_warning

        exception_was_raised = self.actual.exception_was_raised
        value_was_returned = self.actual.value_was_returned
        warning_was_issued = self.actual.warning_was_issued

        if expecting_an_exception and exception_was_raised:
            self._make_exception_mismatch_errmsg_if_necessary()

        if expecting_an_exception and not exception_was_raised:
            self._make_missing_exception_errmsg()

        if not expecting_an_exception and exception_was_raised:
            self._make_unexpected_exception_errmsg()

        if expecting_a_value and value_was_returned:
            self._make_value_mismatch_errmsg_if_necessary()

        if expecting_a_warning and warning_was_issued:
            self._make_warning_mismatch_errmsg_if_necessary()

        if expecting_a_warning and not warning_was_issued:
            self._make_missing_warning_errmsg()

        if not expecting_a_warning and warning_was_issued:
            self._make_unexpected_warnings_errmsg()


# TODO: Create a class to represent actual outcomes with other actual
#       outcomes to make sure that they're consistent.  That would be
#       helpful for cases like the ones run_test_equivalent_calls are
#       covering.  In particular, there are occasional symmetry
#       properties that need to be tested.  This class will be used much
#       less often than CompareActualExpected, but could save some work
#       in the future when we have to write such tests.
