"""Contains test comparators: tools that automate comparisons."""

from numbers import Number
from typing import Union, Tuple, Any, Optional

from astropy import units as u

from plasmapy.utils.pytest_helpers import InvalidTestError

from plasmapy.utils.pytest_helpers.expected import ExpectedTestOutcome
from plasmapy.utils.pytest_helpers.actual import ActualTestOutcome

from plasmapy.utils.pytest_helpers.error_messages import (
    _exc_str,
    _string_together_warnings_for_printing,
)


__all__ = ["CompareActualExpected"]


# TODO: Test that error messages are working correctly!  This is going to
#       be hard to manually check.

class CompareValues:
    """
    Compares the properties of two values.

    Parameters
    ----------
    this
        The first `object` to be compared.

    that
        The second `object` to be compared.

    rtol : number or dimensionless `~astropy.units.Quantity`, optional, keyword-only
        The relative tolerance to be supplied to `~astropy.units.isclose`
        or `~astropy.units.allclose`.  Defaults to ``1e-6``.

    atol : number or `~astropy.units.Quantity`, optional, keyword-only
        The absolute tolerance to be supplied to `astropy.units.isclose`
        or `~astropy.units.allclose`.  If ``atol`` is a
        `~astropy.units.Quantity`, then it must have the same units as
        ``this`` and ``that``.  Defaults to zero in the appropriate units.
    """

    def __init__(
            self,
            this: Any,
            that: Any,
            *,
            rtol: Union[Number, u.Quantity] = 1e-6,
            atol: Optional[Union[Number, u.Quantity]] = None,
    ):

        self._this = this
        self._that = that

        try:
            self.rtol = rtol
            self.atol = atol
        except Exception as exc:  # TODO: make more specific
            raise InvalidTestError(
                "Cannot perform comparison with invalid values for the "
                "relative and/or absolute tolerances for comparison."
            ) from exc

    @property
    def values(self) -> Tuple[Any, Any]:
        """
        The two values that are being compared.
        """

        return (self._this, self._that)

    @property
    def types(self) -> Tuple[type, type]:
        """
        The types of the two values that are being compared.
        """

        return (type(self._this), type(self._that))

    @property
    def rtol(self) -> Union[Number, u.Quantity]:
        """
        The relative tolerance to be used by `~astropy.units.isclose` or
        `~astropy.units.allclose`.

        If ``rtol`` is a `~astropy.units.Quantity`, then it must be
        dimensionless.
        """

        return self._rtol

    @property
    def atol(self) -> Optional[Union[Number, u.Quantity]]:
        """
        The absolute tolerance to be used by `~astropy.units.isclose` or
        `~astropy.units.allclose`.

        If the objects being compared are two `~astropy.units.Quantity`
        instances, then ``atol`` must have consistent dimensions.
        """

        return self._atol

    @rtol.setter
    def rtol(self, relative_tolerance):
        self._rtol = relative_tolerance
        raise NotImplementedError

    @atol.setter
    def atol(self, absolute_tolerance):
        self._atol = absolute_tolerance
        raise NotImplementedError


    @property
    def are_identical(self) -> bool:
        """
        Return `True` if the two values refer to the same identical
        object, and `False` otherwise.
        """

        return self._this is self._that

    @property
    def are_equal(self) -> bool:
        """
        Return `True` if the two values evaluate as equal to each other,
        and `False` otherwise (including if the .
        """

        try:
            return self._this == self._that
        except Exception:
            return False

    @property
    def have_same_types(self) -> bool:
        """
        Return `True` if the two values are of the same `type`, and
        `False` otherwise.
        """

        return self.types[0] is self.types[1]

    @property
    def are_quantities(self) -> bool:
        """
        `True` if each `object` being compared is a `Quantity`, and
        `False` otherwise.
        """

        return isinstance(self._this, u.Quantity) and isinstance(self._that, u.Quantity)

    @property
    def are_quantity_and_unit(self) -> bool:
        """
        `True` if the first `object` being compared is a `~astropy.units.Quantity`
        instance and the second `object` being compared is an instance of
        (a subclass of) `~astropy.units.UnitBase`, and `False` otherwise.
        """

        return isinstance(self._this, u.Quantity) and isinstance(self._that, u.UnitBase)

    @property
    def units_are_consistent(self) -> bool:
        """
        If the objects being compared are both `~astropy.units.Quantity`
        instances, then return `True` if the units are identical and
        `False` otherwise.

        If the first object is a `~astropy.units.Quantity` and the second
        object is a `~astropy.units.UnitBase`

        Return `True` if the two objects being compared are both
        `~astropy.units.Quantity` instances with identical units, or if
        the first object is a `~Quantity.`

        This attribute will raise a `u.UnitsError` if the values being
        compared are a `
        """

        if self.are_quantities:
            return self._this.unit == self._that.unit
        elif self.are_quantity_and_unit:
            return self._this.unit == self._that
        else:
            raise u.UnitsError(
                "The `units_are_consistent` attribute of CompareValues "
                "may only be used when comparing a Quantity with a Quantity, "
                "or a Quantity with a Unit (in that order)."
            )


    @property
    def are_close(self):
        """
        `True` if
        """
        return u.allclose(self._this, self._that, rtol=self.rtol, atol=self.atol)


    def __bool__(self):
        if self.are_identical:
            return True
        elif self.are_equal and self.have_same_types:
            return True
        elif self.are_quantity_and_unit:
            return


class CompareActualExpected:
    """
    A class to compare the actual and expected results of a test, and
    generate an appropriate error message if necessary.


    """


    def __init__(
        self,
        actual: ActualTestOutcome,
        expected: ExpectedTestOutcome,
        *,
        rtol=1e-9,
        atol=0.0,
    ):

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
    def rtol(self) -> Union[Number, u.Quantity]:
        """
        The relative tolerance to be used in `~astropy.units.isclose` or
        `astropy.units.allclose`.  If ``rtol`` is a `~astropy.units.Quantity`
        instance, then it must be dimensionless.
        """
        return self._rtol

    @property
    def atol(self) -> Union[Number, u.Quantity]:
        """
        The absolute tolerance to be used in `~astropy.units.isclose`
        or `~astropy.units.allclose`.
        """
        return self._atol

    @rtol.setter
    def rtol(self, relative_tolerance):
        if isinstance(relative_tolerance, u.Quantity):
            if relative_tolerance.unit == u.dimensionless_unscaled:
                self._rtol = relative_tolerance
            else:
                raise u.UnitsError("`rtol` should be dimensionless if it is a Quantity")
        elif isinstance(relative_tolerance, Number):
            self._rtol = relative_tolerance
        else:
            raise TypeError(f"Invalid value of `rtol`: {relative_tolerance}")

    @atol.setter
    def atol(self, absolute_tolerance):

        # if it is a quantity, then it must have the same units as the expected value,
        # and otherwise it should raise an InvalidTestError
        #
        raise NotImplementedError

#        is_a_quantity = isinstance(absolute_tolerance, u.Quantity)
#        expecting_a_value = self.expected.expecting_a_value
#        expecting_a_quantity = isinstance()

#        if expecting_a_value:


#        expecting_a_quantity = self.expected.expecting_a_value and self.ex

#        expecting_a_quantity = isinstance(self.expected.expecting_a_value)

#        if isinstance(u.Quantity):


#            if isinstance(self.expected.expected_value)

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

        return " ".join(self._error_messages_list)

    def _add_errmsg(self, errmsg):
        """Add an error message to the list of error messages."""
        self._error_messages_list.append(errmsg)

    @property
    def _subject(self):
        """
        Return an appropriate subject for the first sentence in the
        error message.  The call string should be included in each error
        message if and only if it is the first error message.
        """
        is_first_error = not self._error_messages_list
        return (
            f"The command {self.actual.call_string}"
            if is_first_error
            else "This command"
        )


    def _make_unexpected_warnings_errmsg(self):
        """
        Compose an error message for tests where warnings were
        unexpectedly issued.
        """

        _string_together_warnings_for_printing(
            self.actual.warning_types, self.actual.warning_messages,
        )

        is_first_error = not self._error_messages_list
        subject = (
            f"The command {self.actual.call_string}"
            if is_first_error
            else "This command"
        )

        warnings_for_printing = _string_together_warnings_for_printing(
            self.actual.warning_types, self.actual.warning_messages,
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
        both_values = (expected_value, actual_value)

        are_same_object = expected_value is actual_value
        test_is_unit_check = isinstance(expected_value, u.Unit)

        both_are_quantities = all([
            isinstance(expected_value, u.Quantity),
            isinstance(actual_value, u.Quantity),
        ])

        try:
            are_equal = expected_value == actual_value
        except Exception:
            are_equal = False
        finally:
            same_type = type(expected_value) is type(actual_value)

        if are_equal and same_type:
            return






        if actual_value is not expected_value and actual_value != expected_value:
            self._add_errmsg("_make_value_mismatch_errmsg_if_necessary")

    def _make_missing_warning_errmsg(self):
        """
        Compose an error message for tests where a warning was expected
        to be issued, but was not.
        """

        is_first_error = not self._error_messages_list
        subject = (
            f"The command {self.actual.call_string}"
            if is_first_error
            else "This command "
        )

        self._add_errmsg(
            f"{subject} did not raise {_exc_str(self.expected.expected_warning)}"
            f"as expected."
        )

    def _make_warning_mismatch_errmsg_if_necessary(self):
        """
        Compose an error message for tests where the expected warning
        was not issued but different warning(s) were.
        """

        expected_warning = self.expected.expected_warning
        actual_warnings = self.actual.warning_types
        warning_messages = self.actual.warning_messages

        number_of_warnings = len(actual_warnings)

        is_first_error = bool(self._error_messages_list)

        if expected_warning in actual_warnings:
            return

        subject = (
            f"The command {self.actual.call_string}"
            if is_first_error
            else "This command"
        )

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
