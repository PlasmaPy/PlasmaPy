"""Contains test comparators: tools that automate comparisons."""

from numbers import Number
from typing import Union, Tuple, Any, Optional

from astropy import units as u

from plasmapy.tests.helper.expected import ExpectedTestOutcome
from plasmapy.tests.helper.actual import ActualTestOutcome

from plasmapy.utils.formatting.formatting import (
    _name_with_article,
    _string_together_warnings_for_printing,
    _object_name,
)

from plasmapy.tests.helper.exceptions import (
    Failed,
    UnexpectedResultError,
    InconsistentTypeError,
    UnexpectedExceptionError,
    MissingExceptionError,
    UnexpectedWarningError,
    MissingWarningError,
    InvalidTestError,
    ExceptionMismatchError,
    WarningMismatchError,
)

__all__ = ["CompareActualExpected"]


def _get_unit(obj: Any):
    """
    Return the unit corresponding to a unit or Quantity, or return
    `None` when ``obj`` is not a unit or Quantity.
    """

    if isinstance(obj, u.UnitBase):
        return obj
    elif isinstance(obj, u.Quantity):
        return obj.unit
    else:
        return None


def _units_are_compatible(unit1: Optional[u.UnitBase], unit2: Optional[u.UnitBase]):
    """
    Return `True` if ``unit1`` and ``unit2`` are compatible with each
    other. This function considers `None` (representing an `object`
    other than a `~astropy.units.Quantity` which has no units) as
    equivalent to `~astropy.units.dimensionless_unscaled`.
    """

    for unit in [unit1, unit2]:
        if unit is not None and not isinstance(unit, u.UnitBase):
            raise TypeError(f"{unit} is not a unit or None.")

    if unit1 is unit2:
        return True

    new_unit1 = unit1 if isinstance(unit1, u.UnitBase) else u.dimensionless_unscaled
    new_unit2 = unit2 if isinstance(unit2, u.UnitBase) else u.dimensionless_unscaled

    if new_unit1.physical_type != new_unit2.physical_type:
        return False
    elif new_unit1.physical_type == new_unit2.physical_type != "unknown":
        return True

    try:
        new_unit1.to(new_unit2)
    except u.UnitConversionError:
        return False
    else:
        return True


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
        or `~astropy.units.allclose`.  Defaults to ``1e-8``.

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
        rtol: Union[Number, u.Quantity] = 1e-8,
        atol: Optional[Union[Number, u.Quantity]] = None,
    ):

        self._this = this
        self._that = that

        self.rtol = rtol
        self.atol = atol

    @property
    def values(self) -> Tuple[Any, Any]:
        """
        The two values that are being compared.
        """

        return self._this, self._that

    @property
    def types(self) -> Tuple[type, type]:
        """
        The types of the two values that are being compared.
        """

        return type(self._this), type(self._that)

    @property
    def units(self) -> Tuple[Optional[u.UnitBase], Optional[u.UnitBase]]:
        """
        A `tuple` containing the units of ``this`` and ``that`` or `None`
        for each value that does not have units.
        """

        return _get_unit(self._this), _get_unit(self._that)

    @property
    def units_are_identical(self) -> bool:
        """
        `True` if the units of the two values being compared are identical
        to each other or if neither of the values has a unit, and `False` otherwise.
        """

        return self.units[0] is self.units[1]

    @property
    def units_are_compatible(self) -> bool:
        """
        `True` if the units of the two values are dimensionally compatible, and
        `False` otherwise.
        """
        return _units_are_compatible(*self.units)

    @property
    def rtol(self) -> Union[Number, u.Quantity]:
        """
        The relative tolerance to be used by `~astropy.units.isclose` or
        `~astropy.units.allclose`.  If ``rtol`` is a
        `~astropy.units.Quantity`, then it must be dimensionless.
        """

        return self._rtol

    @rtol.setter
    def rtol(self, new_rtol):

        try:
            if 0 <= new_rtol < 1:
                self._rtol = new_rtol
            else:
                raise ValueError
        except Exception:
            raise InvalidTestError(
                "rtol must be a number or dimensionless Quantity with" "0 <= rtol < 0."
            ) from None

    @property
    def atol(self) -> Optional[Union[Number, u.Quantity]]:
        """
        The absolute tolerance to be used by `~astropy.units.isclose` or
        `~astropy.units.allclose`.  If the objects being compared are two
        `~astropy.units.Quantity` instances, then ``atol`` must have
        consistent dimensions.
        """

        return self._atol

    @atol.setter
    def atol(self, new_atol: Optional[Union[Number, u.Quantity]]):

        self._atol = new_atol

    @property
    def are_identical(self) -> bool:
        """
        Return `True` if the two values refer to the same object,
        and `False` otherwise.
        """

        return self._this is self._that

    @property
    def are_equal(self) -> bool:
        """
        Return `True` if the two values evaluate as equal to each other,
        and `False` otherwise.
        """

        if self.are_quantities:
            return u.allclose(*self.values, atol=None, rtol=0, equal_nan=True)

        try:
            equality = self.values[0] == self.values[1]
        except Exception:
            return False

        if isinstance(equality, bool):
            return equality

        try:
            return all(equality)
        except Exception:
            pass

        if not self.units_are_compatible:
            return False

        if self.are_quantity_and_unit:
            return False

        raise InvalidTestError(f"Cannot determine whether or not {self.values} are equal.")

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

        return isinstance(self.values[0], u.Quantity) and isinstance(self.values[1], u.Quantity)

    @property
    def are_quantity_and_unit(self) -> bool:
        """
        `True` if the first `object` being compared is a `~astropy.units.Quantity`
        instance and the second `object` being compared is an instance of
        (a subclass of) `~astropy.units.UnitBase`, and `False` otherwise.
        """

        return isinstance(self.values[0], u.Quantity) and isinstance(self.values[1], u.UnitBase)

    @property
    def are_allclose(self) -> bool:
        """
        `True` if the compared values are element-wise equal to ``that``
        within an absolute tolerance of ``atol`` and a relative tolerance
        of ``rtol``, and `False` otherwise.  This attribute calls
        `~astropy.units.allclose` to make this determination.

        Notes
        -----
        This attribute will return `True` if ``this`` and ``that`` refer
        to the same `object` or are equal to each other.  Otherwise,
        if `~astropy.units.allclose` raises a `TypeError`, then this
        attribute will return `False`.

        Raises
        ------
        InvalidTestError
            If the units of ``atol`` are incompatible with units shared
            by both ``this`` and ``that``.
        """

        if self.are_identical or self.are_equal:
            return True

        try:
            return u.allclose(*self.values, rtol=self.rtol, atol=self.atol, equal_nan=True)
        except u.UnitsError as exc1:
            if self.units_are_compatible and isinstance(self.atol, u.Quantity):
                if not _units_are_compatible(self.units[0], self.atol.unit):
                    raise InvalidTestError(
                        f"The units of atol ({self.atol}) are incompatible with "
                        f"the units of the Quantity instances being compared "
                        f"{self.units}."
                    ) from exc1
            return False
        except TypeError:
            return False

    def __bool__(self):
        """
        Return `True` if the test should pass, and `False` otherwise.
        """

        if self.are_quantity_and_unit:
            return self.units_are_identical

        if not self.have_same_types:
            return False
        elif not self.units_are_identical:
            return False

        if self.are_identical or self.are_allclose or self.are_equal:
            return True
        else:
            return False


class CompareActualExpected:
    """
    A class to compare the actual and expected results of a test, and
    generate an appropriate error message if necessary.

    Parameters
    ----------
    actual : ActualTestOutcome

    expected : ExpectedTestOutcome

    rtol : dimensionless number, optional, keyword-only
        The relative tolerance to be supplied to `~astropy.units.isclose`
        or `~astropy.units.allclose`.  Defaults to ``1e-8``.

    atol : number or `~astropy.units.Quantity`, optional, keyword-only
        The absolute tolerance to be supplied to `astropy.units.isclose`
        or `~astropy.units.allclose`.  If ``atol`` is a
        `~astropy.units.Quantity`, then it must have the same units as
        ``this`` and ``that``.  Defaults to zero in the appropriate units.
    """

    def __init__(
        self,
        actual: ActualTestOutcome,
        expected: ExpectedTestOutcome,
        *,
        rtol: Union[Number, u.Quantity] = 1e-8,
        atol: Optional[Union[Number, u.Quantity]] = 0.0,
    ):

        if not isinstance(actual, ActualTestOutcome):
            raise TypeError("Expecting an instance of ActualTestOutcome")
        if not isinstance(expected, ExpectedTestOutcome):
            raise TypeError("Expecting an instance of ExpectedTestOutcome")

        self._actual = actual
        self._expected = expected

        self._rtol = rtol
        self._atol = atol

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
        The relative tolerance to be used in `astropy.units.allclose`.
        If ``rtol`` is a `~astropy.units.Quantity` instance, then it
        must be dimensionless.
        """

        return self._rtol

    @property
    def atol(self) -> Union[Number, u.Quantity]:
        """
        The absolute tolerance to be used in `~astropy.units.allclose`.
        """

        return self._atol

    @rtol.setter
    def rtol(self, new_rtol):
        self._rtol = new_rtol

    @property
    def test_passed(self) -> bool:
        """
        Return `True` if the actual outcome matches the expected outcome,
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
    def exception(self) -> Exception:
        """
        Return the exception to be raised if the test failed.
        """

        if self.test_passed:
            raise InvalidTestError(
                f"The test of {self.actual.call_string} passed, so no " f"exception is available."
            )
        else:
            return self._exception

    def _add_exception(self, exception):
        """
        Specify the exception associated with the test failure.
        """

        if not issubclass(exception, BaseException):
            raise TypeError("Expecting an exception.")

        more_than_one_thing_is_going_wrong = hasattr(self, "_exception")

        if more_than_one_thing_is_going_wrong:
            self._exception = Failed
        else:
            self._exception = exception

    @property
    def _subject(self):
        """
        Return an appropriate subject for the first sentence in the
        error message.  The call string should be included in each error
        message if and only if it is the first error message.
        """

        is_first_error = not self._error_messages_list
        subject = f"The command {self.actual.call_string}" if is_first_error else "This command"
        return subject

    def _make_exception_mismatch_errmsg_if_necessary(self):
        """
        Compose an error message for tests where a certain type of
        exception was raised, but a different type of exception was
        expected.  If the expected exception was actually raised, do
        nothing.
        """

        expected_exception = self.expected.expected_exception
        actual_exception = self.actual.exception_type

        if actual_exception is expected_exception:
            return

        errmsg = (
            f"{self._subject} raised "
            f"{_name_with_article(self.actual.exception_type)}, instead of "
            f"{_name_with_article(self.expected.expected_exception)} as expected."
        )

        self._add_errmsg(errmsg)
        self._add_exception(ExceptionMismatchError)

    def _make_missing_exception_errmsg(self):
        """
        Compose an error message for tests where an exception should
        have been raised, but was not.
        """

        errmsg = (
            f"The command {self.actual.call_string} did not raise "
            f"{_name_with_article(self.expected.expected_exception)} as expected. "
            f"Instead, this command returned the unexpected value of "
            f"{repr(self.actual.value)}."
        )

        self._add_errmsg(errmsg)
        self._add_exception(MissingExceptionError)

        # TODO: improve representation of the value (as repr doesn't
        #       always result in something particularly readable.)

    def _make_unexpected_exception_errmsg(self):
        """
        Compose an error message for tests where an exception was
        unexpectedly raised.
        """

        errmsg = (
            f"The command {self.actual.call_string} unexpectedly raised "
            f"{_name_with_article(self.actual.exception_type)}."
        )

        self._add_errmsg(errmsg)
        self._add_exception(UnexpectedExceptionError)

    def _make_incompatible_units_errmsg(self):
        """
        Compose an error message for tests where the units of the two
        values being compared cannot be converted to each other.
        """
        actual_unit = _get_unit(self.actual.value)
        expected_unit = _get_unit(self.expected.expected_value)

        incompatible_units_errmsg = (
            f"The units of the returned value ({actual_unit}) are "
            f"incompatible with the units of the expected value "
            f"({expected_unit})."
        )

        self._add_errmsg(incompatible_units_errmsg)
        self._add_exception(u.UnitsError)

    def _make_nonidentical_units_errmsg(self):
        """
        Compose an error message for tests where the units of the two
        values being compared are not identical to each other.
        """

        actual_unit = _get_unit(self.actual.value)
        expected_unit = _get_unit(self.expected.expected_value)

        unit_mismatch_errmsg = (
            f"The units of the returned value ({actual_unit}) are not "
            f"identical to the units of the expected value "
            f"({expected_unit})."
        )

        self._add_errmsg(unit_mismatch_errmsg)
        self._add_exception(u.UnitsError)

    def _make_different_types_errmsg(self):
        """
        Compose an error message for tests where the two values have
        different types.
        """

        actual_type = type(self.actual.value)
        expected_type = type(self.expected.expected_value)

        actual_type_name = _object_name(actual_type, showmodule=True)
        expected_type_name = _object_name(expected_type, showmodule=True)

        errmsg = (
            f"The type of the returned value ({actual_type_name})"
            f" is different than the type of the expected value "
            f"({expected_type_name})."
        )

        self._add_errmsg(errmsg)
        self._add_exception(InconsistentTypeError)

    def _make_value_mismatch_errmsg_if_necessary(self):
        """
        Compose an error message for tests where the expected value does
        not match the value that was actually returned.  If the expected
        and actual values match, then do nothing.
        """

        comparison = CompareValues(
            self.actual.value, self.expected.expected_value, rtol=self.rtol, atol=self.atol,
        )

        if comparison:
            return

        value_mismatch_errmsg = (
            f"{self._subject} returned a value of {self.actual.value}, "
            f"which differs from the expected value of {self.expected.expected_value}."
        )

        self._add_errmsg(value_mismatch_errmsg)

        if not comparison.have_same_types and not comparison.are_quantity_and_unit:
            self._make_different_types_errmsg()
        elif not comparison.units_are_compatible:
            self._make_incompatible_units_errmsg()
        elif not comparison.units_are_identical:
            self._make_nonidentical_units_errmsg()
        else:
            self._add_exception(UnexpectedResultError)

        # TODO: Should we add a method to check whether the len(...) of
        #       the expected and actual outcomes matches or not?  That
        #       could potentially be useful for debugging.

    def _make_missing_warning_errmsg(self):
        """
        Compose an error message for tests where a warning was expected
        to be issued, but was not.
        """

        missing_warning_errmsg = (
            f"{self._subject} did not raise "
            f"{_name_with_article(self.expected.expected_warning)}"
            f"as expected."
        )

        self._add_errmsg(missing_warning_errmsg)
        self._add_exception(MissingWarningError)

    def _make_unexpected_warnings_errmsg(self):
        """
        Compose an error message for tests where warnings were
        unexpectedly issued.
        """

        warnings_for_printing = _string_together_warnings_for_printing(
            self.actual.warning_types, self.actual.warning_messages
        )

        number_of_warnings = len(self.actual.warning_types)

        unexpected_warnings_errmsg = (
            f"{self._subject} unexpectedly issued the following warnings"
            f"{'s' if number_of_warnings > 1 else ''}:"
            f"\n\n"
            f"{warnings_for_printing}"
        )

        self._add_errmsg(unexpected_warnings_errmsg)
        self._add_exception(UnexpectedWarningError)

    def _make_warning_mismatch_errmsg_if_necessary(self):
        """
        Compose an error message for tests where the expected warning
        was not issued but different warning(s) were.
        """

        expected_warning = self.expected.expected_warning
        actual_warnings = self.actual.warning_types
        warning_messages = self.actual.warning_messages

        number_of_warnings = len(actual_warnings)

        if expected_warning in actual_warnings:
            return

        warnings_for_printing = _string_together_warnings_for_printing(
            actual_warnings, warning_messages,
        )

        warning_mismatch_errmsg = (
            f"{self._subject} was expected to issue {_name_with_article(expected_warning)}, "
            f"but instead issued the following warning"
            f"{'s' if number_of_warnings > 1 else ''}:"
            f"\n\n"
            f"{warnings_for_printing}"
        )

        self._add_errmsg(warning_mismatch_errmsg)
        self._add_exception(WarningMismatchError)

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
