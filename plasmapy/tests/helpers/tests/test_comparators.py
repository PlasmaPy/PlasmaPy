"""Test the classes that compare the actual and expected outcomes."""

import collections

import numpy as np
import pytest
from astropy import units as u
from plasmapy.tests.helpers.actual import ActualTestOutcome
from plasmapy.tests.helpers.comparators import (
    CompareActualExpected,
    CompareValues,
    _get_unit,
    _units_are_compatible,
)
from plasmapy.tests.helpers.exceptions import (
    ExceptionMismatchError,
    Failed,
    InconsistentTypeError,
    InvalidTestError,
    MissingExceptionError,
    MissingWarningError,
    UnexpectedExceptionError,
    UnexpectedResultError,
    UnexpectedWarningError,
    WarningMismatchError,
)
from plasmapy.tests.helpers.expected import ExpectedTestOutcome
from plasmapy.tests.helpers.inputs import FunctionTestInputs
from plasmapy.tests.helpers.tests.sample_functions import (
    SampleException,
    SampleExceptionSubclass,
    SampleWarning,
    SampleWarningSubclass,
    issue_warning_return_42,
    raise_exception,
    return_42,
    return_42_meters,
    return_none,
    return_np_array,
    sum_of_args_and_kwargs,
)
from plasmapy.utils.formatting.formatting import _name_with_article

noargs = ()
nokwargs = {}


table_of_function_args_kwargs_expected_and_exception = [
    (raise_exception, noargs, nokwargs, SampleException, None),
    (raise_exception, noargs, nokwargs, Exception, ExceptionMismatchError),
    (
        raise_exception,
        noargs,
        nokwargs,
        SampleExceptionSubclass,
        ExceptionMismatchError,
    ),
    (raise_exception, noargs, nokwargs, 42, UnexpectedExceptionError),
    (raise_exception, noargs, nokwargs, SampleWarning, Failed),
    (return_42, noargs, nokwargs, 42, None),
    (return_42, noargs, nokwargs, 6 * 9, UnexpectedResultError),
    (return_42, noargs, nokwargs, 42.0, InconsistentTypeError),
    (return_42, noargs, nokwargs, SampleException, MissingExceptionError),
    (return_42, noargs, nokwargs, (SampleWarning, 42), MissingWarningError),
    (return_42, noargs, nokwargs, (42, SampleWarning), MissingWarningError),
    (return_42_meters, noargs, nokwargs, 42.0 * u.m, None),
    (return_42_meters, noargs, nokwargs, 42 * (1 + 1e-15) * u.m, None),
    (return_42_meters, noargs, nokwargs, 4200.0 * u.cm, u.UnitsError),
    (return_42_meters, noargs, nokwargs, 42 * u.kg, u.UnitsError),
    (return_42_meters, noargs, nokwargs, 6 * 9 * u.m, UnexpectedResultError),
    (return_42_meters, noargs, nokwargs, u.m, None),
    (return_42_meters, noargs, nokwargs, u.cm, u.UnitsError),
    (return_none, noargs, nokwargs, None, None),
    (return_none, noargs, nokwargs, False, InconsistentTypeError),
    (issue_warning_return_42, noargs, nokwargs, SampleWarning, None),
    (
        issue_warning_return_42,
        noargs,
        nokwargs,
        SampleWarningSubclass,
        WarningMismatchError,
    ),
    (issue_warning_return_42, noargs, nokwargs, Warning, WarningMismatchError),
    (issue_warning_return_42, noargs, nokwargs, (SampleWarning, 42), None),
    (issue_warning_return_42, noargs, nokwargs, (42, SampleWarning), None),
    (
        issue_warning_return_42,
        noargs,
        nokwargs,
        (SampleWarning, 6 * 9),
        UnexpectedResultError,
    ),
    (
        issue_warning_return_42,
        noargs,
        nokwargs,
        (6 * 9, SampleWarning),
        UnexpectedResultError,
    ),
    (
        issue_warning_return_42,
        noargs,
        nokwargs,
        (SampleWarningSubclass, 42),
        WarningMismatchError,
    ),
    (
        issue_warning_return_42,
        noargs,
        nokwargs,
        (42, SampleWarningSubclass),
        WarningMismatchError,
    ),
    (issue_warning_return_42, noargs, nokwargs, (Warning, 42), WarningMismatchError,),
    (issue_warning_return_42, noargs, nokwargs, (42, Warning), WarningMismatchError,),
    (issue_warning_return_42, noargs, nokwargs, 42, UnexpectedWarningError),
    (sum_of_args_and_kwargs, (2, 3), {"kw1": 5, "kw2": 7}, 17, None),
    (sum_of_args_and_kwargs, (2, 3), {"kw1": 5, "kw2": 7}, 42, UnexpectedResultError,),
    (
        sum_of_args_and_kwargs,
        (2, 3),
        {"kw1": 5, "kw2": 7},
        Warning,
        MissingWarningError,
    ),
    (
        sum_of_args_and_kwargs,
        (2, 3),
        {"kw1": 5, "kw2": 7},
        Exception,
        MissingExceptionError,
    ),
    (
        sum_of_args_and_kwargs,
        (2, 3),
        {"kw1": 5, "kw2": 7},
        (Warning, 17),
        MissingWarningError,
    ),
    (
        sum_of_args_and_kwargs,
        (2, 3),
        {"kw1": 5, "kw2": 7},
        (17, Warning),
        MissingWarningError,
    ),
    (
        sum_of_args_and_kwargs,
        (2, 3),
        {"kw1": 5, "kw2": 7},
        "17",
        InconsistentTypeError,
    ),
]


@pytest.mark.parametrize(
    "function, args, kwargs, expected, expected_exception",
    table_of_function_args_kwargs_expected_and_exception,
)
def test_compare_actual_expected(function, args, kwargs, expected, expected_exception):
    """
    Test that `CompareActualExpected` correctly determines whether or
    not a test passes.

    Parameters
    ----------
    function
        The sample function to be tested

    args
        A `tuple` or `list` of positional arguments

    kwargs
        A dictionary of keyword arguments

    expected
        The expected outcome of the test to be provided to
        `CompareActualExpected`.  The value of this ``expected`` will
        sometimes be correct and sometimes  be incorrect in order to
        perform test this class.

    test_should_pass
        `True` if ``function(*args, **kwargs)`` should result in
        ``expected``, and `False` otherwise.  Equivalently, this
        argument should be `True` if ``expected`` is correct and
        `False` if ``expected`` is incorrect.
    """

    test_should_pass = expected_exception is None

    test_information = (
        "\n\n"
        f"function = {function.__name__}\n"
        f"args = {repr(args)}\n"
        f"kwargs = {repr(kwargs)}\n"
        f"expected = {repr(expected)}\n"
        f"test_should_pass = {test_should_pass}\n"
        f"expected_exception = {expected_exception}"
    )

    try:
        inputs = FunctionTestInputs(function, args, kwargs)
        actual = ActualTestOutcome(inputs)
        expected = ExpectedTestOutcome(expected)
    except Exception as exc:
        raise Failed(
            f"The classes needed to instantiate CompareActualExpected "
            f"were not able to be instantiated." + test_information
        ) from exc

    try:
        comparison = CompareActualExpected(actual, expected, atol=None, rtol=1e-7)
    except Exception as exc:
        raise Failed(
            f"Problem instantiating CompareActualExpected." + test_information
        ) from exc

    if test_should_pass is True and comparison.test_passed is False:
        raise Failed(
            "This test was expected to pass, but instead failed." + test_information
        )
    elif test_should_pass is False and comparison.test_passed is True:
        raise Failed(
            "This test was expected to fail, but instead passed." + test_information
        )
    elif test_should_pass is False:
        if comparison.exception is not expected_exception:
            actual_exception_from_test = _name_with_article(comparison.exception)
            expected_exception_for_test = _name_with_article(expected_exception)
            raise ExceptionMismatchError(
                f"This should should raise {expected_exception_for_test}, "
                f"but instead raised {actual_exception_from_test}." + test_information
            )


compared_values_attr_expected = [
    (1, 1, "values", (1, 1)),
    ("1", "1", "types", (str, str)),
    (1, 1, "units", (None, None)),
    (1, 1, "units_are_identical", True),
    (1, 1, "units_are_compatible", True),
    (None, None, "are_identical", True),
    (1, 1, "are_equal", True),
    (0, 0.0, "are_equal", True),
    (1, 1, "are_quantities", False),
    (1, 1, "have_same_types", True),
    (5 * u.m, u.m, "units", (u.m, u.m)),
    (5 * u.m, u.cm, "are_quantity_and_unit", True),
    (5 * u.m, 6 * u.m, "units_are_identical", True),
    (5 * u.m, 6 * u.cm, "units_are_identical", False),
    (5 * u.m, u.cm, "units_are_compatible", True),
    (5 * u.m, 5.000001 * u.m, "are_allclose", True),
    (5 * u.m, 5.01 * u.m, "are_allclose", False),
    ([5.0, 6.0] * u.m, [5.0, 6.0] * u.m, "are_allclose", True),
    ([5.0, 6.0] * u.m, [5.0, 6.00000001] * u.m, "are_allclose", True),
    ([5.0, 6.0] * u.m, [5.0, 6.01] * u.m, "are_allclose", False),
    (5 * u.m, 5.000001 * u.m, "are_equal", False),
    (5 * u.m, 5.01 * u.m, "are_equal", False),
    ([5.0, 6.0] * u.m, [5.0, 6.0] * u.m, "are_equal", True),
    ([5.0, 6.0] * u.m, [5.0, 6.00000001] * u.m, "are_equal", False),
    ([5.0, 6.0] * u.m, [5.0, 6.01] * u.m, "are_equal", False),
]


@pytest.mark.parametrize(
    "this, that, attribute, expected", compared_values_attr_expected
)
def test_compare_values_attributes(this, that, attribute, expected):
    """
    Test that the attributes of `CompareValues` return the correct
    results for different pairs of compared values.
    """

    try:
        comparison = CompareValues(this, that, rtol=1e-4)
    except Exception:
        pytest.fail(f"Unable to instantiate CompareValues for {this} and {that}.")

    try:
        value_of_attribute = getattr(comparison, attribute)
    except Exception:
        pytest.fail(
            f"Unable to access attribute {attribute} for the CompareValues "
            f"instance for {this} and {that}."
        )

    if value_of_attribute != expected and value_of_attribute is not expected:
        pytest.fail(
            f"For the CompareValues instance for {this} and {that}, "
            f"the {repr(attribute)} attribute returns a value of "
            f"{value_of_attribute}, which differs from the expected "
            f"value of {expected}."
        )


compared_values_and_boolean_value = [
    (1, 1, True),
    (1, 2, False),
    ("1", "1", True),
    (5 * u.m, 6 * u.m, False),
    (None, None, True),
    (0, 0.0, False),
    (5 * u.m, u.m, True),
    (5 * u.m, u.cm, False),
    (5 * u.m, 6 * u.m, False),
    (5 * u.m, 6 * u.cm, False),
    (5 * u.m, 500 * u.cm, False),
    (5 * u.m, 5.000001 * u.m, True),
    (5 * u.m, 5.01 * u.m, False),
    ([5.0, 6.0] * u.m, [5.0, 6.0] * u.m, True),
    ([5.0, 6.0] * u.m, [5.0, 6.00000001] * u.m, True),
    ([5.0, 6.0] * u.m, [5.0, 6.01] * u.m, False),
    (np.nan, np.nan, True),
    (np.inf, np.inf, True),
    (np.inf, np.nan, False),
    (1, np.int32(1), False),
]


@pytest.mark.parametrize(
    "this, that, when_made_boolean", compared_values_and_boolean_value
)
def test_compare_values_bool(this, that, when_made_boolean):
    """
    Test that CompareValues.__bool__ returns the expected results.
    """

    try:
        comparison = CompareValues(this, that, rtol=1e-4)
    except Exception:
        pytest.fail(f"Unable to instantiate CompareValues for {this} and {that}.")

    try:
        made_boolean = bool(comparison)
    except Exception:
        pytest.fail(
            f"The CompareValues instance for {this} and {that} cannot "
            f"be made boolean."
        )

    if made_boolean is not when_made_boolean:
        pytest.fail(
            f"The CompareValues instance for {this} and {that} corresponds "
            f"to a boolean value of {made_boolean}, when it should "
            f"actually correspond to {when_made_boolean}."
        )


@pytest.mark.parametrize("rtol", [-1e-14, 1 + 1e-14, 5 * u.m, 1, "..."])
def test_compare_values_rtol_exceptions(rtol):
    """
    Test that bad values of ``rtol`` raise appropriate exceptions
    in CompareValues.
    """

    with pytest.raises(InvalidTestError):
        CompareValues(1, 1, rtol=rtol)
        pytest.fail(
            f"CompareValues with rtol = {rtol} is not raising an "
            f"InvalidTestError as expected."
        )


@pytest.mark.parametrize("rtol", [1e-14, 0.999999, 0.5 * u.dimensionless_unscaled])
def test_compare_values_rtol(rtol):
    """Test that good values of ``rtol`` get passed through okay."""

    comparison = CompareValues(1, 1, rtol=rtol)
    if comparison.rtol != rtol:
        pytest.fail(
            f"rtol attribute of CompareValues is not the expected " f"value of {rtol}."
        )


inputs_and_expected_units = [(u.m, u.m), (5 * u.kg * u.s, u.kg * u.s), (1, None)]


@pytest.mark.parametrize("input, expected_unit", inputs_and_expected_units)
def test_get_unit(input, expected_unit):
    """Test that `_get_unit` returns the expected unit."""

    gotten_unit = _get_unit(input)

    if gotten_unit != expected_unit:
        pytest.fail(f"_get_unit({input}) is not {expected_unit}.")


units_and_expected_compatibility = [
    (u.m, u.m, True),
    (u.m, u.cm, True),
    (u.m, u.kg, False),
    (u.g * u.m * u.s ** 7, u.kg * u.km * u.s ** 7, True),
    (None, u.dimensionless_unscaled, True),
    (None, None, True),
]


@pytest.mark.parametrize(
    "unit1, unit2, expected_compatibility", units_and_expected_compatibility
)
def test_units_are_compatible(unit1, unit2, expected_compatibility):
    """
    Test that `_units_are_compatible` correctly compares the
    compatibility of two units.
    """

    actual_compatibility = _units_are_compatible(unit1, unit2)
    if actual_compatibility is not expected_compatibility:
        pytest.fail(
            f"_units_are_compatible({unit1}, {unit2}) resulted "
            f"in {actual_compatibility}, not the expected value of "
            f"{expected_compatibility}."
        )


case = collections.namedtuple("case", ["func", "args", "kwargs", "expected", "errmsg"])

func_args_kwargs_expected_errmsg = [
    case(
        return_42,
        noargs,
        nokwargs,
        41,
        "The command return_42() returned a value of 42, which differs "
        "from the expected value of 41",
    ),
    case(
        return_42,
        noargs,
        nokwargs,
        SyntaxError,
        "The command return_42() did not raise a SyntaxError as expected. "
        "Instead, this command returned the unexpected value of 42.",
    ),
    case(
        return_42,
        noargs,
        nokwargs,
        np.int32(42),
        "The command return_42() returned a value of 42, which differs "
        "from the expected value of 42. The type of the returned value "
        "(int) is different than the type of the expected value (numpy.int32).",
    ),
    case(
        return_42_meters,
        noargs,
        nokwargs,
        42 * u.cm,
        "The command return_42_meters() returned a value of 42.0 m, which differs "
        "from the expected value of 42.0 cm. The units of the returned value (m) "
        "are not identical to the units of the expected value (cm).",
    ),
    case(
        return_np_array,
        (1.0, 2.0, 3.0),
        nokwargs,
        np.array([1.0, 2.0, 3.1]),
        "The command return_np_array(1.0, 2.0, 3.0) returned a value of "
        "[1. 2. 3.], which differs from the expected value of [1.  2.  3.1].",
    ),
    case(
        issue_warning_return_42,
        noargs,
        nokwargs,
        43,
        "The command issue_warning_return_42() returned a value "
        "of 42, which differs from the expected value of 43. This command "
        "unexpectedly issued the following warnings:\n\n"
        "SampleWarning: warning message",
    ),
    case(
        issue_warning_return_42,
        noargs,
        nokwargs,
        (SampleWarning, 43),
        "The command issue_warning_return_42() returned a "
        "value of 42, which differs from the expected value of 43.",
    ),
    case(
        issue_warning_return_42,
        noargs,
        nokwargs,
        Warning,
        "The command issue_warning_return_42() was expected to "
        "issue a Warning, but instead issued the following warning:\n\n"
        "SampleWarning: warning message",
    ),
    case(
        issue_warning_return_42,
        noargs,
        nokwargs,
        (SampleWarningSubclass, 43),
        "The command issue_warning_return_42() returned a value "
        "of 42, which differs from the expected value of 43. This command "
        "was expected to issue a SampleWarningSubclass, but instead "
        "issued the following warning:\n\n"
        "SampleWarning: warning message",
    ),
    case(
        raise_exception,
        noargs,
        nokwargs,
        Exception,
        "The command raise_exception() raised a SampleException, "
        "instead of an Exception as expected.",
    ),
    case(
        raise_exception,
        noargs,
        nokwargs,
        SampleExceptionSubclass,
        "The command raise_exception() raised a SampleException, "
        "instead of a SampleExceptionSubclass as expected.",
    ),
    case(
        raise_exception,
        noargs,
        nokwargs,
        42,
        "The command raise_exception() unexpectedly raised a SampleException.",
    ),
    case(
        sum_of_args_and_kwargs,
        (5.3, 2.42),
        {"kw1": 1.56 * u.dimensionless_unscaled, "kw2": 4.2},
        42,
        "The command sum_of_args_and_kwargs(5.3, 2.42, kw1=1.56, "
        "kw2=4.2) returned a value of 13.48, which differs from the "
        "expected value of 42. The type of the returned value "
        "(astropy.units.quantity.Quantity) "
        "is different than the type of the expected value (int).",
    ),
]


@pytest.mark.parametrize(
    "func, args, kwargs, expected, errmsg", func_args_kwargs_expected_errmsg
)
def test_compare_actual_expected_errmsg(func, args, kwargs, expected, errmsg):
    """
    Test that `CompareActualExpected` generates the appropriate error messages.

    Parameters
    ----------
    func
        The sample function for the test.

    args
        The positional arguments to be passed to ``func``.

    kwargs
        The keyword arguments to be passed to ``func``.

    expected
        The incorrect expected outcome of a test.

    errmsg : str
        The error message that is expected to be generated (or a portion
        thereof).
    """

    try:
        inputs = FunctionTestInputs(func, args, kwargs)
        actual = ActualTestOutcome(inputs)
        expected = ExpectedTestOutcome(expected)
        comparison = CompareActualExpected(actual, expected, rtol=1e-6)
    except Exception as exc:
        raise InvalidTestError("Unable to instantiate preconditions.") from exc

    if not comparison.error_message:
        pytest.fail("No error message was created.")

    if errmsg not in comparison.error_message:
        pytest.fail(
            f"The instance of CompareActualExpected created from:\n\n"
            f"   func = {func.__name__}\n"
            f"   args = {args}\n"
            f" kwargs = {kwargs}\n\n"
            f"resulted in an error message of:"
            f"\n\n  {comparison.error_message}\n\n"
            f"which differs from the expected error message which "
            f"should contain:\n\n"
            f"  {errmsg}"
        )
