import collections
import pytest
import numpy as np
from astropy import units as u

from plasmapy.tests.helper.actual import ActualTestOutcome
from plasmapy.tests.helper.expected import ExpectedTestOutcome
from plasmapy.tests.helper.inputs import FunctionTestInputs

from plasmapy.tests.helper.comparators import (
    CompareActualExpected,
    CompareValues,
    _get_unit,
    _units_are_compatible,
)
from plasmapy.tests.helper.exceptions import InvalidTestError

from plasmapy.tests.helper.tests.sample_functions import (
    return_42,
    return_42_meters,
    return_np_array,
    return_sum_of_two_args_and_kwargs,
    return_none,
    raise_sample_exception,
    issue_sample_warning_and_return_42,
    SampleException,
    SampleExceptionSubclass,
    SampleWarning,
    SampleWarningSubclass,
)

noargs = ()
nokwargs = {}


@pytest.mark.parametrize(
    "function, args, kwargs, expected, test_should_pass",
    [
        (raise_sample_exception, noargs, nokwargs, SampleException, True),
        (raise_sample_exception, noargs, nokwargs, Exception, False),
        (raise_sample_exception, noargs, nokwargs, SampleExceptionSubclass, False),
        (raise_sample_exception, noargs, nokwargs, 42, False),
        (raise_sample_exception, noargs, nokwargs, SampleWarning, False),
        (return_42, noargs, nokwargs, 42, True),
        (return_42, noargs, nokwargs, 6 * 9, False),
        (return_42, noargs, nokwargs, 42.0, False),
        (return_42, noargs, nokwargs, SampleException, False),
        (return_42, noargs, nokwargs, (SampleWarning, 42), False),
        (return_42, noargs, nokwargs, (42, SampleWarning), False),
        (return_42_meters, noargs, nokwargs, 42.0 * u.m, True),
        (return_42_meters, noargs, nokwargs, 42 * (1 + 1e-15) * u.m, True),
        (return_42_meters, noargs, nokwargs, 4200.0 * u.cm, False),
        (return_42_meters, noargs, nokwargs, 42 * u.kg, False),
        (return_42_meters, noargs, nokwargs, 6 * 9 * u.m, False),
        (return_42_meters, noargs, nokwargs, u.m, True),
        (return_42_meters, noargs, nokwargs, u.cm, False),
        (return_none, noargs, nokwargs, None, True),
        (return_none, noargs, nokwargs, False, False),
        (issue_sample_warning_and_return_42, noargs, nokwargs, SampleWarning, True),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            SampleWarningSubclass,
            False,
        ),
        (issue_sample_warning_and_return_42, noargs, nokwargs, Warning, False),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (SampleWarning, 42),
            True,
        ),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (42, SampleWarning),
            True,
        ),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (SampleWarning, 6 * 9),
            False,
        ),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (6 * 9, SampleWarning),
            False,
        ),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (SampleWarningSubclass, 42),
            False,
        ),
        (
            issue_sample_warning_and_return_42,
            noargs,
            nokwargs,
            (42, SampleWarningSubclass),
            False,
        ),
        (issue_sample_warning_and_return_42, noargs, nokwargs, (Warning, 42), False),
        (issue_sample_warning_and_return_42, noargs, nokwargs, (42, Warning), False),
        # if unexpected warnings should cause a failure
        (issue_sample_warning_and_return_42, noargs, nokwargs, 42, False),
        (return_sum_of_two_args_and_kwargs, (2, 3), {"kw1": 5, "kw2": 7}, 17, True),
        (return_sum_of_two_args_and_kwargs, (2, 3), {"kw1": 5, "kw2": 7}, 42, False),
        (
            return_sum_of_two_args_and_kwargs,
            (2, 3),
            {"kw1": 5, "kw2": 7},
            Warning,
            False,
        ),
        (
            return_sum_of_two_args_and_kwargs,
            (2, 3),
            {"kw1": 5, "kw2": 7},
            Exception,
            False,
        ),
        (
            return_sum_of_two_args_and_kwargs,
            (2, 3),
            {"kw1": 5, "kw2": 7},
            (Warning, 17),
            False,
        ),
        (
            return_sum_of_two_args_and_kwargs,
            (2, 3),
            {"kw1": 5, "kw2": 7},
            (17, Warning),
            False,
        ),
        (return_sum_of_two_args_and_kwargs, (2, 3), {"kw1": 5, "kw2": 7}, "17", False),
    ],
)
def test_comparator_actual_expected(function, args, kwargs, expected, test_should_pass):
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

    test_information = (
        "\n\n"
        f"function = {function.__name__}\n"
        f"args = {repr(args)}\n"
        f"kwargs = {repr(kwargs)}\n"
        f"expected = {repr(expected)}\n"
        f"test_should_pass = {test_should_pass}"
    )

    try:
        inputs = FunctionTestInputs(function, args, kwargs)
        actual = ActualTestOutcome(inputs)
        expected = ExpectedTestOutcome(expected)
    except Exception:
        pytest.fail(
            "Problem instantiating classes needed to instantiate "
            "CompareActualExpected" + test_information
        )

    try:
        comparison = CompareActualExpected(actual, expected, atol=None, rtol=1e-7)
    except Exception:
        raise InvalidTestError(
            f"Unable to instantiate CompareActualExpected." f"{test_information}"
        )

    test_information += f"\ncomparison.test_passed = {comparison.test_passed}"

    if comparison.error_message:
        test_information += (
            f"\n\n" f"Presumably incorrect error message: " + comparison.error_message
        )

    if test_should_pass and not comparison.test_passed:
        pytest.fail("Test should have passed, but did not." + test_information)

    if not test_should_pass and comparison.test_passed:
        pytest.fail("Test should not have passed, but did." + test_information)


@pytest.mark.parametrize(
    "this, that, attribute, expected",
    [
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
    ],
)
def test_compare_values_attributes(this, that, attribute, expected):
    """
    Test that the attributes of CompareValues return the correct results.
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


@pytest.mark.parametrize(
    "this, that, when_made_boolean",
    [
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
    ],
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
    """Test that good values of rtol get passed through okay."""
    comparison = CompareValues(1, 1, rtol=rtol)
    if comparison.rtol != rtol:
        pytest.fail(
            f"rtol attribute of CompareValues is not the expected " f"value of {rtol}."
        )


inputs_and_expected_units = [(u.m, u.m), (5 * u.kg * u.s, u.kg * u.s), (1, None)]


@pytest.mark.parametrize("input, expected_unit", inputs_and_expected_units)
def test__get_unit(input, expected_unit):
    gotten_unit = _get_unit(input)
    if gotten_unit != expected_unit:
        pytest.fail(f"_get_unit({input}) is not ")


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
def test__units_are_compatible(unit1, unit2, expected_compatibility):
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
        issue_sample_warning_and_return_42,
        noargs,
        nokwargs,
        43,
        "The command issue_sample_warning_and_return_42() returned a value "
        "of 42, which differs from the expected value of 43. This command "
        "unexpectedly issued the following warnings:\n\n"
        "SampleWarning: warning message",
    ),
    case(
        issue_sample_warning_and_return_42,
        noargs,
        nokwargs,
        (SampleWarning, 43),
        "The command issue_sample_warning_and_return_42() returned a "
        "value of 42, which differs from the expected value of 43.",
    ),
    case(
        issue_sample_warning_and_return_42,
        noargs,
        nokwargs,
        Warning,
        "The command issue_sample_warning_and_return_42() was expected to "
        "issue a Warning, but instead issued the following warning:\n\n"
        "SampleWarning: warning message",
    ),
    case(
        issue_sample_warning_and_return_42,
        noargs,
        nokwargs,
        (SampleWarningSubclass, 43),
        "The command issue_sample_warning_and_return_42() returned a value "
        "of 42, which differs from the expected value of 43. This command "
        "was expected to issue a SampleWarningSubclass, but instead "
        "issued the following warning:\n\n"
        "SampleWarning: warning message",
    ),
    case(
        raise_sample_exception,
        noargs,
        nokwargs,
        Exception,
        "The command raise_sample_exception() raised a SampleException, "
        "instead of an Exception as expected.",
    ),
    case(
        raise_sample_exception,
        noargs,
        nokwargs,
        SampleExceptionSubclass,
        "The command raise_sample_exception() raised a SampleException, "
        "instead of a SampleExceptionSubclass as expected.",
    ),
    case(
        raise_sample_exception,
        noargs,
        nokwargs,
        42,
        "The command raise_sample_exception() unexpectedly raised a "
        "SampleException.",
    ),
    case(
        return_sum_of_two_args_and_kwargs,
        (5.3, 2.42),
        {"kw1": 1.56 * u.dimensionless_unscaled, "kw2": 4.2},
        42,
        "The command return_sum_of_two_args_and_kwargs(5.3, 2.42, kw1=1.56, "
        "kw2=4.2) returned a value of 13.48, which differs from the expected "
        "value of 42. The type of the returned value "
        "(astropy.units.quantity.Quantity) "
        "is different than the type of the expected value (int).",
    ),
]


@pytest.mark.parametrize(
    "func, args, kwargs, expected, errmsg", func_args_kwargs_expected_errmsg
)
def test_compare_actual_expected_errmsg(func, args, kwargs, expected, errmsg):
    """
    Test that CompareActualExpected generates the appropriate error
    messages.

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

    if not errmsg in comparison.error_message:
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
