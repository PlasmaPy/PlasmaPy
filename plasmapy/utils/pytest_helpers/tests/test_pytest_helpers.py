import astropy.units as u
import pytest
import warnings

from typing import Any

from plasmapy.particles import Particle
from plasmapy.tests.helpers import (
    MissingExceptionFail,
    MissingWarningFail,
    TypeMismatchFail,
    UnexpectedExceptionFail,
    UnexpectedResultFail,
)
from plasmapy.utils.code_repr import call_string
from plasmapy.utils.exceptions import PlasmaPyError, PlasmaPyWarning
from plasmapy.utils.pytest_helpers import run_test, run_test_equivalent_calls


def generic_function(*args, **kwargs):
    return None


def raise_exception(*args, **kwargs):
    raise PlasmaPyError(
        f"This exception was raised by raise_exception with:\n\n"
        f"  args = {args}\n"
        f"kwargs = {kwargs}\n"
    )


def issue_warning(*args, **kwargs) -> int:
    warnings.warn(f"\n{args}\n{kwargs}", PlasmaPyWarning)
    return 42


def adams_number(*args, **kwargs):
    return 42


def return_quantity(*args, should_warn: bool = False):
    if should_warn:
        warnings.warn("", UserWarning)
    return 5 * u.m / u.s


def return_arg(arg: Any, should_warn: bool = False) -> Any:
    if should_warn:
        warnings.warn("", UserWarning)
    return arg


f_args_kwargs_expected_whaterror = [
    [adams_number, 0, {"y": 1}, 42, None],
    [adams_number, (1,), {"y": 1}, 42, None],
    [adams_number, (2, 1), {}, 42, None],
    [adams_number, 3, {"y": 1}, 6 * 9, UnexpectedResultFail],
    [adams_number, (4,), {"y": 1}, 6 * 9, UnexpectedResultFail],
    [adams_number, (5, 1), {}, 6 * 9, UnexpectedResultFail],
    [raise_exception, 6, {"y": 1}, PlasmaPyError, None],
    [raise_exception, (7,), {"y": 1}, PlasmaPyError, None],
    [raise_exception, (8, 1), {}, PlasmaPyError, None],
    [raise_exception, 9, {"y": 1}, TypeError, UnexpectedExceptionFail],
    [raise_exception, (10,), {"y": 1}, TypeError, UnexpectedExceptionFail],
    [raise_exception, (11, 1), {}, Exception, UnexpectedExceptionFail],
    [issue_warning, 12, {"y": 1}, PlasmaPyWarning, None],
    [issue_warning, (13,), {"y": 1}, PlasmaPyWarning, None],
    [issue_warning, (14, 1), {}, PlasmaPyWarning, None],
    [issue_warning, 15, {"y": 1}, (42, UserWarning), MissingWarningFail],
    [issue_warning, (16,), {"y": 1}, (42, UserWarning), MissingWarningFail],
    [issue_warning, (17, 1), {}, (42, UserWarning), MissingWarningFail],
    [return_quantity, (18), {}, 5 * u.m / u.s, None],
    [return_quantity, (19), {}, u.m / u.s, None],
    [return_quantity, (20), {}, u.barn * u.Mpc, u.UnitsError],
    [return_quantity, (21), {}, 4 * u.m / u.s, UnexpectedResultFail],
    [return_quantity, (22), {}, 5 * u.kg / u.s, u.UnitsError],
    [return_quantity, (23), {"should_warn": True}, (5 * u.m / u.s, UserWarning), None],
    [
        return_quantity,
        (24),
        {"should_warn": False},
        (5 * u.m / u.s, UserWarning),
        MissingWarningFail,
    ],
    [return_arg, u.kg / u.K, {}, u.kg / u.K, None],
    [return_arg, u.kg / u.K, {}, u.kg / u.N, u.UnitsError],
    [return_arg, u.kg, {}, u.g, u.UnitsError],
    [return_arg, u.C, {"should_warn": True}, (u.C, UserWarning), None],
    [adams_number, 1, {"x": 1}, u.pc, u.UnitsError],
    [return_arg, Particle("p+"), {}, Particle("proton"), None],
    [return_arg, Particle("e+"), {}, Particle("e-"), UnexpectedResultFail],
    [return_arg, Particle("mu+"), {}, type, TypeMismatchFail],
    [return_arg, (2,), {}, IOError, MissingExceptionFail],
]


@pytest.mark.parametrize(
    "f, args, kwargs, expected, whaterror", f_args_kwargs_expected_whaterror
)
def test_run_test(f, args, kwargs, expected, whaterror):
    """
    Test the behavior of the test helper function.

    The arguments `f`, `args`, `kwargs`, and `expected` are to be
    passed directly to `~plasmapy.utils.run_test`.  If the test is
    expected to pass, then `whaterror` should be set to `None`.  If the
    test is expected to fail, then `whaterror` should be set to the
    exception that should be raised during the test failure.

    This function does not test that the exception messages are correct;
    it only checks that the exception messages do not raise any
    unexpected exceptions.

    Parameters
    ----------
    f : callable
        The function or method to be tested.

    args : tuple
        The positional arguments to be sent to f.

    kwargs : dict
        The keyword arguments to be sent to f.

    expected : object
        The expected result from the test, which can be an `object`, an
        exception, a warning, or a `tuple` containing the expected
        result and the type of warning.

    whaterror : exception or None
        The type of exception that the test helper function is supposed
        to raise, or `None` if the test helper function is not supposed
        to raise an exception.

    """

    try:
        if whaterror is None:
            run_test(f, args, kwargs, expected)
        else:
            with pytest.raises(whaterror):
                run_test(f, args, kwargs, expected)
                pytest.fail(
                    f"run_test did not raise an exception for "
                    f"{call_string(f, args, kwargs, color=None)} "
                    f"with expected = {repr(expected)} and "
                    f"whaterror = {repr(whaterror)}."
                )
    except Exception as spectacular_exception:
        raise Exception(
            f"An unexpected exception was raised while running "
            f"{call_string(f, args, kwargs, color=None)} with "
            f"expected = {repr(expected)} and "
            f"whaterror = {repr(whaterror)}."
        ) from spectacular_exception


def test_run_test_rtol():
    run_test(return_arg, 1.0, {}, 0.999999, rtol=1.1e-6)


def test_run_test_rtol_failure():
    with pytest.raises(UnexpectedResultFail):
        run_test(return_arg, 1.0, {}, 0.999999, rtol=1e-7)
        pytest.fail("No exception raised for rtol test.")


def test_run_test_atol():
    run_test(return_arg, 1.0, {}, 0.999999, atol=1.1e-6, rtol=0.0)


def test_run_test_atol_failure():
    with pytest.raises(UnexpectedResultFail):
        run_test(return_arg, (1.0,), {}, 0.999999, atol=1e-7)
        pytest.fail("No exception raised for atol test.")


def func(x, raise_exception=False, issue_warning=False):
    if raise_exception:
        raise ValueError("")
    elif issue_warning:
        warnings.warn("", UserWarning)
    return x


inputs_table = [
    (func, 1, 1),
    (func, (2,), {}, 2),
    (func, 3, {"raise_exception": True}, ValueError),
    (func, 4, {"issue_warning": True}, UserWarning),
    (func, 5, {"issue_warning": True}, (5, UserWarning)),
]


@pytest.mark.parametrize("inputs", inputs_table)
def test_func(inputs):
    """Test cases originally put in the docstring."""
    run_test(inputs)


# TODO: organize this in a namedtuple to improve readability?

run_test_equivalent_calls_table = [
    # cases like inputs = (func, (args, kwargs), (args, kwargs), (args, kwargs))
    [(return_arg, [(1,), {}], [(1,), {}]), None],
    [(return_arg, [(1,), {}], [(86,), {}]), UnexpectedResultFail],
    [(return_arg, [(1,), {}], [(1,), {}], [(1,), {}], [(1,), {}], [(1,), {}]), None],
    # cases like inputs = [(func, args, kwargs), (func, args, kwargs))
    (((return_arg, (1,), {}), (return_arg, (1,), {})), None),
    (
        (
            (return_arg, (1,), {}),
            (return_arg, (1,), {}),
            (return_arg, (1,), {}),
            (return_arg, (1,), {}),
        ),
        None,
    ),
    (
        (
            (return_arg, (1,), {}),
            (return_arg, (2,), {}),
            (return_arg, (3,), {}),
            (return_arg, (4,)),
        ),
        UnexpectedResultFail,
    ),
    # cases where there are no kwargs
    ((return_arg, [1], [1]), None),
    ((return_arg, ["this"], ["that"]), UnexpectedResultFail),
    ([(return_arg, [1], [1])], None),
    # cases where there are no kwargs and the args are not in tuples or lists
    ((return_arg, 1, 1, 1, 1), None),
    ((return_arg, 1, 1, 1, 8794), UnexpectedResultFail),
    (((lambda x, y: x + y, (1, 0), {}), (lambda x, y: x * y, (1, 1), {})), None),
    (((lambda x, y: x + y, (1, 0), {}), (lambda x, y: x * y, (1, 1), {})), None),
    (
        ((lambda x, y: x + y, (1, 0), {}), (lambda x, y: x * y, (1, 59), {})),
        UnexpectedResultFail,
    ),
]


@pytest.mark.parametrize("inputs, error", run_test_equivalent_calls_table)
def test_run_test_equivalent_calls(inputs, error):
    if error is None:
        try:
            run_test_equivalent_calls(*inputs)
        except Exception as exc:
            raise Exception(
                f"Unexpected exception for run_tests_equivalent_calls with "
                f"the following inputs =\n\n  {inputs}"
            ) from exc
    elif issubclass(error, Exception):
        with pytest.raises(error):
            run_test_equivalent_calls(inputs)


def test_run_test_equivalent_calls_types():
    run_test_equivalent_calls(return_arg, 1, 1.0, require_same_type=False)
    with pytest.raises(UnexpectedResultFail):
        run_test_equivalent_calls(return_arg, 1, 1.0)
