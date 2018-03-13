import pytest
import warnings
import astropy.units as u
from typing import Any

from ..pytest_helpers import (
    call_string,
    run_test,
    RunTestError,
    UnexpectedExceptionError,
    MissingWarningError,
    MissingExceptionError,
)
from ..exceptions import PlasmaPyWarning, PlasmaPyError


def f(*args, **kwargs):
    return None


def raise_exception(*args, **kwargs):
    raise PlasmaPyError(
        f"This exception was raised by raise_exception with:\n\n"
        f"  args = {args}\n"
        f"kwargs = {kwargs}\n")


def issue_warning(*args, **kwargs) -> int:
    warnings.warn(f"\n{args}\n{kwargs}", PlasmaPyWarning)
    return 42


def adams_number(*args, **kwargs):
    return 42


def return_quantity(should_warn: bool = False):
    if should_warn:
        warnings.warn("It was the best of times, it was the worst of times", UserWarning)
    return 5 * u.m / u.s


def return_arg(arg: Any, should_warn: bool = False) -> Any:
    if should_warn:
        warnings.warn("It was the age of wisdom, it was the age of foolishness", UserWarning)
    return arg


# f, args, kwargs, expected
call_string_table = [
    (f, (), {}, "f()"),
    (f, (1), {}, "f(1)"),
    (f, ('x'), {}, "f('x')"),
    (f, (1, 'b', {}), {}, "f(1, 'b', {})"),
    (f, (), {'kw': 1}, "f(kw=1)"),
    (f, (), {'x': 'c'}, "f(x='c')"),
    (f, (1, 'b'), {'b': 42, 'R2': 'D2'}, "f(1, 'b', b=42, R2='D2')"),
    (run_test, run_test, {run_test: run_test},
     'run_test(run_test, run_test=run_test)'),
]


@pytest.mark.parametrize("f,args,kwargs,expected", call_string_table)
def test_call_string(f, args, kwargs, expected):
    """Tests that call_string returns a string that is
    equivalent to the function call."""
    assert expected == call_string(f, args, kwargs)


f_args_kwargs_expected_whaterror = [

    (adams_number, 0, {'y': 1}, 42, None),
    (adams_number, (1,), {'y': 1}, 42, None),
    (adams_number, (2, 1), {}, 42, None),

    (adams_number, 3, {'y': 1}, 6 * 9, RunTestError),
    (adams_number, (4,), {'y': 1}, 6 * 9, RunTestError),
    (adams_number, (5, 1), {}, 6 * 9, RunTestError),

    (raise_exception, 6, {'y': 1}, PlasmaPyError, None),
    (raise_exception, (7,), {'y': 1}, PlasmaPyError, None),
    (raise_exception, (8, 1), {}, PlasmaPyError, None),

    (raise_exception, 9, {'y': 1}, TypeError, UnexpectedExceptionError),
    (raise_exception, (10,), {'y': 1}, TypeError, UnexpectedExceptionError),
    (raise_exception, (11, 1), {}, TypeError, UnexpectedExceptionError),

    (issue_warning, 12, {'y': 1}, PlasmaPyWarning, None),
    (issue_warning, (13,), {'y': 1}, PlasmaPyWarning, None),
    (issue_warning, (14, 1), {}, PlasmaPyWarning, None),

    (issue_warning, 0, {'y': 1}, (42, UserWarning), MissingWarningError),
    (issue_warning, (0,), {'y': 1}, (42, UserWarning), MissingWarningError),
    (issue_warning, (0, 1), {}, (42, UserWarning), MissingWarningError),

    (return_quantity, (), {}, 5 * u.m / u.s, None),
    (return_quantity, (), {}, u.m / u.s, None),
    (return_quantity, (), {}, u.barn * u.Mpc, u.UnitsError),
    (return_quantity, (), {}, 4 * u.m / u.s, RunTestError),
    (return_quantity, (), {}, 5 * u.kg / u.s, u.UnitsError),
    (return_quantity, (), {'should_warn': True}, (5 * u.m / u.s, UserWarning), None),
    (return_quantity, (), {'should_warn': False}, (5 * u.m / u.s, UserWarning), MissingWarningError),

    (return_arg, u.kg / u.K, {}, u.kg / u.K, None),
    (return_arg, u.kg / u.K, {}, u.kg / u.N, u.UnitsError),
    (return_arg, u.kg, {}, u.g, u.UnitsError),
    (return_arg, u.C, {'should_warn': True}, (u.C, UserWarning), None),

]



@pytest.mark.parametrize(
    "f, args, kwargs, expected, whaterror",
    f_args_kwargs_expected_whaterror,
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
            with pytest.raises(whaterror, message = (
                    f"run_test did not raise an exception for:\n\n"
                    f"  {call_string(f, args, kwargs)}\n\n"
                    f"with expected = {str(expected)}")):
                run_test(f, args, kwargs, expected)
    except Exception as spectacular_exception:
        raise Exception(
            f"An unexpected exception was raised for:\n"
            f"f = {f.__name__}\n"
            f"args = {args}\n"
            f"kwargs = {kwargs}\n"
            f"expected = {expected}\n"
            f"whaterror = {whaterror}") from spectacular_exception


def test_run_test_tolerances():

    with pytest.raises(RunTestError, message="No exception raised for rtol test."):
        run_test(return_arg, 1.0, {}, 0.999999, rtol=1e-7)

    run_test(return_arg, 1.0, {}, 0.999999, rtol=1.1e-6)

    with pytest.raises(RunTestError, message="No exception raised for atol test."):
        run_test(return_arg, (1.0,), {}, 0.999999, atol=1e-7)

    run_test(return_arg, 1.0, {}, 0.999999, atol=1.1e-6)
