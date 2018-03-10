import pytest

from ..pytest_helpers import (
    call_string,
    run_test,
    RunTestError,
    UnexpectedExceptionError,
    MissingWarningError,
    MissingExceptionError,
)

from ..exceptions import PlasmaPyWarning, PlasmaPyError

import warnings

def f(*args, **kwargs):
    return None


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


f_args_kwargs_expected_shouldpass = [

    (adams_number, 0, {'y': 1}, 42, True),
    (adams_number, (1,), {'y': 1}, 42, True),
    (adams_number, (2, 1), {}, 42, True),

    (adams_number, 3, {'y': 1}, 6 * 9, False),
    (adams_number, (4,), {'y': 1}, 6 * 9, False),
    (adams_number, (5, 1), {}, 6 * 9, False),

    (raise_exception, 6, {'y': 1}, PlasmaPyError, True),
    (raise_exception, (7,), {'y': 1}, PlasmaPyError, True),
    (raise_exception, (8, 1), {}, PlasmaPyError, True),

    (raise_exception, 9, {'y': 1}, TypeError, False),
    (raise_exception, (10,), {'y': 1}, TypeError, False),
    (raise_exception, (11, 1), {}, TypeError, False),

    (issue_warning, 12, {'y': 1}, PlasmaPyWarning, True),
    (issue_warning, (13,), {'y': 1}, PlasmaPyWarning, True),
    (issue_warning, (14, 1), {}, PlasmaPyWarning, True),

    (issue_warning, 0, {'y': 1}, (42, UserWarning), False),
    (issue_warning, (0,), {'y': 1}, (42, UserWarning), False),
    (issue_warning, (0, 1), {}, (42, UserWarning), False),

]

@pytest.mark.parametrize(
    "f, args, kwargs, expected, should_pass",
    f_args_kwargs_expected_shouldpass,
)
def test_run_test(f, args, kwargs, expected, should_pass):
    if should_pass:
        run_test(f, args, kwargs, expected)
    else:
        with pytest.raises(RunTestError, message = (
                f"run_test did not raise an exception for:\n\n"
                f"  {call_string(f, args, kwargs)}\n\n"
                f"with expected = {expected}")):
            run_test(f, args, kwargs, expected)


