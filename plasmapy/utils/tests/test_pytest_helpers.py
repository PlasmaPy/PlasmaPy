import pytest

from ..pytest_helpers import call_string, run_test

from ..exceptions import PlasmaPyWarning, PlasmaPyError


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
    """Tests that _function_call_string returns a string that is
    equivalent to the function call."""
    assert expected == call_string(f, args, kwargs)


def raise_exception(x, y=None):
    raise PlasmaPyError("I'm sorry, Dave. I'm afraid I can't do that.")


def warn_warning(x, y=None):
    warnings.warn("Danger, Will Robinson!", PlasmaPyWarning)
    return 42


def adams_number(x, y=None):
    return 42


f_args_kwargs_expected_shouldpass = [

    (adams_number, 0, {'y': 1}, 42, True),
    (adams_number, (0,), {'y': 1}, 42, True),
    (adams_number, (0, 1), {}, 42, True),

    (adams_number, 0, {'y': 1}, 6 * 9, False),
    (adams_number, (0,), {'y': 1}, 6 * 9, False),
    (adams_number, (0, 1), {}, 6 * 9, False),

    (raise_exception, 0, {'y': 1}, PlasmaPyError, True),
    (raise_exception, (0,), {'y': 1}, PlasmaPyError, True),
    (raise_exception, (0, 1), {}, PlasmaPyError, True),

    (raise_exception, 0, {'y': 1}, TypeError, False),
    (raise_exception, (0,), {'y': 1}, TypeError, False),
    (raise_exception, (0, 1), {}, TypeError, False),

    (warn_warning, 0, {'y': 1}, PlasmaPyWarning, True),
    (warn_warning, (0,), {'y': 1}, PlasmaPyWarning, True),
    (warn_warning, (0, 1), {}, PlasmaPyWarning, True),

    (warn_warning, 0, {'y': 1}, (42, UserWarning), False),
    (warn_warning, (0,), {'y': 1}, (42, UserWarning), False),
    (warn_warning, (0, 1), {}, (42, UserWarning), False),

]

@pytest.mark.parametrize(
    "f, args, kwargs, expected, should_pass",
    f_args_kwargs_expected_shouldpass,
)
def test_passing_tests(f, args, kwargs, expected, should_pass):
    if should_pass:
        run_test(f, args, kwargs, expected)
    else:
        with pytest.raises(Exception, message=
            f"run_test did not raise an exception with:\n"
            f"f        = {f.__name__}\n"
            f"args     = {args}\n"
            f"kwargs   = {kwargs}\n"
            f"expected = {expected}"
        ):
            run_test(f, args, kwargs, expected)


