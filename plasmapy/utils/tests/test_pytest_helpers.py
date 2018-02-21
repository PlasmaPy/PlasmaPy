import pytest

from ..pytest_helpers import (
    _function_call_string,
    run_test_of_function,
)


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
]


@pytest.mark.parametrize("f,args,kwargs,expected", call_string_table)
def test__function_call_string(f, args, kwargs, expected):
    """Tests that _function_call_string returns a string that is
    equivalent to the function call."""
    assert expected == _function_call_string(f, args, kwargs)
