import pytest
from astropy import units as u

from plasmapy.utils.formatting import (
    call_string, method_call_string, attribute_call_string,
)


def generic_function(*args, **kwargs):
    return None


def adams_number():
    return 42


# function, args, kwargs, expected
call_string_table = [
    (generic_function, (), {}, "generic_function()"),
    (generic_function, (1), {}, "generic_function(1)"),
    (generic_function, ("x"), {}, "generic_function('x')"),
    (generic_function, (1, "b", {}), {}, "generic_function(1, 'b', {})"),
    (generic_function, (), {"kw": 1}, "generic_function(kw=1)"),
    (generic_function, (), {"x": "c"}, "generic_function(x='c')"),
    (
        generic_function,
        (1, "b"),
        {"b": 42, "R2": "D2"},
        "generic_function(1, 'b', b=42, R2='D2')",
    ),
]


@pytest.mark.parametrize("function,args,kwargs,expected", call_string_table)
def test_call_string(function, args, kwargs, expected):
    """Tests that call_string returns a string that is
    equivalent to the function call."""
    assert expected == call_string(function, args, kwargs)


class SampleClass:
    def method(self, *args, **kwargs):
        pass

    @property
    def attr(self):
        pass


class_method_call_string_table = [
    ((), {}, (), {}, "SampleClass().method()"),
    (
        (1, 2),
        {"a": "x", "b": "y"},
        (3, 4.2),
        {"q": UserWarning},
        "SampleClass(1, 2, a='x', b='y').method(3, 4.2, q=UserWarning)",
    ),
    (
        (Exception),
        {},
        (SampleClass, adams_number),
        {"r": 5.3 * u.m},
        "SampleClass(Exception).method(SampleClass, adams_number, r=5.3*u.m)",
    ),
    (1, {}, 2, {}, "SampleClass(1).method(2)"),
]


@pytest.mark.parametrize(
    "c_args, c_kwargs, m_args, m_kwargs, expected", class_method_call_string_table
)
def test_class_method_call_string(c_args, c_kwargs, m_args, m_kwargs, expected):
    """Test that `method_call_string` returns the expected results."""
    actual = method_call_string(
        SampleClass, "method", c_args, c_kwargs, m_args, m_kwargs
    )
    assert expected == actual


class_attribute_call_string_table = [
    ((), {}, "SampleClass().attr"),
    ((1, 2), {"a": "x", "b": "y"}, "SampleClass(1, 2, a='x', b='y').attr"),
    (1, {}, "SampleClass(1).attr"),
    ({"dict": "ionary"}, {}, "SampleClass({'dict': 'ionary'}).attr"),
]


@pytest.mark.parametrize(
    "c_args, c_kwargs, expected", class_attribute_call_string_table
)
def test_class_attribute_call_string(c_args, c_kwargs, expected):
    """Test that `attribute_call_string` returns the expected results."""
    actual = attribute_call_string(SampleClass, "attr", c_args, c_kwargs)
    assert expected == actual
