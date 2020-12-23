import numpy as np
import pytest

from astropy import units as u

from plasmapy.utils.code_repr import (
    _format_quantity,
    _string_together_warnings_for_printing,
    attribute_call_string,
    call_string,
    method_call_string,
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
    """
    Tests that call_string returns a string that is equivalent to the
    function call.
    """
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


@pytest.mark.parametrize(
    "quantity, expected",
    [
        (5.3 * u.m, "5.3*u.m"),
        (5.4 / u.m, "5.4/u.m"),
        (5.5 * u.m ** -2, "5.5*u.m**-2"),
        (5.0 * u.dimensionless_unscaled, "5.0*u.dimensionless_unscaled"),
        (np.array([3.5, 4.2]) * u.m ** -2.5, "np.array([3.5, 4.2])*u.m**-2.5"),
    ],
)
def test_format_quantity(quantity, expected):
    """
    Test that `~astropy.units.Quantity` objects get converted into a
    string as expected.
    """
    assert _format_quantity(quantity) == expected


@pytest.mark.parametrize("not_a_quantity", [3.5, "1.2"])
def test_format_quantity_typeerror(not_a_quantity):
    """"""
    with pytest.raises(TypeError):
        _format_quantity(not_a_quantity)


def test_stringing_together_warning():
    """Test that warning names and messages get strung together correctly."""
    warnings = [UserWarning, DeprecationWarning]
    warning_messages = ["msg1", "msg2"]
    expected = "UserWarning: msg1\n\nDeprecationWarning: msg2"
    assert (
        _string_together_warnings_for_printing(warnings, warning_messages) == expected
    )
