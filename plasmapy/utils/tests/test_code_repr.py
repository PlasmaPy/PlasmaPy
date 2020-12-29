import numpy as np
import pytest

from astropy import units as u

from plasmapy.utils.code_repr import (
    _code_repr_of_arg,
    _code_repr_of_args_and_kwargs,
    _code_repr_of_keyword_name,
    _code_repr_of_ndarray,
    _code_repr_of_quantity,
    _name_with_article,
    _object_name,
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
    ({"key": "value"}, {}, "SampleClass({'key': 'value'}).attr"),
]


@pytest.mark.parametrize(
    "c_args, c_kwargs, expected", class_attribute_call_string_table
)
def test_class_attribute_call_string(c_args, c_kwargs, expected):
    """Test that `attribute_call_string` returns the expected results."""
    actual = attribute_call_string(SampleClass, "attr", c_args, c_kwargs)
    assert expected == actual


@pytest.mark.parametrize(
    "array_inputs, max_items, expected",
    [
        ([0], np.inf, 'np.array([0])'),
        ([0., 1.], np.inf, 'np.array([0., 1.])'),
        ([0.0, 1.142], np.inf, 'np.array([0., 1.142])'),
        ([[0, 1, 2], [3, 4, 5]], np.inf, "np.array([[0, 1, 2], [3, 4, 5]])"),
        ([np.inf], np.inf, 'np.array([np.inf])'),
        ([np.nan], np.inf, 'np.array([np.nan])'),
        ([np.nan, np.inf, -np.inf], np.inf, "np.array([np.nan, np.inf, -np.inf])"),
        ([1], 1, "np.array([1])"),
        ([1, 2], 1, "np.array([1, ...])"),
        ([1, 2, 3, 4], 2, "np.array([1, 2, ...])"),
        ([[1, 2, 3], [4, 5, 6]], 5, "np.array([[1, 2, 3], [4, 5, ...]])"),
        ([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], 1, "np.array([[[1, ...]]])"),
    ],
)
def test__code_repr_of_ndarray(array_inputs, max_items, expected):
    """
    Test that `numpy.ndarray` objects get converted into a string as
    expected.  Subsequently, test that evaluating the code representation
    of an ndarray returns an array equal to the array created from
    `array_inputs`.
    """
    array = np.array(array_inputs)
    actual = _code_repr_of_ndarray(array, max_items=max_items)
    if actual != expected:
        pytest.fail(
            f"The representation of an ndarray for array_inputs = {array_inputs} "
            f"and max_items = {max_items} is not the expected result:\n"
            f"expected: {repr(expected)}\n"
            f"actual:   {repr(actual)}\n"
        )

    if max_items >= array.size:
        recreated_array = eval(actual)
        if not np.allclose(array, recreated_array, atol=1e-15, equal_nan=True):
            pytest.fail(
                f"Evaluating the representation of an ndarray for array_inputs = "
                f"{array_inputs} and max_items = {max_items} does not recreate "
                f"the original array, as expected."
                f"array:           {array}\n"
                f"recreated_array: {recreated_array}"
            )

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
def test_code_repr_of_quantity(quantity, expected):
    """
    Test that `astropy.units.Quantity` objects get converted into a
    string as expected.
    """
    assert _code_repr_of_quantity(quantity) == expected


@pytest.mark.parametrize("not_a_quantity", [3.5, "1.2"])
def test_code_repr_of_quantity_typeerror(not_a_quantity):
    """
    Test that `_code_repr_of_quantity` raises a `TypeError` when
    supplied with an object other than a `astropy.units.Quantity`.
    """
    with pytest.raises(TypeError):
        _code_repr_of_quantity(not_a_quantity)


def test_string_together_warnings():
    """Test that warning names and messages get strung together correctly."""
    warnings = [UserWarning, DeprecationWarning]
    warning_messages = ["msg1", "msg2"]
    expected = "UserWarning: msg1\n\nDeprecationWarning: msg2"
    assert (
        _string_together_warnings_for_printing(warnings, warning_messages) == expected
    )
