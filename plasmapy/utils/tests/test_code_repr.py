"""Tests for code representation functions."""

import numpy as np
import pytest

from astropy import units as u

from plasmapy.utils.code_repr import (
    _code_repr_of_arg,
    _code_repr_of_args_and_kwargs,
    _code_repr_of_ndarray,
    _code_repr_of_quantity,
    _exception_name_with_indef_article,
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
        ([0], np.inf, "np.array([0])"),
        ([0.0, 1.0], np.inf, "np.array([0., 1.])"),
        ([0.0, 1.142], np.inf, "np.array([0., 1.142])"),
        ([[0, 1, 2], [3, 4, 5]], np.inf, "np.array([[0, 1, 2], [3, 4, 5]])"),
        ([np.inf], np.inf, "np.array([np.inf])"),
        ([np.nan], np.inf, "np.array([np.nan])"),
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


@pytest.mark.parametrize(
    "arg, expected",
    [
        (1, "1"),
        ("asdf", "'asdf'"),
        (3.42, "3.42"),
        ([3.42, 3.84], "[3.42, 3.84]"),
        (4.2 * u.m, "4.2*u.m"),
        (np.array([1, 2, 3]), "np.array([1, 2, 3])"),
        (np.array([[1, 2], [3, 4]]), "np.array([[1, 2], [3, 4]])"),
        (np.array([1.0, 2.0, 3.0]) * u.m / u.s, "np.array([1., 2., 3.])*u.m/u.s"),
    ],
)
def test__code_repr_of_arg(arg, expected):
    """
    Test that _code_repr_of_arg correctly transforms arguments into a
    `str` the represents how the arg would appear in code.
    """
    code_repr_of_arg = _code_repr_of_arg(arg)
    if code_repr_of_arg != expected:
        pytest.fail(
            f"_code_repr_of_arg is not returning the expected result.\n"
            f"arg = {arg}\n"
            f"expected:         {expected}\n"
            f"code_repr_of_arg: {code_repr_of_arg}"
        )


@pytest.mark.parametrize(
    "args, kwargs, expected",
    [
        ((), {}, ""),
        (1, {"a": "b"}, "1, a='b'"),
        ((1, 2, 3), {"a": "A", "b": "B"}, "1, 2, 3, a='A', b='B'"),
    ],
)
def test__code_repr_of_args_and_kwargs(args, kwargs, expected):
    """
    Test that `_code_repr_of_args_and_kwargs` returns a string containing
    the positional and keyword arguments, as they would appear in a
    function call.
    """
    args_and_kwargs = _code_repr_of_args_and_kwargs(args, kwargs)
    if args_and_kwargs != expected:
        pytest.fail(
            f"_code_repr_of_args_and_kwargs with the following arguments:\n"
            f"  args:   {args}\n"
            f"  kwargs: {kwargs}\n"
            f"is not returning the expected string:\n"
            f"  expected: {repr(expected)}\n"
            f"  actual:   {repr(args_and_kwargs)}"
        )


@pytest.mark.parametrize(
    "obj, expected",
    [
        (ArithmeticError, "an ArithmeticError"),
        (EOFError, "an EOFError"),
        (IndexError, "an IndexError"),
        (OSError, "an OSError"),
        (TypeError, "a TypeError"),
        (UserWarning, "a UserWarning"),
        (UnicodeError, "a UnicodeError"),
        (ValueError, "a ValueError"),
    ],
)
def test__exception_name_with_indef_article(obj, expected):
    """
    Test that `_exception_name_with_indef_article` returns the expected string, which
    contains ``"a "`` or ``"an "`` followed by the name of ``obj``.
    """
    name_with_article = _exception_name_with_indef_article(obj)
    if name_with_article != expected:
        pytest.fail(
            f"For calling _exc_name_with_indef_article for {obj}, expecting "
            f"{repr(expected)} but got {repr(name_with_article)}."
        )


@pytest.mark.parametrize(
    "obj, showmodule, expected_name",
    [
        (ValueError, False, "ValueError"),
        (TypeError, True, "TypeError"),
        (np.ndarray, False, "ndarray"),
        (np.ndarray, True, "np.ndarray"),
        (super, False, "super"),
        (super, True, "super"),  # hide module for builtins even if showmodule==True
        (u.Unit, False, "Unit"),
        (u.Unit, True, "u.Unit"),
    ],
)
def test__object_name(obj, showmodule, expected_name):
    actual_name = _object_name(obj, showmodule=showmodule)
    if actual_name != expected_name:
        pytest.fail(
            f"For obj = {repr(obj)}, the expected output of _obj_name with "
            f"showmodule = {showmodule} does not match the actual output:\n"
            f"  expected: {expected_name}\n"
            f"  actual:   {actual_name}"
        )
