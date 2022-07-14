"""Tests for code representation functions."""

import numpy as np
import pytest

from astropy import units as u
from collections import namedtuple

from plasmapy.utils.code_repr import (
    _code_repr_of_arg,
    _code_repr_of_args_and_kwargs,
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


function_case = namedtuple("function_case", ("func", "args", "kwargs", "expected"))


@pytest.mark.parametrize(
    "func, args, kwargs, expected",
    [
        function_case(
            func=generic_function, args=(), kwargs={}, expected="generic_function()"
        ),
        function_case(
            func=generic_function, args=(1), kwargs={}, expected="generic_function(1)"
        ),
        function_case(
            func=generic_function,
            args=("x"),
            kwargs={},
            expected="generic_function('x')",
        ),
        function_case(
            func=generic_function,
            args=(1, "b", {}),
            kwargs={},
            expected="generic_function(1, 'b', {})",
        ),
        function_case(
            func=generic_function,
            args=(),
            kwargs={"kw": 1},
            expected="generic_function(kw=1)",
        ),
        function_case(
            func=generic_function,
            args=(),
            kwargs={"x": "c"},
            expected="generic_function(x='c')",
        ),
        function_case(
            func=generic_function,
            args=(1, 2, 3),
            kwargs={"a": 4, "b": 5, "c": 6},
            expected="generic_function(1, 2, 3, a=4, b=5, c=6)",
        ),
        function_case(
            func=generic_function,
            args=np.array((1.0, 2.0)),
            kwargs={},
            expected="generic_function(np.array([1., 2.]))",
        ),
        function_case(
            func=generic_function,
            args=np.array((3.0, 4.0)) * u.m,
            kwargs={},
            expected="generic_function(np.array([3., 4.])*u.m)",
        ),
        function_case(
            func=generic_function,
            args=(ValueError,),
            kwargs={"x": TypeError},
            expected="generic_function(ValueError, x=TypeError)",
        ),
        function_case(
            func=generic_function,
            args=(),
            kwargs={"bowl": super},
            expected="generic_function(bowl=super)",
        ),
        function_case(
            func=generic_function,
            args=(1, "b"),
            kwargs={"b": 42, "R2": "D2"},
            expected="generic_function(1, 'b', b=42, R2='D2')",
        ),
        function_case(
            func=generic_function,
            args=np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
            kwargs={},
            expected="generic_function(np.array([1, 2, 3, 4, 5, 6, 7, 8, ...]))",
        ),
        function_case(
            func=generic_function,
            args=[(1, "a")],
            kwargs={},
            expected="generic_function((1, 'a'))",
        ),
        function_case(
            func=generic_function,
            args=([1, "a"],),
            kwargs={},
            expected="generic_function([1, 'a'])",
        ),
    ],
)
def test_call_string(func, args, kwargs, expected):
    """
    Tests that call_string returns a string that is equivalent to the
    function call.
    """
    actual = call_string(func, args, kwargs, max_items=8)
    assert actual == expected, (
        "When call_string is called with:\n"
        f"  function: {func.__name__}\n"
        f"  args:     {args}\n"
        f"  kwargs:   {kwargs}\n"
        "the actual result does not match the expected result:\n"
        f"  expected: {repr(expected)}\n"
        f"  actual:   {repr(actual)}"
    )


class SampleClass:
    def method(self, *args, **kwargs):
        pass

    @property
    def attr(self):
        pass


method_case = namedtuple(
    "method_case",
    (
        "args_to_cls",
        "kwargs_to_cls",
        "args_to_method",
        "kwargs_to_method",
        "expected",
    ),
)


@pytest.mark.parametrize(
    "args_to_cls, kwargs_to_cls, args_to_method, kwargs_to_method, expected",
    [
        method_case(
            args_to_cls=(),
            kwargs_to_cls={},
            args_to_method=(),
            kwargs_to_method={},
            expected="SampleClass().method()",
        ),
        method_case(
            args_to_cls=(1, 2),
            kwargs_to_cls={"a": "x", "b": "y"},
            args_to_method=(3, 4.2),
            kwargs_to_method={"q": UserWarning},
            expected="SampleClass(1, 2, a='x', b='y').method(3, 4.2, q=UserWarning)",
        ),
        method_case(
            args_to_cls=(Exception),
            kwargs_to_cls={},
            args_to_method=(SampleClass, adams_number),
            kwargs_to_method={"r": 5.3 * u.m},
            expected="SampleClass(Exception).method(SampleClass, adams_number, r=5.3*u.m)",
        ),
        method_case(
            args_to_cls=1,
            kwargs_to_cls={},
            args_to_method=2,
            kwargs_to_method={},
            expected="SampleClass(1).method(2)",
        ),
    ],
)
def test_method_call_string(
    args_to_cls, kwargs_to_cls, args_to_method, kwargs_to_method, expected
):
    """Test that `method_call_string` returns the expected results."""
    actual = method_call_string(
        cls=SampleClass,
        method="method",
        args_to_cls=args_to_cls,
        kwargs_to_cls=kwargs_to_cls,
        args_to_method=args_to_method,
        kwargs_to_method=kwargs_to_method,
    )
    assert actual == expected, (
        "When method_call_string is called with:\n"
        f"  cls:      SampleClass\n"
        f"  method:   'method'\n"
        f"  c_args:   {args_to_cls}\n"
        f"  c_kwargs: {kwargs_to_cls}\n"
        f"  m_args:   {args_to_method}\n"
        f"  m_kwargs: {kwargs_to_method}\n"
        "the actual outcome does not match the expected result:\n"
        f"  expected: {repr(expected)}\n"
        f"  actual:   {repr(actual)}"
    )


attribute_case = namedtuple(
    "attribute_case",
    (
        "args_to_cls",
        "kwargs_to_cls",
        "expected",
    ),
)


@pytest.mark.parametrize(
    "args_to_cls, kwargs_to_cls, expected",
    [
        attribute_case(args_to_cls=(), kwargs_to_cls={}, expected="SampleClass().attr"),
        attribute_case(
            args_to_cls=(1, 2),
            kwargs_to_cls={"a": "x", "b": "y"},
            expected="SampleClass(1, 2, a='x', b='y').attr",
        ),
        attribute_case(args_to_cls=1, kwargs_to_cls={}, expected="SampleClass(1).attr"),
        attribute_case(
            args_to_cls={"key": "value"},
            kwargs_to_cls={},
            expected="SampleClass({'key': 'value'}).attr",
        ),
    ],
)
def test_attribute_call_string(args_to_cls, kwargs_to_cls, expected):
    """Test that `attribute_call_string` returns the expected results."""
    actual = attribute_call_string(SampleClass, "attr", args_to_cls, kwargs_to_cls)
    assert actual == expected, (
        "When attribute_call_string is called with:\n"
        f"  cls:      SampleClass\n"
        f"  attr:     'attr'\n"
        f"  c_args:   {args_to_cls}\n"
        f"  c_kwargs: {kwargs_to_cls}\n"
        "the actual outcome does not match the expected result:\n"
        f"  expected: {repr(expected)}\n"
        f"  actual:   {repr(actual)}"
    )


ndarray_case = namedtuple("ndarray_case", ("array_inputs", "max_items", "expected"))


@pytest.mark.parametrize(
    "array_inputs, max_items, expected",
    [
        ndarray_case(array_inputs=[0], max_items=np.inf, expected="np.array([0])"),
        ndarray_case(
            array_inputs=[0.0, 1.0], max_items=np.inf, expected="np.array([0., 1.])"
        ),
        ndarray_case(
            array_inputs=[0.0, 1.142],
            max_items=np.inf,
            expected="np.array([0., 1.142])",
        ),
        ndarray_case(
            array_inputs=[[0, 1, 2], [3, 4, 5]],
            max_items=np.inf,
            expected="np.array([[0, 1, 2], [3, 4, 5]])",
        ),
        ndarray_case(
            array_inputs=[np.inf], max_items=np.inf, expected="np.array([np.inf])"
        ),
        ndarray_case(
            array_inputs=[np.nan], max_items=np.inf, expected="np.array([np.nan])"
        ),
        ndarray_case(
            array_inputs=[np.nan, np.inf, -np.inf],
            max_items=np.inf,
            expected="np.array([np.nan, np.inf, -np.inf])",
        ),
        ndarray_case(array_inputs=[1], max_items=1, expected="np.array([1])"),
        ndarray_case(array_inputs=[1, 2], max_items=1, expected="np.array([1, ...])"),
        ndarray_case(
            array_inputs=[1, 2, 3, 4], max_items=2, expected="np.array([1, 2, ...])"
        ),
        ndarray_case(
            array_inputs=[[1, 2, 3], [4, 5, 6]],
            max_items=5,
            expected="np.array([[1, 2, 3], [4, 5, ...]])",
        ),
        ndarray_case(
            array_inputs=[[[1, 2], [3, 4]], [[5, 6], [7, 8]]],
            max_items=1,
            expected="np.array([[[1, ...]]])",
        ),
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
    assert actual == expected, (
        f"The representation of an ndarray for array_inputs = {array_inputs} "
        f"and max_items = {max_items} is not the expected result:\n"
        f"  expected: {repr(expected)}\n"
        f"  actual:   {repr(actual)}\n"
    )

    if max_items >= array.size:
        recreated_array = eval(actual)
        assert np.allclose(array, recreated_array, atol=1e-15, equal_nan=True), (
            f"Evaluating the representation of an ndarray for array_inputs = "
            f"{array_inputs} and max_items = {max_items} does not recreate "
            f"the original array, as expected."
            f"array:           {array}\n"
            f"recreated_array: {recreated_array}"
        )


quantity_case = namedtuple("QuantityTestCases", ("quantity", "expected"))


@pytest.mark.parametrize(
    "quantity, expected",
    [
        quantity_case(quantity=5.3 * u.m, expected="5.3*u.m"),
        quantity_case(quantity=5.4 / u.m, expected="5.4/u.m"),
        quantity_case(quantity=5.5 * u.m**-2, expected="5.5*u.m**-2"),
        quantity_case(
            quantity=u.Quantity(5.0), expected="5.0*u.dimensionless_unscaled"
        ),
        quantity_case(
            quantity=np.array([3.5, 4.2]) * u.m**-2.5,
            expected="np.array([3.5, 4.2])*u.m**-2.5",
        ),
    ],
)
def test__code_repr_of_quantity(quantity, expected):
    """
    Test that `astropy.units.Quantity` objects get converted into a
    string as expected.
    """
    actual = _code_repr_of_quantity(quantity)
    assert actual == expected, (
        f"_code_repr_of_quantity for {quantity} is not producing the "
        f"expected result:\n"
        f"expected: {repr(expected)}\n"
        f"actual:   {repr(actual)}"
    )


def test__string_together_warnings_for_printing():
    """Test that warning names and messages get strung together correctly."""
    warnings = [UserWarning, DeprecationWarning]
    warning_messages = ["msg1", "msg2"]
    expected = "UserWarning: msg1\n\nDeprecationWarning: msg2"
    actual = _string_together_warnings_for_printing(warnings, warning_messages)
    assert actual == expected, (
        f"_string_together_warnings_for_printing is not producing the "
        f"expected result:\n"
        f"  expected: {repr(expected)}\n"
        f"  actual: {repr(actual)}"
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
    assert code_repr_of_arg == expected, (
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
        (np.array([1.0, 2.0]), {"k": 3}, "np.array([1., 2.]), k=3"),
        (np.array([1.0, 2.0]) * u.m / u.s, {}, "np.array([1., 2.])*u.m/u.s"),
    ],
)
def test__code_repr_of_args_and_kwargs(args, kwargs, expected):
    """
    Test that `_code_repr_of_args_and_kwargs` returns a string containing
    the positional and keyword arguments, as they would appear in a
    function call.
    """
    args_and_kwargs = _code_repr_of_args_and_kwargs(args, kwargs)
    assert args_and_kwargs == expected, (
        f"_code_repr_of_args_and_kwargs with the following arguments:\n"
        f"  args:   {args}\n"
        f"  kwargs: {kwargs}\n"
        f"is not returning the expected string:\n"
        f"  expected: {repr(expected)}\n"
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
def test__name_with_article(obj, expected):
    """
    Test that `_name_with_article` returns the expected string, which
    contains ``"a "`` or ``"an "`` followed by the name of ``obj``.
    """
    name_with_article = _name_with_article(obj)
    assert name_with_article == expected, (
        f"When calling _name_with_article for {obj}, expecting "
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
    """Test that `_object_name` produces the expected output."""
    actual_name = _object_name(obj, showmodule=showmodule)
    assert actual_name == expected_name, (
        f"For obj = {repr(obj)}, the expected output of _obj_name with "
        f"showmodule = {showmodule} does not match the actual output:\n"
        f"  expected: {expected_name}\n"
        f"  actual:   {actual_name}"
    )
