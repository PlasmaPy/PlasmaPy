"""
Test module for `plasmapy.utils.decorators.lite_func.bind_lite_func`.
"""
import pytest

from numba import jit, njit

from plasmapy.utils.decorators.lite_func import bind_lite_func


def foo(x):
    """Test function used for decoration."""
    if not isinstance(x, float):
        raise ValueError
    return x


def foo_lite(x):
    """Test Lite-Function to be bound to foo."""
    return x


def bar():
    """
    Test support function for the Lite-Function framework.  To be bound
    to foo.
    """
    print("I am a helper function that support the Lite-Function 'foo_lite'.")


@pytest.mark.parametrize(
    "lite_func, attrs, _error",
    [
        # conditions on attrs kwarg
        (foo_lite, "not a dictionary", TypeError),
        (foo_lite, {"lite": lambda x: x}, ValueError),
        (foo_lite, {"bar": "not a functions"}, ValueError),
        #
        # conditions on lite_func arg
        (6, None, ValueError),  # not a function
        (print, None, ValueError),  # can not be builtin
    ],
)
def test_raises(lite_func, attrs, _error):
    """Test scenarios that will raise an Exception."""
    with pytest.raises(_error):
        bind_lite_func(lite_func, attrs=attrs)(foo)


@pytest.mark.parametrize(
    "lite_func, attrs",
    [
        (foo_lite, None),
        (jit(foo_lite), None),
        (njit(foo_lite), None),
        (foo_lite, {"bar": bar}),
    ],
)
def test_binding(lite_func, attrs):
    """Test that the expected members are bound to the decorated function."""
    dfoo = bind_lite_func(lite_func, attrs=attrs)(foo)

    if attrs is None:
        attrs = {"lite": lite_func}
    elif isinstance(attrs, dict):
        attrs["lite"] = lite_func
    else:
        pytest.fail(
            "Parametrization was not setup correctly!  The value for 'attrs' must"
            "be None or a dictionary."
        )

    for name, func in attrs.items():
        assert hasattr(dfoo, name)
        assert getattr(dfoo, name) == func

    assert dfoo(5.0) == 5.0
    assert dfoo.lite(5.0) == 5.0


@pytest.mark.parametrize(
    "lite_func, attrs",
    [
        (foo_lite, None),
        (foo_lite, {"bar": bar}),
    ],
)
def test_lite_func_dunder(lite_func, attrs):
    """Test that the ``__bound_lite_func__`` dunder is properly defined."""
    dfoo = bind_lite_func(lite_func, attrs=attrs)(foo)

    if attrs is None:
        attrs = {"lite": lite_func}
    elif isinstance(attrs, dict):
        attrs["lite"] = lite_func
    else:
        pytest.fail(
            "Parametrization was not setup correctly!  The value for 'attrs' must"
            "be None or a dictionary."
        )

    for name, func in attrs.items():
        assert hasattr(dfoo, "__bound_lite_func__")
        assert name in dfoo.__bound_lite_func__

        origin = f"{func.__module__}.{func.__name__}"
        assert dfoo.__bound_lite_func__[name] == origin
