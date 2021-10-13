"""
Test module  for `plasmapy.utils.decorators.lite_func.bind_lite_func`.
"""
import pytest

from plasmapy.utils.decorators.lite_func import bind_lite_func
from plasmapy.utils.exceptions import PlasmaPyWarning


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
    print(
        "I am a helper function that support the Lite-Function 'foo_lite'."
    )


@pytest.mark.parametrize(
    "lite_func, attrs, _error",
    [
        # conditions on attrs kwarg
        (foo_lite, "not a dictionary", TypeError),
        (foo_lite, {"lite": lambda x: x}, ValueError),
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
    "lite_func, attrs, _warning",
    [
        (foo_lite, {"bar": "not a functions"}, PlasmaPyWarning),
    ],
)
def test_warns(lite_func, attrs, _warning):
    """Test scenarios that will issue a Warning."""
    with pytest.warns(_warning):
        bind_lite_func(lite_func, attrs=attrs)(foo)
