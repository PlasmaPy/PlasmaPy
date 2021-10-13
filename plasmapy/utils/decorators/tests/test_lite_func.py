import pytest

from plasmapy.utils.decorators.lite_func import bind_lite_func
from plasmapy.utils.exceptions import PlasmaPyWarning


class TestBindLiteFunc:
    """
    Test class for `plasmapy.utils.decorators.lite_func.bind_lite_func`.
    """

    @staticmethod
    def foo(x):
        if not isinstance(x, float):
            raise ValueError
        return x

    @staticmethod
    def foo_lite(x):
        return x

    @staticmethod
    def bar():
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
    def test_raises(self, lite_func, attrs, _error):
        """Test scenarios that will raise an Exception."""
        with pytest.raises(_error):
            bind_lite_func(lite_func, attrs=attrs)(self.foo)

    @pytest.mark.skip
    def test_warns(self):
        ...
