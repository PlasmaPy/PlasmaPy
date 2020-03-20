"""Test the classes responsible for holding test cases."""

import pytest

from plasmapy.tests.helpers.cases import (
    AttrTestCase,
    FunctionTestCase,
    MethodTestCase,
)
from plasmapy.tests.helpers.tests.sample_functions import (
    SampleClass2,
    return_42,
)


class TestFunctionTestCase:
    """Test that `FunctionTestCase` correctly stores its attributes."""

    def setup_class(self):

        self.function = return_42
        self.args = (52, 24)
        self.kwargs = {"x": 38}
        self.expected = 42
        self.atol = 1e-7
        self.rtol = 1e-9

        self.instance = FunctionTestCase(
            expected=self.expected,
            function=self.function,
            args=self.args,
            kwargs=self.kwargs,
            atol=self.atol,
            rtol=self.rtol,
        )

    @pytest.mark.parametrize("attr", ["expected", "args", "kwargs", "atol", "rtol"])
    def test_that_attribute_is_stored_correctly(self, attr):
        assert getattr(self, attr) is getattr(self.instance, attr)

    def test_that_stored_function_has_correct_name(self):
        assert self.function.__name__ == self.instance.function.__name__ == "return_42"


class TestMethodTestCase:
    """Test that `MethodTestCase` correctly stores its attributes."""

    def setup_class(self):

        self.expected = Exception
        self.cls = SampleClass2
        self.method = "method_name"
        self.cls_args = (23, 32)
        self.cls_kwargs = {"x": 835}
        self.method_args = (32, 22)
        self.method_kwargs = {"y": 53}
        self.atol = 1e-23
        self.rtol = 1e-22

        self.instance = MethodTestCase(
            expected=self.expected,
            cls=self.cls,
            method=self.method,
            cls_args=self.cls_args,
            cls_kwargs=self.cls_kwargs,
            method_args=self.method_args,
            method_kwargs=self.method_kwargs,
            atol=self.atol,
            rtol=self.rtol,
        )

    @pytest.mark.parametrize(
        "attr",
        [
            "expected",
            "cls",
            "method",
            "cls_args",
            "cls_kwargs",
            "method_args",
            "method_kwargs",
            "atol",
            "rtol",
        ],
    )
    def test_that_attribute_is_stored_correctly(self, attr):
        assert getattr(self, attr) is getattr(self.instance, attr)


class TestAttrTestCase:
    """Test that `AttrTestCase` correctly stores its attributes."""

    def setup_class(self):

        self.expected = Exception
        self.cls = SampleClass2
        self.attribute = "attr_name"
        self.cls_args = (23, 32)
        self.cls_kwargs = {"x": 835}
        self.atol = 1e-23
        self.rtol = 1e-22

        self.instance = AttrTestCase(
            expected=self.expected,
            cls=self.cls,
            attribute=self.attribute,
            cls_args=self.cls_args,
            cls_kwargs=self.cls_kwargs,
            atol=self.atol,
            rtol=self.rtol,
        )

    @pytest.mark.parametrize(
        "attr",
        ["expected", "cls", "attribute", "cls_args", "cls_kwargs", "atol", "rtol"],
    )
    def test_that_attribute_is_stored_correctly(self, attr):
        assert getattr(self, attr) is getattr(self.instance, attr)
