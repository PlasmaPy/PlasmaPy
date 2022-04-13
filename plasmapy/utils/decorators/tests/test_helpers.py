"""Tests for module :mod:`plasmapy.utils.decorators`."""
import inspect
import pytest

from unittest import mock

from plasmapy.utils.decorators.helpers import modify_docstring, preserve_signature


# --------------------------------------------------------------------------------------
# Test Decorator `modify_docstring`
# --------------------------------------------------------------------------------------
class TestModifyDocstring:
    # create function to test on
    @staticmethod
    def func_simple_docstring(x: float, y: float) -> float:
        """A simple docstring."""
        return x + y

    @staticmethod
    def func_complex_docstring(x: float, y: float) -> float:
        """
        A simple docstring.

        Parameters
        ----------
        x: float
            first argument

        y: float
            second argument

        Returns
        -------
        float
            addition or arguments
        """
        return x + y

    def test_no_arg_exception(self):
        """Raise exception if decorator is used without modifying docstring."""
        with pytest.raises(TypeError):
            modify_docstring(self.func_simple_docstring)

    def test_save_original_doc(self):
        original_doc = self.func_simple_docstring.__doc__
        wfoo = modify_docstring(prepend="Hello")(self.func_simple_docstring)
        assert hasattr(wfoo, "__original_doc__")
        assert wfoo.__original_doc__ == original_doc

    @pytest.mark.parametrize(
        "prepend, append, expected",
        [(5, None, TypeError), (None, 5, TypeError)],
    )
    def test_raises(self, prepend, append, expected):
        with pytest.raises(expected):
            modify_docstring(
                prepend=prepend, append=append, func=self.func_simple_docstring
            )

    def test_preserve_signature(self):
        wfoo = modify_docstring(prepend="Hello")(self.func_simple_docstring)
        assert hasattr(wfoo, "__signature__")
        assert wfoo.__signature__ == inspect.signature(self.func_simple_docstring)

    @pytest.mark.parametrize(
        "prepend, append, func_name, additions",
        [
            (
                "Hello",
                "Goodbye",
                "func_simple_docstring",
                (["Hello", ""], ["", "Goodbye"]),
            ),
            (
                "Hello",
                "Goodbye",
                "func_complex_docstring",
                (["Hello", ""], ["", "Goodbye"]),
            ),
            ("Hello", None, "func_simple_docstring", (["Hello", ""], [])),
            (None, "Goodbye", "func_simple_docstring", ([], ["", "Goodbye"])),
            (
                "\n".join(
                    ["    Hello", "    ", "        * item 1", "            * item 2"]
                ),
                None,
                "func_simple_docstring",
                (["Hello", "", "* item 1", "    * item 2", ""], []),
            ),
            (
                None,
                "\n".join(
                    [
                        "    Notes",
                        "    -----",
                        "    ",
                        "        * item 1",
                        "            * item 2",
                    ]
                ),
                "func_complex_docstring",
                ([], ["", "Notes", "-----", "", "    * item 1", "        * item 2"]),
            ),
        ],
    )
    def test_modification(self, prepend, append, func_name, additions):
        func = getattr(self, func_name)

        expected = "\n".join(
            additions[0] + inspect.cleandoc(func.__doc__).splitlines() + additions[1]
        )

        wfunc = modify_docstring(prepend=prepend, append=append)(func)
        assert wfunc.__doc__ == expected

    def test_arguments_passed(self):
        wfunc = modify_docstring(prepend="Hello")(self.func_simple_docstring)
        assert wfunc(5, 4) == 9


# --------------------------------------------------------------------------------------
# Test Decorator `preserve_signature`
# --------------------------------------------------------------------------------------
def test_preserve_signature():
    # create function to mock
    def foo(x: float, y: float) -> float:
        return x + y

    mock_foo = mock.Mock(side_effect=foo, name="mock_foo", autospec=True)
    mock_foo.__name__ = "mock_foo"
    mock_foo.__signature__ = inspect.signature(foo)

    # -- '__signature__' attribute DNE yet --
    # decorate
    wfoo = preserve_signature(foo)

    # test
    assert hasattr(wfoo, "__signature__")
    assert wfoo.__signature__ == inspect.signature(foo)
    assert wfoo(2, 3) == 5

    # -- '__signature__' attribute is already defined --
    # decorate
    wfoo = preserve_signature(mock_foo)
    assert hasattr(wfoo, "__signature__")
    assert wfoo.__signature__ == inspect.signature(foo)
    assert wfoo(2, 3) == 5
    assert mock_foo.called

    # reset
    mock_foo.reset_mock()
