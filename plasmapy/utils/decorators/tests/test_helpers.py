"""Tests for module :mod:`plasmapy.utils.decorators`."""
import inspect
from unittest import mock

from ..helpers import preserve_signature


# ----------------------------------------------------------------------------------------
# Test Decorator `preserve_signature`
# ----------------------------------------------------------------------------------------
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
