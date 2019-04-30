from ..pytest_helpers_base import Expected


@pytest.mark.parametrize("exception", [ValueError, Exception, TypeError])
def test_exception(exception):
    expected = Expected(exception)
    assert expected.exception is exception
    assert expected.warning is None
    assert expected.unit is None
    assert expected.value is None


@pytest.mark.parametrize("warning", [UserWarning, Warning])
def test_warning(warning):
    expected = Expected(warning)
    assert expected.warning is Warning
    assert expected.exception is None
    assert expected.unit is None
    assert expected.value is None

import pytest

@pytest.fixture
def set_something_up():
    return 42

def test_something(set_something_up):
    assert set_something_up == 42