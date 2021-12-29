"""Tests of generic decorator functionality."""

from plasmapy.utils.decorators.generic import GenericDecorator

generic_decorator = GenericDecorator().as_decorator


@generic_decorator
def decorated_function(x, *, y, z=5):
    return x, y, z


def test_decorated_function():
    assert decorated_function(2, y=3) == (2, 3, 5)
