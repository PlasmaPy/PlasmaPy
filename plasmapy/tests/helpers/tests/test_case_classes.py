"""Test the classes responsible for holding test cases."""

from plasmapy.tests.helpers.cases import (FunctionTestCase, MethodTestCase, AttrTestCase)
from plasmapy.tests.helpers.tests.sample_functions import return_42


def test_instantiation_and_storage():
    """..."""

    function = return_42
    args = (52, 24)
    kwargs = {"x": 38}
    expected = 42
    atol = 1e-7
    rtol = 1e-8

    case = FunctionTestCase(
        function=function,
        args=args,
        kwargs=kwargs,
        expected=expected,
        rtol=rtol,
        atol=atol,
    )

    assert case.function is function
    assert case.args == args
    assert case.kwargs == kwargs
    assert case.expected == expected
    assert case.atol == atol
    assert case.rtol == rtol
