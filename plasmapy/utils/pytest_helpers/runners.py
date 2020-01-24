from plasmapy.utils.pytest_helpers.inputs import (
    InvalidTestError,
    FunctionTestInputs,
    ClassAttributeTestInputs,
    ClassMethodTestInputs,
)

__all__ =


def function_test_runner(
        expected,
        function,
        args=(),
        kwargs=None,
        *,
        rtol=1e-6,
        atol=None,
):
    """
    Test that calling a function with particular arguments results in
    the expected outcome.

    Parameters
    ----------
    function

    args

    kwargs

    expected

    rtol

    atol

    Raises
    ------

    Examples
    --------

    """

    ...


def method_test_runner(
        cls,
        cls_args=None,
        cls_kwargs=None,
        method=None,
        method_args=None,
        method_kwargs=None,
        expected=None,
        *
        rtol=1e-6,
        atol=None,
):
    """
    Test that calling a class method results in the expected outcome.

    Parameters
    ----------
    cls

    cls_args

    cls_kwargs

    method

    method_args

    method_kwargs

    expected

    Raises
    ------

    Examples
    --------

    """

    ...

def attr_test_runner(
        cls,
        cls_args=None,
        cls_kwargs=None,
        attr=None,
        expected=None,
        *,
        rtol=1e-6,
        atol=None,
):
    """
    Test that accessing a class attribute results in the expected outcome.

    Parameters
    ----------
    cls

    cls_args

    cls_kwargs

    attr

    expected

    rtol

    atol

    Examples
    --------

    """

    ...