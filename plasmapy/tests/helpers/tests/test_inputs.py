"""Test the classes that store and check the inputs for a test."""

import pytest

from plasmapy.tests.helpers.exceptions import InvalidTestError
from plasmapy.tests.helpers.inputs import (
    _validate_args,
    _validate_kwargs,
    ClassAttributeTestInputs,
    ClassMethodTestInputs,
    FunctionTestInputs,
)

sample_args = (2, 3)
sample_args_list = [5, 7]
sample_kwargs = {"kwarg1": 13, "kwarg2": 17}
bad_kwargs_nonstring_key = {"kwarg1": 19, 23: 29}
bad_kwargs_wrong_type = [31, 37]
sample_args_to_method = (41, 43)
sample_kwargs_to_method = {"method_kwarg1": 47, "method_kwarg2": 53}
sample_value = 59


def add_args_and_kwargs(arg1, arg2, kwarg1=None, kwarg2=None):
    """Sum up the positional and keyword arguments."""
    return sum([arg1, arg2, kwarg1, kwarg2])


def double_argument(arg):
    """Double the argument."""
    return arg * 2


def accept_no_args_and_return_value():
    """Return a sample value."""
    return sample_value


class SampleClassNoArgs:
    """A sample class that requires no positional or keyword arguments to instantiate."""

    @property
    def sample_attribute(self):
        """Return a sample value."""
        return sample_value

    def sample_method(
        self, method_arg1, method_arg2, method_kwarg1=None, method_kwarg2=None
    ):
        """Return the sum of the positional and keyword arguments to the method."""
        return method_arg1 + method_arg2 + method_kwarg1 + method_kwarg2


class SampleClass(SampleClassNoArgs):
    """A sample class that takes positional and keyword arguments."""

    def __init__(self, arg1, arg2, kwarg1=None, kwarg2=None):
        pass


@pytest.mark.parametrize(
    "function, args, kwargs, attribute, expected",
    [
        (
            add_args_and_kwargs,
            sample_args,
            sample_kwargs,
            "function",
            add_args_and_kwargs,
        ),
        (add_args_and_kwargs, sample_args, sample_kwargs, "args", sample_args),
        (
            add_args_and_kwargs,
            sample_args_list,
            sample_kwargs,
            "args",
            sample_args_list,
        ),
        (add_args_and_kwargs, sample_args, sample_kwargs, "kwargs", sample_kwargs),
        (double_argument, sample_value, None, "args", (sample_value,)),
        (accept_no_args_and_return_value, None, None, "args", ()),
        (accept_no_args_and_return_value, None, None, "kwargs", {}),
        (
            accept_no_args_and_return_value,
            None,
            None,
            "call_string",
            "accept_no_args_and_return_value()",
        ),
        (
            add_args_and_kwargs,
            sample_args,
            sample_kwargs,
            "call_string",
            "add_args_and_kwargs(2, 3, kwarg1=13, kwarg2=17)",
        ),
    ],
)
def test_function_test_inputs(function, args, kwargs, attribute, expected):
    """Test the attributes of FunctionTestInputs."""
    instance = FunctionTestInputs(function, args, kwargs)
    value_of_attribute = getattr(instance, attribute)
    if value_of_attribute != expected:
        pytest.fail(
            f"FunctionTestInputs({args}, {kwargs}).{attribute} returns "
            f"{value_of_attribute}, instead of the expected value of "
            f"{expected}."
        )


@pytest.mark.parametrize(
    "function, args, kwargs, expected",
    [
        (
            add_args_and_kwargs,
            sample_args,
            sample_kwargs,
            add_args_and_kwargs(*sample_args, **sample_kwargs),
        ),
        (double_argument, sample_value, None, double_argument(sample_value),),
        (accept_no_args_and_return_value, None, None, sample_value),
    ],
)
def test_function_test_inputs_call(function, args, kwargs, expected):
    """
    Test that the ``call`` method on FunctionTestInputs returns the
    value expected from calling the function with the supplied args
    and kwargs.
    """
    instance = FunctionTestInputs(function, args, kwargs)
    result_of_call = instance.call()
    if result_of_call != expected:
        pytest.fail(
            f"FunctionTestInputs({function.__name__}, {args}, {kwargs}).call() "
            f"returns {result_of_call}, instead of the expected value of "
            f"{expected}."
        )


@pytest.mark.parametrize(
    "function, args, kwargs",
    [
        (add_args_and_kwargs, sample_args, bad_kwargs_nonstring_key),
        (add_args_and_kwargs, sample_args, bad_kwargs_wrong_type),
        ("not a function", None, None),
    ],
)
def test_function_test_inputs_errors(function, args, kwargs):
    """
    Test that FunctionTestInputs raises an InvalidTestError when passed
    an invalid function or keyword arguments.
    """
    with pytest.raises(InvalidTestError):
        FunctionTestInputs(function, args, kwargs)
        name = function.__name__ if hasattr(function, "__name__") else repr(function)
        pytest.fail(
            f"FunctionTestInputs({name}, {args}, {kwargs}) did not "
            f"raise an InvalidTestError, as expected."
        )


@pytest.mark.parametrize(
    "cls, args, kwargs, attribute_being_tested, expected",
    [
        (SampleClass, sample_args, sample_kwargs, "cls", SampleClass),
        (SampleClass, sample_args, sample_kwargs, "args_to_cls", sample_args),
        (SampleClass, sample_args_list, sample_kwargs, "args_to_cls", sample_args_list),
        (SampleClass, None, None, "args_to_cls", ()),
        (SampleClass, 42, None, "args_to_cls", (42,)),
        (SampleClass, sample_args, sample_kwargs, "kwargs_to_cls", sample_kwargs),
        (SampleClass, sample_args, None, "kwargs_to_cls", {}),
        (SampleClass, sample_args, None, "attribute", "sample_attribute"),
        (SampleClass, None, None, "call_string", "SampleClass().sample_attribute"),
        (
            SampleClass,
            sample_args,
            None,
            "call_string",
            "SampleClass(2, 3).sample_attribute",
        ),
        (
            SampleClass,
            1,
            {"a": 2},
            "call_string",
            "SampleClass(1, a=2).sample_attribute",
        ),
    ],
)
def test_class_attr_test_inputs(cls, args, kwargs, attribute_being_tested, expected):
    """
    Test that the attributes of a ``ClassAttributeTestInputs`` instance
    return the expected values for a sample class with an attribute
    named ``'sample_attribute'``.

    Parameters
    ----------
    cls
        The sample class being passed to ``ClassAttributeTestInputs``

    args
        Positional arguments being passed to ``ClassAttributeTestInputs``

    kwargs : `dict` or `None`
        Keyword arguments being passed to ``ClassAttributeTestInputs``

    attribute_being_tested : str
        An attribute of ``ClassAttributeTestInputs``

    expected
        The expected value corresponding to ``attribute_being_tested``
        for a ``ClassAttributeTestInputs`` instance.

    """
    attribute_of_cls = "sample_attribute"
    instance = ClassAttributeTestInputs(cls, attribute_of_cls, args, kwargs)
    value_of_attribute = getattr(instance, attribute_being_tested)
    if value_of_attribute != expected:
        pytest.fail(
            f"SampleClass({cls.__name__}, {attribute_of_cls}, {args}, {kwargs})"
            f".{attribute_being_tested} returns {value_of_attribute} "
            f"instead of the expected value of {expected}."
        )


@pytest.mark.parametrize(
    "cls, args, kwargs, expected",
    [
        (SampleClass, sample_args, sample_kwargs, sample_value),
        (SampleClassNoArgs, None, None, sample_value),
    ],
)
def test_class_attr_test_inputs_call(cls, args, kwargs, expected):
    """
    Test that ``ClassAttributeTestInputs.call`` instance
    return the expected values for a sample class with an attribute
    named ``'sample_attribute'``.

    Parameters
    ----------
    cls
        The sample class being passed to ``ClassAttributeTestInputs``

    args
        Positional arguments being passed to ``ClassAttributeTestInputs``

    kwargs : `dict` or `None`
        Keyword arguments being passed to ``ClassAttributeTestInputs``

    expected
        The expected value corresponding to calling the ``call`` method
        for a ``ClassAttributeTestInputs`` instance.

    """
    attribute_of_cls = "sample_attribute"
    instance = ClassAttributeTestInputs(cls, attribute_of_cls, args, kwargs)
    result_of_call = instance.call()
    if result_of_call != expected:
        pytest.fail(
            f"{cls.__name__}(*{args}, **{kwargs}).call() is not "
            f"returning the expected value."
        )


@pytest.mark.parametrize(
    "cls, args, kwargs, attribute_of_cls",
    [
        (SampleClass, sample_args, bad_kwargs_wrong_type, "sample_attribute"),
        (SampleClass, sample_args, bad_kwargs_nonstring_key, "sample_attribute"),
        (SampleClass, sample_args, sample_kwargs, sample_value),
    ],
)
def test_class_attr_test_inputs_errors(cls, args, kwargs, attribute_of_cls):
    """Test that ``ClassAttributeTestInputs`` raises appropriate exceptions."""
    with pytest.raises(InvalidTestError):
        ClassAttributeTestInputs(cls, attribute_of_cls, args, kwargs)
        name = cls.__name__ if hasattr(cls, "__name__") else repr(cls)
        pytest.fail(
            f"ClassAttributeTestInputs is not raising an exception for "
            f"cls = {name}, args = {args}, kwargs = {kwargs}, and "
            f"attribute = {repr(attribute_of_cls)}."
        )


common_inputs = (
    SampleClass,
    sample_args,
    sample_kwargs,
    sample_args_to_method,
    sample_kwargs_to_method,
)


@pytest.mark.parametrize(
    "cls, args_to_cls, kwargs_to_cls, args_to_method, kwargs_to_method, attribute_being_tested, expected",
    [
        (*common_inputs, "cls", SampleClass),
        (*common_inputs, "args_to_cls", sample_args),
        (*common_inputs, "kwargs_to_cls", sample_kwargs),
        (*common_inputs, "args_to_method", sample_args_to_method),
        (*common_inputs, "kwargs_to_method", sample_kwargs_to_method),
        (*common_inputs, "method", "sample_method"),
        (SampleClassNoArgs, None, None, None, None, "cls", SampleClassNoArgs),
        (SampleClassNoArgs, None, None, None, None, "args_to_cls", ()),
        (SampleClassNoArgs, None, None, None, None, "kwargs_to_cls", {}),
        (SampleClassNoArgs, None, None, None, None, "method", "sample_method"),
        (SampleClassNoArgs, None, None, None, None, "args_to_method", ()),
        (SampleClassNoArgs, None, None, None, None, "kwargs_to_method", {}),
        (SampleClassNoArgs, 42, None, None, None, "args_to_cls", (42,)),
        (SampleClassNoArgs, None, None, 42, None, "args_to_method", (42,)),
        (
            SampleClassNoArgs,
            None,
            None,
            None,
            None,
            "call_string",
            "SampleClassNoArgs().sample_method()",
        ),
        (
            *common_inputs,
            "call_string",
            (
                "SampleClass(2, 3, kwarg1=13, kwarg2=17)."
                "sample_method(41, 43, method_kwarg1=47, method_kwarg2=53)"
            ),
        ),
    ],
)
def test_class_method_test_inputs(
    cls,
    args_to_cls,
    kwargs_to_cls,
    args_to_method,
    kwargs_to_method,
    attribute_being_tested,
    expected,
):
    """
    Test that the attributes of a ``ClassMethodTestInputs`` instance
    return the expected values for a sample class with a method
    named ``'sample_method'``.

    Parameters
    ----------
    cls
        The sample class being passed to ``ClassMethodTestInputs``

    args_to_cls
        Positional arguments to be used when instantiating ``cls``

    kwargs_to_cls : `dict` or `None`
        Keyword arguments to be used when instantiating ``cls``

    args_to_method
        Positional arguments to be used when calling the method
        ``sample_method`` in an instance of ``cls``

    kwargs_to_method : `dict` or `None`
        Keyword arguments to be used when calling the method
        ``sample_method`` in an instance of ``cls``

    attribute_being_tested
        Name of the attribute of ``ClassMethodTestInputs`` that is
        being tested

    expected
        The expected value of the attribute of ``ClassMethodTestInputs``
        that is being tested

    """
    method_name = "sample_method"
    instance = ClassMethodTestInputs(
        cls, method_name, args_to_cls, kwargs_to_cls, args_to_method, kwargs_to_method
    )
    value_of_attribute = getattr(instance, attribute_being_tested)
    if value_of_attribute != expected:
        pytest.fail(
            f"Expecting {expected} but got {value_of_attribute} "
            f"for attribute {attribute_being_tested}."
        )


@pytest.mark.parametrize(
    "cls, args_to_cls, kwargs_to_cls, args_to_method, kwargs_to_method, expected",
    [
        (
            *common_inputs,
            sum(sample_args_to_method) + sum(sample_kwargs_to_method.values()),
        )
    ],
)
def test_class_method_test_inputs_call(
    cls, args_to_cls, kwargs_to_cls, args_to_method, kwargs_to_method, expected
):
    """
    Test that the ``call`` method on ClassMethodTestInputs returns
    the expected values for a sample class with a method named
    ``sample_method``.
    """
    method_name = "sample_method"
    instance = ClassMethodTestInputs(
        cls, method_name, args_to_cls, kwargs_to_cls, args_to_method, kwargs_to_method
    )
    value_of_call = instance.call()
    if value_of_call != expected:
        pytest.fail(
            f"The call method on a ClassMethodTestInputs instance is "
            f"returning {value_of_call}, rather than the expected value "
            f"of {expected}.  This case has cls = {cls.__name__}, "
            f"args_to_cls = {args_to_cls}, kwargs_to_cls = {kwargs_to_cls}, "
            f"args_to_method = {args_to_method}, and kwargs_to_method = "
            f"{kwargs_to_method}."
        )


@pytest.mark.parametrize(
    "cls, args_to_cls, kwargs_to_cls, args_to_method, kwargs_to_method, method_name",
    [
        (*common_inputs, 42),
        (*common_inputs, "sample_attribute"),
        (SampleClass, None, (), None, None, "sample_method"),
        (SampleClass, None, None, None, (), "sample_method"),
        (add_args_and_kwargs, None, None, None, None, "sample_method"),
    ],
)
def test_class_method_test_inputs_errors(
    cls, args_to_cls, kwargs_to_cls, args_to_method, kwargs_to_method, method_name
):
    """
    Test that ``ClassMethodTestInputs`` raises exceptions upon
    instantiation, as appropriate.  This includes cases where
    the method name is not a string or refers to an attribute,
    and where keyword arguments are not a dictionary.
    """
    with pytest.raises(InvalidTestError):
        ClassMethodTestInputs(
            cls,
            method_name,
            args_to_cls,
            kwargs_to_cls,
            args_to_method,
            kwargs_to_method,
        )


def test_validate_args_value():
    """
    Test that passing a value into _validate_args returns a tuple
    containing that value.
    """
    assert (sample_value,) == _validate_args(sample_value)


@pytest.mark.parametrize("input", [sample_args, sample_args_list])
def test_validate_args_collections(input):
    """
    Test that passing a tuple or list into _validate_args returns the
    same tuple or list.
    """
    assert input == _validate_args(input)


def test_validate_kwargs():
    """
    Test that _validate_kwargs returns the input dictionary if it is a
    valid keyword arguments dictionary.
    """
    assert sample_kwargs == _validate_kwargs(sample_kwargs)


@pytest.mark.parametrize(
    "bad_kwargs, exception",
    [(bad_kwargs_wrong_type, TypeError), (bad_kwargs_nonstring_key, ValueError)],
)
def test_validate_kwargs_errors(bad_kwargs, exception):
    """
    Test that _validate_kwargs raises appropriate exceptions when
    passed invalid keyword arguments.
    """
    with pytest.raises(exception):
        _validate_kwargs(bad_kwargs)
