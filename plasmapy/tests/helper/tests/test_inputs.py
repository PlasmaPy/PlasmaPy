import pytest

from plasmapy.tests.helper.inputs import (
    ClassMethodTestInputs,
    ClassAttributeTestInputs,
    FunctionTestInputs,
    _validate_args,
    _validate_kwargs,
)

from plasmapy.tests.helper.exceptions import InvalidTestError

sample_args = (2, 3)
sample_args_list = [5, 7]
sample_kwargs = {"kwarg1": 13, "kwarg2": 17}
bad_kwargs_nonstring_key = {"kwarg1": 19, 23: 29}
bad_kwargs_wrong_type = [31, 37]
sample_method_args = (41, 43)
sample_method_kwargs = {"method_kwarg1": 47, "method_kwarg2": 53}
sample_value = 59


def sample_function(arg1, arg2, kwarg1=None, kwarg2=None):
    """Sum up the positional and keyword arguments."""
    return sum([arg1, arg2, kwarg1, kwarg2])


def sample_function_one_arg(arg):
    """Double the argument."""
    return arg * 2


def sample_function_no_args():
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
        (sample_function, sample_args, sample_kwargs, "function", sample_function),
        (sample_function, sample_args, sample_kwargs, "args", sample_args),
        (sample_function, sample_args_list, sample_kwargs, "args", sample_args_list),
        (sample_function, sample_args, sample_kwargs, "kwargs", sample_kwargs),
        (sample_function_one_arg, sample_value, None, "args", (sample_value,)),
        (sample_function_no_args, None, None, "args", ()),
        (sample_function_no_args, None, None, "kwargs", {}),
        (
            sample_function_no_args,
            None,
            None,
            "call_string",
            "sample_function_no_args()",
        ),
        (
            sample_function,
            sample_args,
            sample_kwargs,
            "call_string",
            "sample_function(2, 3, kwarg1=13, kwarg2=17)",
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
            sample_function,
            sample_args,
            sample_kwargs,
            sample_function(*sample_args, **sample_kwargs),
        ),
        (
            sample_function_one_arg,
            sample_value,
            None,
            sample_function_one_arg(sample_value),
        ),
        (sample_function_no_args, None, None, sample_value),
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
        (sample_function, sample_args, bad_kwargs_nonstring_key),
        (sample_function, sample_args, bad_kwargs_wrong_type),
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
        (SampleClass, sample_args, sample_kwargs, "cls_args", sample_args),
        (SampleClass, sample_args_list, sample_kwargs, "cls_args", sample_args_list),
        (SampleClass, None, None, "cls_args", ()),
        (SampleClass, 42, None, "cls_args", (42,)),
        (SampleClass, sample_args, sample_kwargs, "cls_kwargs", sample_kwargs),
        (SampleClass, sample_args, None, "cls_kwargs", {}),
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
    sample_method_args,
    sample_method_kwargs,
)


@pytest.mark.parametrize(
    "cls, cls_args, cls_kwargs, method_args, method_kwargs, attribute_being_tested, expected",
    [
        (*common_inputs, "cls", SampleClass),
        (*common_inputs, "cls_args", sample_args),
        (*common_inputs, "cls_kwargs", sample_kwargs),
        (*common_inputs, "method_args", sample_method_args),
        (*common_inputs, "method_kwargs", sample_method_kwargs),
        (*common_inputs, "method", "sample_method"),
        (SampleClassNoArgs, None, None, None, None, "cls", SampleClassNoArgs),
        (SampleClassNoArgs, None, None, None, None, "cls_args", ()),
        (SampleClassNoArgs, None, None, None, None, "cls_kwargs", {}),
        (SampleClassNoArgs, None, None, None, None, "method", "sample_method"),
        (SampleClassNoArgs, None, None, None, None, "method_args", ()),
        (SampleClassNoArgs, None, None, None, None, "method_kwargs", {}),
        (SampleClassNoArgs, 42, None, None, None, "cls_args", (42,)),
        (SampleClassNoArgs, None, None, 42, None, "method_args", (42,)),
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
    cls_args,
    cls_kwargs,
    method_args,
    method_kwargs,
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

    cls_args
        Positional arguments to be used when instantiating ``cls``

    cls_kwargs : `dict` or `None`
        Keyword arguments to be used when instantiating ``cls``

    method_args
        Positional arguments to be used when calling the method
        ``sample_method`` in an instance of ``cls``

    method_kwargs : `dict` or `None`
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
        cls, method_name, cls_args, cls_kwargs, method_args, method_kwargs
    )
    value_of_attribute = getattr(instance, attribute_being_tested)
    if value_of_attribute != expected:
        pytest.fail(
            f"Expecting {expected} but got {value_of_attribute} "
            f"for attribute {attribute_being_tested}."
        )


@pytest.mark.parametrize(
    "cls, cls_args, cls_kwargs, method_args, method_kwargs, expected",
    [(*common_inputs, sum(sample_method_args) + sum(sample_method_kwargs.values()))],
)
def test_class_method_test_inputs_call(
    cls, cls_args, cls_kwargs, method_args, method_kwargs, expected
):
    """
    Test that the ``call`` method on ClassMethodTestInputs returns
    the expected values for a sample class with a method named
    ``sample_method``.
    """
    method_name = "sample_method"
    instance = ClassMethodTestInputs(
        cls, method_name, cls_args, cls_kwargs, method_args, method_kwargs
    )
    value_of_call = instance.call()
    if value_of_call != expected:
        pytest.fail(
            f"The call method on a ClassMethodTestInputs instance is "
            f"returning {value_of_call}, rather than the expected value "
            f"of {expected}.  This case has cls = {cls.__name__}, "
            f"cls_args = {cls_args}, cls_kwargs = {cls_kwargs}, "
            f"method_args = {method_args}, and method_kwargs = "
            f"{method_kwargs}."
        )


@pytest.mark.parametrize(
    "cls, cls_args, cls_kwargs, method_args, method_kwargs, method_name",
    [
        (*common_inputs, 42),
        (*common_inputs, "sample_attribute"),
        (SampleClass, None, (), None, None, "sample_method"),
        (SampleClass, None, None, None, (), "sample_method"),
        (sample_function, None, None, None, None, "sample_method"),
    ],
)
def test_class_method_test_inputs_errors(
    cls, cls_args, cls_kwargs, method_args, method_kwargs, method_name
):
    """
    Test that ``ClassMethodTestInputs`` raises exceptions upon
    instantiation, as appropriate.  This includes cases where
    the method name is not a string or refers to an attribute,
    and where keyword arguments are not a dictionary.
    """
    with pytest.raises(InvalidTestError):
        ClassMethodTestInputs(
            cls, method_name, cls_args, cls_kwargs, method_args, method_kwargs
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
