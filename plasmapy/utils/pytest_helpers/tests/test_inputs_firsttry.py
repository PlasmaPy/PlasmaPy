import pytest

from plasmapy.utils.pytest_helpers.inputs import (
    ClassMethodTestInputs,
    ClassAttributeTestInputs,
    FunctionTestInputs,
    _validate_args,
    _validate_kwargs
)

from plasmapy.utils.pytest_helpers.exceptions import InvalidTestError

sample_args_tuple = (2, 3)
sample_args_list = [5, 7, 11]
sample_kwargs = {'kwarg1': 13, 'kwarg2': 17}
sample_bad_kwargs_nonstring_key = {'x': 19, 23: 'y', '29': 'z'}
sample_bad_kwargs_wrong_type = [23, 27]
sample_method_args = (31, 37)
sample_method_kwargs = {'method_kwarg1': 41, 'method_kwarg2': 43}
sample_value = 47

bad_kwargs_list = [sample_bad_kwargs_wrong_type, sample_bad_kwargs_nonstring_key]


def sample_function(arg1, arg2, kwarg1=None, kwarg2=None):
    return arg1 + arg2 + kwarg1 + kwarg2


def sample_function_no_args():
    return sample_value


class SampleClass:
    def __init__(self, arg1, arg2, kwarg1=None, kwarg2=None):
        self.arg1 = arg1
        self.arg2 = arg2
        self.kwarg1 = kwarg1
        self.kwarg2 = kwarg2

    def sample_method(self, method_arg1, method_arg2, method_kwarg1=None, method_kwarg2=None):
        return method_arg1 + method_arg2 + method_kwarg1 + method_kwarg2

    @property
    def sample_attribute(self):
        return sample_value


@pytest.fixture
def sample_test_function_instance():
    return FunctionTestInputs(sample_function, sample_args_tuple, sample_kwargs)



def test_of_function_test_inputs():
    instance = FunctionTestInputs()


@pytest.fixture
def sample_test_function_instance_no_args_or_kwargs():
    return FunctionTestInputs(sample_function_no_args)


@pytest.fixture
def sample_test_attribute_instance():
    return ClassAttributeTestInputs(
        SampleClass,
        'sample_attribute',
        sample_args_tuple,
        sample_kwargs,
    )


@pytest.fixture
def sample_test_attribute_instance_no_args_or_kwargs():
    return ClassAttributeTestInputs(SampleClass, 'sample_attribute')


@pytest.fixture
def sample_test_method_instance():
    return ClassMethodTestInputs(
        SampleClass,
        'sample_method',
        sample_args_tuple,
        sample_kwargs,
        sample_method_args,
        sample_method_kwargs,
    )


@pytest.fixture
def sample_test_method_instance_no_args_or_kwargs():
    return ClassMethodTestInputs(SampleClass, 'sample_method')


def test_func_inputs_function(sample_test_function_instance):
    """Test that FunctionTestInputs correctly stores the function to be tested."""
    if sample_function is not sample_test_function_instance.function:
        pytest.fail("FunctionTestInputs.function is not the inputted function.")


def test_func_inputs_args(sample_test_function_instance):
    """Test that FunctionTestInputs correctly stores args provided in the form of a tuple."""
    if sample_args_tuple != sample_test_function_instance.args:
        pytest.fail("FunctionTestInputs.args does not contain the provided arguments.")


def test_func_inputs_kwargs(sample_test_function_instance):
    """Test that FunctionTestInputs correctly stores kwargs."""
    if sample_kwargs != sample_test_function_instance.kwargs:
        pytest.fail("FunctionTestInputs.kwargs does not contain the provided keyword arguments.")


def test_func_inputs_using_a_value_as_args():
    """
    Test that FunctionTestInputs correctly treats a value inputted as
    the sole positional argument by containing it within a tuple.
    """
    instance = FunctionTestInputs(sample_function, sample_value)
    if instance.args != (sample_value,):
        pytest.fail("FunctionTestInputs.args is not a tuple containing the correct value.")


def test_func_inputs_no_args(sample_test_function_instance_no_args_or_kwargs):
    """
    Test that an instance of FunctionTestInputs created without any
    positional arguments has an empty tuple for the args attribute.
    """
    if sample_test_function_instance_no_args_or_kwargs.args != ():
        pytest.fail(
            "The args attribute should be an empty tuple if no args are "
            "being provided for the function call."
        )


def test_func_inputs_no_kwargs(sample_test_function_instance_no_args_or_kwargs):
    """
    Test that an instance of FunctionTestInputs created without any
    keyword arguments has an empty dictionary for the kwargs attribute.
    """
    if sample_test_function_instance_no_args_or_kwargs.kwargs != {}:
        pytest.fail(
            "The kwargs attribute should be an empty dictionary if no "
            "kwargs are being provided for the function call."
        )


def test_func_inputs_call(sample_test_function_instance):
    """
    Test that the call method on FunctionTestInputs works correctly
    when args and kwargs are provided.
    """
    result = sample_test_function_instance.call()
    expected = sample_function(*sample_args_tuple, **sample_kwargs)
    if result != expected:
        pytest.fail(
            "The call method on FunctionTestInputs is not returning the "
            "expected result."
        )


def test_func_inputs_call_no_args_or_kwargs(sample_test_function_instance_no_args_or_kwargs):
    """
    Test that the call method on FunctionTestInputs works correctly
    when args and kwargs are not provided.
    """
    result = sample_test_function_instance_no_args_or_kwargs.call()
    expected = sample_function_no_args()
    if result != expected:
        pytest.fail(
            f"The call method for a FunctionTestInputs with no arguments "
            f"returned {result} instead of the expected value of {expected}."
        )


@pytest.mark.parametrize("bad_function", [SampleClass, '', SampleClass.sample_method])
def test_func_inputs_bad_function(bad_function):
    with pytest.raises(InvalidTestError):
        FunctionTestInputs(bad_function)


@pytest.mark.parametrize("bad_kwargs", bad_kwargs_list)
def test_func_inputs_bad_kwargs(bad_kwargs):
    """Test that an appropriate exception is raised for bad inputs to FunctionTestInputs."""
    with pytest.raises(InvalidTestError):
        FunctionTestInputs(sample_function, kwargs=bad_kwargs)


def test_cls_inputs_attr_cls(sample_test_attribute_instance):
    if sample_test_attribute_instance.cls is not SampleClass:
        pytest.fail(
            "ClassAttributeTestInputs is not storing the correct class "
            "in its cls attribute."
        )


def test_cls_inputs_attr_args(sample_test_attribute_instance):
    if sample_test_attribute_instance.cls_args != sample_args_tuple:
        pytest.fail(
            "The args attribute of ClassAttributeTestInputs does not "
            "contain a tuple containing the supplied args."
        )


def test_cls_inputs_attr_kwargs(sample_test_attribute_instance):
    if sample_test_attribute_instance.cls_kwargs is not sample_kwargs:
        pytest.fail(
            "The kwargs attribute of ClassAttributeTestInputs does not "
            "contain a dict of the supplied kwargs."
        )


def test_cls_inputs_attr_attribute(sample_test_attribute_instance):
    if sample_test_attribute_instance.attribute != 'sample_attribute':
        pytest.fail(
            "The attribute named 'attribute' of ClassAttributeTestInputs "
            "is not the name of the attribute that was passed in during "
            "instantiation."
        )


def test_cls_inputs_attr_call(sample_test_attribute_instance):
    result = sample_test_attribute_instance.call()
    if result != sample_value:
        pytest.fail()






def test_validate_args_value():
    """
    Test that passing a value into _validate_args returns a tuple
    containing that value.
    """
    assert (sample_value,) == _validate_args(sample_value)


@pytest.mark.parametrize("input", [sample_args_tuple, sample_args_list])
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
    "bad_kwargs, exception", [
        (sample_bad_kwargs_wrong_type, TypeError),
        (sample_bad_kwargs_nonstring_key, ValueError),
    ])
def test_validate_kwargs_errors(bad_kwargs, exception):
    """
    Test that _validate_kwargs raises appropriate exceptions when
    passed invalid keyword arguments.
    """
    with pytest.raises(exception):
        _validate_kwargs(bad_kwargs)

