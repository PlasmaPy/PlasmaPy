import inspect
from typing import Callable, Optional, Dict

from .expected import ExpectedTestOutcome
from .inputs import FunctionTestInputs, ClassMethodTestInputs, ClassAttributeTestInputs


#def run_test(
#    expected,
#    func_or_class: Callable,
#    args=None,
#    kwargs: Optional[Dict[str, Any]] = None,
#    method=None,
#    method_args=None,
#    method_kwargs: Optional[Dict[str, Any]] = None,
#    attribute: Optional[str] = None,
#    rtol = 0.0,
#    atol = 0.0,
#) -> None:

    # The goal here is to keep only the highest level functionality,
    # and delegate all of the lower level stuff to more specialized
    # functions and classes.  The class ExpectedTestOutcome will take
    # care of processing the test inputs.  The various TestInputs
    # classes will process the different args & kwargs, and each have
    # a method that will call the thing that needs to be called (be
    # it a function, class method, or class attribute).

#    __tracebackhide__ = True
#    expected_test_outcome = ExpectedTestOutcome(expected)
#    test_inputs = _test_inputs_factory(
#        func_or_class, args, kwargs, method, method_args, method_kwargs, attribute
#    )
#    _perform_test(expected_test_outcome, test_inputs, rtol, atol)


#def _test_inputs_factory(
#        func_or_class, args, kwargs, method, method_args, method_kwargs, attribute
#):
#    """Select the class """
#    __tracebackhide__ = True
#    if inspect.isfunction(func_or_class):
#        return FunctionTestInputs(func_or_class, args, kwargs)
#    elif inspect.isclass(func_or_class):
#
#    else:
#        raise TypeError("Expecting a function or a class.")


def test_a_function(
        expected,
        func,
        args,
        kwargs: Optional[Dict[str, Any]],
) -> None:
    __tracebackhide__ = True
    expected_test_outcome = ExpectedTestOutcome(expected)
    test_inputs = FunctionTestInputs(function, args, kwargs)
    _perform_test(expected_test_outcome, test_inputs)


def test_a_class_method(
        expected,
        cls,
        cls_args = None,
        cls_kwargs: Optional[Dict[str, Any]] = None,
        method = None,
        method_args = None,
        method_kwargs: Optional[Dict[str, Any]] = None,
) -> None:
    __tracebackhide__ = True
    expected_test_outcome = ExpectedTestOutcome(expected)
    test_inputs = ClassMethodTestInputs(
        cls, cls_args, cls_kwargs, method, method_args, method_kwargs,
    )
    _perform_test(expected_test_outcome, )


def test_a_class_attribute(
        expected,
        cls,
        cls_args = None,
        cls_kwargs: Optional[Dict[str, Any]] = None,
        attribute: str = None,
) -> None:
    __tracebackhide__ = True
    expected_test_outcome = expected
    test_inputs = ClassAttributeTestInputs(
        cls, cls_args, cls_kwargs, attribute
    )


def _perform_test(expected, test_inputs, rtol, atol):
    __tracebackhide__ = True
    if expected.expecting_an_exception:
        _test_for_exception(expected, test_inputs)
    elif expected.expecting_a_warning or :
        _test_for_warning_and_or_result(expected, test_inputs, rtol, atol)


def _test_for_exception(expected, test_inputs):
    __tracebackhide__ = True


def _test_for_warning_and_or_result(expected, test_inputs, rtol, atol):
    __tracebackhide__ = True



