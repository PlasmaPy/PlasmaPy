import abc
from typing import Tuple, List, Dict, Union, Callable
import inspect


class TestInputs:
    def __init__(
            self,
            cls_or_func,
            args,
            kwargs,
            attribute,
            method,
            method_args,
            method_kwargs,
    ):
        self.factory = choose_factory





def _validate_args(args):
    """
    If the argument is a `tuple` or `list`, then return the argument.
    Otherwise, return a `tuple` containing the argument.
    """
    return args if isinstance(args, (tuple, list)) else (args,)


def _validate_kwargs(kwargs: Dict[str]) -> Dict[str]:
    """
    Check that the provided keyword arguments are valid, and if so
    """

    if not isinstance(kwargs, dict):
        raise TypeError("kwargs must be represented as a dictionary.")
    for key in kwargs.keys():
        if not isinstance(key, str):
            raise ValueError(f"Invalid key for keyword arguments dict: {key}")
    return kwargs


@abc.ABC
class AbstractTestInputs:
    @abc.abstractmethod
    def call(self):
        ...


class FunctionTestInputs(AbstractTestInputs):
    def __init__(self, function, args=tuple(), kwargs = {}):
        self._info = {}
        self.function = function
        self.args = args
        self.kwargs = kwargs

    @property
    def function(self) -> callable:
        __tracebackhide__ = False
        return self._info['function']

    @property
    def args(self):
        return self._info['args']

    @property
    def kwargs(self):
        return self._info['kwargs']

    @function.setter
    def function(self, provided_function):
        if not inspect.isfunction(provided_function):
            raise TypeError("The provided function is not a function.")

    @args.setter
    def args(self, provided_args):
        is_tuple_or_list = isinstance(provided_args, (tuple, list))
        self._args = provided_args if is_tuple_or_list else (provided_args,)

    @kwargs.setter
    def kwargs(self, provided_kwargs):
        if not isinstance(kwargs, dict):
            raise TypeError("Keyword arguments should be a dictionary with strings as keys")

    def __call__(self):
        __tracebackhide__ = True
        return self.function(*self.args, **self.kwargs)


class AbstractClassTestInputs(AbstractTestInputs):

    @property
    def cls(self):
        return self._info['cls']

    @property
    def cls_args(self):
        return self._info['cls_args']

    @property
    def cls_kwargs(self):
        return self._info['cls_kwargs']

    @cls.setter
    def cls(self, provided_cls):
        if inspect.isclass(provided_cls):
            self.cls = provided_cls
        else:
            raise TypeError("Expecting a class.")

    @cls_args.setter
    def cls_args(self, provided_cls_args):
        self._info['cls_args'] = _validate_args(provided_cls_args)

    @cls_kwargs.setter
    def cls_kwargs(self, provided_cls_kwargs):
        self._info['cls_kwargs'] = _validate_kwargs(provided_cls_kwargs)


class ClassAttributeTestInputs(ClassTestInputs):
    def __init__(self, cls, cls_args, cls_kwargs, attribute):
        self._info = {}
        self.cls = cls
        self.cls_args = cls_args
        self.cls_kwargs = cls_kwargs
        self.attribute = attribute

    @property
    def attribute(self):
        return self._info['attribute']

    @attribute.setter
    def attribute(self, provided_attribute):
        if not isinstance(provided_attribute, str):
            raise TypeError("Expecting a string")
        if not hasattr(self.cls, provided_attribute):
            raise ValueError(f"No attribute named {provided_attribute}.")
        else:
            self._info['attribute'] = provided_attribute

    @property
    def call(self):
        __tracebackhide__ = True
        return self.cls(*self.cls_args, **self.cls_kwargs)


class ClassMethodTestInputs(AbstractTestInputs):
    def __init__(self, cls, cls_args, cls_kwargs, method, method_args, method_kwargs):
        self.info = {}
        self.cls = cls
        self.cls_args = cls_args
        self.cls_kwargs = cls_kwargs
        self.method = method
        self.method_args = method_args
        self.method_kwargs = method_kwargs
