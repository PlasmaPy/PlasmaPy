"""Classes and functions that validate, store, and call inputs for test cases."""

__all__ = [
    "AbstractTestInputs",
    "FunctionTestInputs",
    "GenericClassTestInputs",
    "ClassAttributeTestInputs",
    "ClassMethodTestInputs",
]

import inspect
from abc import ABC, abstractmethod
from typing import Any, Callable, Dict, List, NoReturn, Optional, Tuple, Union

from plasmapy.tests.helpers.exceptions import InvalidTestError
from plasmapy.utils.formatting.formatting import (
    _object_name,
    call_string,
    class_attribute_call_string,
    class_method_call_string,
)


def _validate_args(args: Any) -> Union[Tuple, List]:
    """
    If the argument is a `tuple` or `list`, then return the argument.
    Otherwise, return a `tuple` containing the argument.  If the input
    is `None`, then this function will return an empty `tuple`.
    """

    if isinstance(args, (tuple, list)):
        return args
    elif args is None:
        return tuple()
    else:
        return (args,)


def _validate_kwargs(kwargs: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Check that the argument is a valid `dict` to represent keyword
    arguments to be provided to some function.  If the input is `None`,
    then this function will return an empty `dict`.  This function does
    not check that

    Raises
    ------
    TypeError
        If the input is not a `dict`.

    ValueError
        If the `dict` contains a key that is not a string.
    """

    if kwargs is None:
        return dict()
    if not isinstance(kwargs, dict):
        raise TypeError("kwargs must be represented as a dictionary.")
    for key in kwargs.keys():
        if not isinstance(key, str):
            raise ValueError(f"Invalid key for keyword arguments dict: {key}")
    return kwargs


class AbstractTestInputs(ABC):
    @abstractmethod
    def call(self) -> NoReturn:
        pass

    @property
    @abstractmethod
    def call_string(self) -> str:
        pass


class FunctionTestInputs(AbstractTestInputs):
    """
    A class to contain a function to be tested, along with the
    positional and keyword arguments to be provided to the function
    during the test.

    Parameters
    ----------
    function
        The function to be tested.

    args : optional
        The positional arguments to be provided to the function that is
        being tested.  If this is a `tuple` or `list`, then the
        arguments will be each of the items in the collection. If there
        is only one positional argument, then it may be inputted as is
        without being put into a `tuple` or `list`.

    kwargs : `dict`, optional
        The keyword arguments to be provided to the function that is
        being tested. If provided, the keys of this `dict` must be
        strings.
    """

    def __init__(
        self,
        function: Callable,
        args: Optional[Union[Tuple, List, Any]] = None,
        kwargs: Optional[Dict[str, Any]] = None,
    ):

        try:
            self._info = dict()
            self.function = function
            self.args = args
            self.kwargs = kwargs
        except Exception as exc:
            # TODO: Make sure that repr is what we want here and below
            raise InvalidTestError(
                f"Unable to instantiate with function = {repr(function)}, "
                f"args = {repr(args)}, and kwargs = {repr(kwargs)}."
            ) from exc

    @property
    def function(self) -> Callable:
        """The function to be tested."""

        __tracebackhide__ = False
        return self._info["function"]

    @property
    def args(self) -> Union[Tuple, List]:
        """The positional arguments to be provided to the function that is being tested."""

        return self._info["args"]

    @property
    def kwargs(self) -> Dict[str, Any]:
        """The keyword arguments to be provided to the function that is being tested."""

        return self._info["kwargs"]

    @function.setter
    def function(self, provided_function: Callable):

        if not callable(provided_function):
            function_name = _object_name(provided_function)
            raise TypeError(f"The provided function ({function_name}) is not callable.")
        self._info["function"] = provided_function

    @args.setter
    def args(self, provided_args: Any):

        self._info["args"] = _validate_args(provided_args)

    @kwargs.setter
    def kwargs(self, provided_kwargs: Optional[Dict[str, Any]]):

        self._info["kwargs"] = _validate_kwargs(provided_kwargs)

    def call(self) -> Any:
        """
        Call the function to be tested with the provided positional and
        keyword arguments.
        """

        __tracebackhide__ = True
        return self.function(*self.args, **self.kwargs)

    @property
    def call_string(self) -> str:
        """
        Return the string corresponding to calling a function or class
        method or accessing a class attribute.
        """

        return call_string(self.function, args=self.args, kwargs=self.kwargs)


class GenericClassTestInputs(AbstractTestInputs):

    _info = dict()

    @property
    def cls(self):
        """The class to be tested."""

        return self._info["cls"]

    @property
    def cls_args(self) -> Union[Tuple, List]:
        """The positional arguments to be passed to the class during instantiation."""

        return self._info["cls_args"]

    @property
    def cls_kwargs(self) -> Dict[str, Any]:
        """The keyword arguments to be passed to the class during instantiation."""

        return self._info["cls_kwargs"]

    @cls.setter
    def cls(self, provided_cls):

        if inspect.isclass(provided_cls):
            self._info["cls"] = provided_cls
        else:
            raise TypeError("Expecting a class.")

    @cls_args.setter
    def cls_args(self, provided_cls_args: Any):

        self._info["cls_args"] = _validate_args(provided_cls_args)

    @cls_kwargs.setter
    def cls_kwargs(self, provided_cls_kwargs: Optional[Dict[str, Any]]):

        self._info["cls_kwargs"] = _validate_kwargs(provided_cls_kwargs)


class ClassAttributeTestInputs(GenericClassTestInputs):
    """
    Contains a class to be tested, positional and keyword arguments to
    be provided to the class during instantiation, the name of the
    class method to be tested, and positional and keyword arguments to
    be provided to the method during the test.

    Parameters
    ----------
    cls
        The class containing the method to be tested.

    attribute : str
        The name of the attribute contained in ``cls`` to be tested.

    cls_args : optional
        The positional arguments to be provided to the class that is
        being tested during instantiation. If this is a `tuple` or
        `list`, then the arguments will be each of the items in the
        collection. If there is only one positional argument, then it
        may be inputted as is without being put into a `tuple` or `list`.

    cls_kwargs : `dict`, optional
        The keyword arguments to be provided to the class that is
        being tested during instantiation. If provided, the keys of
        this `dict` must be strings.

    Raises
    ------
    ~plasmapy.utils.pytest_helpers.InvalidTestError
        If this class cannot be instantiated.
    """

    def __init__(
        self,
        cls,
        attribute: str,
        cls_args: Any = None,
        cls_kwargs: Optional[Dict[str, Any]] = None,
    ):

        try:
            self._info = {}
            self.cls = cls
            self.cls_args = cls_args
            self.cls_kwargs = cls_kwargs
            self.attribute = attribute
        except Exception as exc:
            # TODO: Make sure that the input variables are being represented correctly here
            raise InvalidTestError(
                f"Unable to instantiate ClassAttributeTestInputs for "
                f"class {cls.__name__} with cls_args = {repr(cls_args)}, "
                f"cls_kwargs = {repr(cls_kwargs)}, and attribute = "
                f"{repr(attribute)}."
            ) from exc

    @property
    def attribute(self) -> str:
        """The name of the attribute to be tested."""

        return self._info["attribute"]

    @attribute.setter
    def attribute(self, attribute_name: str):

        if not isinstance(attribute_name, str):
            raise TypeError("Expecting the name of a class attribute as a string.")
        if not hasattr(self.cls, attribute_name):
            raise ValueError(f"No attribute named {attribute_name}.")
        else:
            self._info["attribute"] = attribute_name

    def call(self) -> Any:
        """Instantiate the class and access the attribute."""

        __tracebackhide__ = True
        instance = self.cls(*self.cls_args, **self.cls_kwargs)
        return getattr(instance, self.attribute)

    @property
    def call_string(self) -> str:
        """
        Return a string that is designed to reproduce class instantiation
        and accessing the class attribute.
        """

        return class_attribute_call_string(
            cls=self.cls,
            attr=self.attribute,
            cls_args=self.cls_args,
            cls_kwargs=self.cls_kwargs,
        )


class ClassMethodTestInputs(GenericClassTestInputs):
    """
    Contains a class to be tested, positional and keyword arguments to
    be provided to the class during instantiation, the name of the
    class method to be tested, and positional and keyword arguments to
    be provided to the method during the test.

    Parameters
    ----------
    cls
        The class containing the method to be tested.

    method : str
        The name of the method contained in ``cls`` to be tested.

    cls_args : optional
        The positional arguments to be provided to the class that is
        being tested during instantiation.  If this is a `tuple` or
        `list`, then the arguments will be each of the items in the
        collection. If there is only one positional argument, then it
        may be inputted as is without being put into a `tuple` or `list`.

    cls_kwargs : `dict`, optional
        The keyword arguments to be provided to the function that is
        being tested. If provided, the keys of this `dict` must be
        strings.

    method_args : optional
        The positional arguments to be provided to the method that is
        being tested.  If this is a `tuple` or `list`, then the
        arguments will be each of the items in the collection. If there
        is only one positional argument, then it may be inputted as is
        without being put into a `tuple` or `list`.

    method_kwargs : `dict`, optional
        The keyword arguments to be provided to the method that is being
        tested. If provided, the keys of this `dict` must be strings.

    Raises
    ------
    ~plasmapy.utils.pytest_helpers.InvalidTestError
        If this class cannot be instantiated.

    """

    def __init__(
        self,
        cls,
        method: str,
        cls_args=None,
        cls_kwargs: Optional[Dict[str, Any]] = None,
        method_args=None,
        method_kwargs: Optional[Dict[str, Any]] = None,
    ):

        try:
            self._info = {}
            self.cls = cls
            self.cls_args = cls_args
            self.cls_kwargs = cls_kwargs
            self.method = method
            self.method_args = method_args
            self.method_kwargs = method_kwargs
        except Exception as exc:
            # TODO: Test that error message is representing things correctly
            raise InvalidTestError(
                f"Unable to instantiate ClassMethodTestInputs for class "
                f"{cls.__name__} with cls_args = {cls_args}, cls_kwargs = "
                f"{cls_kwargs}, method = {method}, method_args = "
                f"{method_args}, and method_kwargs = {method_kwargs}."
            ) from exc

    @property
    def method(self) -> str:
        """The name of the method to be tested."""

        return self._info["method"]

    @property
    def cls_args(self) -> Union[Tuple, List]:
        """The positional arguments to be provided to the class upon instantiation."""

        return self._info["cls_args"]

    @property
    def cls_kwargs(self) -> Dict[str, Any]:
        """The keyword arguments to be provided to the class upon instantiation."""

        return self._info["cls_kwargs"]

    @property
    def method_args(self) -> Union[Tuple, List]:
        """The positional arguments to be provided to the method being tested."""

        return self._info["method_args"]

    @property
    def method_kwargs(self) -> Dict[str, Any]:
        """The keyword arguments to be provided to the method being tested."""

        return self._info["method_kwargs"]

    @method.setter
    def method(self, method_name: str):

        if not isinstance(method_name, str):
            raise TypeError("Expecting the name of a method as a string.")
        elif not hasattr(self.cls, method_name):
            raise ValueError(
                f"No method named {method_name} in class {self.cls.__name__}."
            )

        actual_method = getattr(self.cls, method_name)

        if callable(actual_method):
            self._info["method"] = method_name
        else:
            raise ValueError(
                f"Expecting the attribute {method_name} in "
                f"{self.cls.__name__} to be callable."
            )

    @cls_args.setter
    def cls_args(self, provided_args):

        self._info["cls_args"] = _validate_args(provided_args)

    @cls_kwargs.setter
    def cls_kwargs(self, provided_kwargs: Optional[Dict[str, Any]]):

        self._info["cls_kwargs"] = _validate_kwargs(provided_kwargs)

    @method_args.setter
    def method_args(self, provided_args):

        self._info["method_args"] = _validate_args(provided_args)

    @method_kwargs.setter
    def method_kwargs(self, provided_kwargs: Optional[Dict[str, Any]]):

        self._info["method_kwargs"] = _validate_kwargs(provided_kwargs)

    def call(self) -> Any:
        """Instantiate the class and call the appropriate method."""

        __tracebackhide__ = True
        instance = self.cls(*self.cls_args, **self.cls_kwargs)
        method = getattr(instance, self.method)
        return method(*self.method_args, **self.method_kwargs)

    @property
    def call_string(self) -> str:
        """
        Return a string that is designed to reproduce class instantiation
        and accessing the class attribute.
        """

        return class_method_call_string(
            cls=self.cls,
            method=self.method,
            cls_args=self.cls_args,
            cls_kwargs=self.cls_kwargs,
            method_args=self.method_args,
            method_kwargs=self.method_kwargs,
        )
