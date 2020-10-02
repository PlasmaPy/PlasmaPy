"""
Tests for the fitting function classes defined in `plasmapy.analysis.fit_functions`.
"""
import inspect
import numpy as np
import pytest

from abc import ABC, abstractmethod

import plasmapy.analysis.fit_functions as ffuncs


class BaseFFTests(ABC):
    abc = ffuncs.AbstractFitFunction

    @property
    @abstractmethod
    def ff_class(self):
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def func(x, *args):
        raise NotImplementedError

    def test_basics(self):
        assert hasattr(self.ff_class, "_param_names")
        if self.ff_class._param_names == NotImplemented:
            pytest.fail(
                f"{self.ff_class} class attribute '_param_names' needs to "
                f" be defined as a tuple of strings representing the names of "
                f"the fit parameters."
            )

        foo = self.ff_class()

        assert hasattr(foo, "func")
        assert hasattr(foo, "func_err")

        assert hasattr(foo, "FitParamTuple")
        assert issubclass(foo.FitParamTuple, tuple)

        assert hasattr(foo, "curve_fit_results")
        assert foo.curve_fit_results is None

        assert hasattr(foo, "param_names")
        assert isinstance(foo.param_names, tuple)
        assert all(isinstance(val, str) for val in foo.param_names)

        for attr in ("params", "param_errors"):
            assert hasattr(foo, attr)
            assert getattr(foo, attr) is None
            params = [1] * len(foo.param_names)
            setattr(foo, attr, params)
            assert getattr(foo, attr) == foo.FitParamTuple(*params)
            assert all(hasattr(getattr(foo, attr), name)
                       for name in foo.param_names)

        assert hasattr(foo, "latex_str")
        assert isinstance(foo.latex_str, str)

        assert hasattr(foo, "rsq")
        assert foo.rsq is None

    def test_inheritance(self):
        assert issubclass(self.ff_class, self.abc)

    def test_iscallable(self):
        assert callable(self.ff_class())


# class TestAbstractFitFunction(BaseFFTests):
#     @staticmethod
#     def func(x, a, b, c):
#         return a * x ** 2 + b * x + c
#
#     class FooFitFunc(ffuncs.AbstractFitFunction):
#         _param_names = ("a", "b", "c")
#
#         @property
#         def latex_str(self) -> str:
#             return fr"a \, x^2 + b \, x + c"
#
#         @staticmethod
#         def func(x, a, b, c):
#             return a * x ** 2 + b * x + c
#
#         def func_err(self, x, x_err=None, rety=False):
#             a, b, c = self.parameters
#             a_err, b_err, c_err = self.parameters_err
#
#             a_term = ((x ** 2) * a_err) ** 2
#             b_term = (x * b_err) ** 2
#             c_term = c_err ** 2
#             err = a_term + b_term + c_term
#
#             if x_err is not None:
#                 x_term = ((2 * a * x + b) * x_err) ** 2
#                 err += x_term
#             err = np.sqrt(err)
#
#             if rety:
#                 y = self.func(x, a, b, c)
#                 return err, y
#
#             return err
#
#     ff_class = FooFitFunc
#
#     def test_basics(self):
#         pass


class TestFFLinear(BaseFFTests):
    """
    Tests for fit function class `plasmapy.analysis.fit_functions.Linear`.
    """
    ff_class = ffuncs.Linear

    @staticmethod
    def func(x, m, b):
        return m * x + b

    def test_basics(self):
        super().test_basics()

        foo = self.ff_class()

        assert foo.param_names == ("m", "b")
        assert foo.latex_str == fr"m \, x + b"
        assert foo.__str__() == f"f(x) = m x + b"
