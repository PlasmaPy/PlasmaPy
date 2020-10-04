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
    _test_params = NotImplemented
    _test_param_errors = NotImplemented

    @property
    @abstractmethod
    def ff_class(self):
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def func(x, *args):
        raise NotImplementedError

    def test_inheritance(self):
        assert issubclass(self.ff_class, self.abc)

    def test_iscallable(self):
        assert callable(self.ff_class())

    def test_basics(self):
        assert hasattr(self.ff_class, "_param_names")
        if self.ff_class._param_names == NotImplemented:
            pytest.fail(
                f"{self.ff_class} class attribute '_param_names' needs to "
                f" be defined as a tuple of strings representing the names of "
                f"the fit parameters."
            )

        # required attributes/methods
        for name in ("curve_fit", "func", "func_err", "root_solve"):
            assert hasattr(self.ff_class, name)

        # required properties
        for name in ("curve_fit_results", "latex_str", "params", "param_errors",
                     "param_names", "rsq"):
            assert hasattr(self.ff_class, name)
            assert isinstance(getattr(self.ff_class, name), property)

    def test_instantiation(self):
        # default
        foo = self.ff_class()

        assert isinstance(foo.param_names, tuple)
        assert len(foo.param_names) != 0
        assert all(isinstance(val, str) for val in foo.param_names)

        assert hasattr(foo, "FitParamTuple")
        assert issubclass(foo.FitParamTuple, tuple)
        for name in foo.param_names:
            assert hasattr(foo.FitParamTuple, name)

        assert foo.curve_fit_results is None
        assert foo.params is None
        assert foo.param_errors is None
        assert foo.rsq is None

        assert isinstance(foo.latex_str, str)

        # assign at instantiation
        params = [1] * len(foo.param_names)
        foo = self.ff_class(params=params, param_errors=params)
        assert foo.params == foo.FitParamTuple(*params)
        assert foo.param_errors == foo.FitParamTuple(*params)

        with pytest.raises(ValueError):
            self.ff_class(params=5)
        with pytest.raises(ValueError):
            self.ff_class(param_errors=5)

        params = [2] * len(foo.param_names)
        params[0] = "let me in"
        with pytest.raises(ValueError):
            self.ff_class(params=params)
        with pytest.raises(ValueError):
            self.ff_class(param_errors=params)

        params = [2] * len(foo.param_names)
        params += [5]
        with pytest.raises(ValueError):
            self.ff_class(params=params)
        with pytest.raises(ValueError):
            self.ff_class(param_errors=params)

    def test_param_assignment(self):
        foo = self.ff_class()

        # setting params property
        params = [2] * len(foo.param_names)
        for val in (params, foo.FitParamTuple(*params)):
            foo.params = val
            assert foo.params == foo.FitParamTuple(*params)

        with pytest.raises(ValueError):
            foo.params = 5

        params = [2] * len(foo.param_names)
        params[0] = "let me in"
        with pytest.raises(ValueError):
            foo.params = params

        params = [2] * len(foo.param_names)
        params += [5]
        with pytest.raises(ValueError):
            foo.params = params

        # setting param_errors property
        params = [2] * len(foo.param_names)
        for val in (params, foo.FitParamTuple(*params)):
            foo.param_errors = val
            assert foo.param_errors == foo.FitParamTuple(*params)

        with pytest.raises(ValueError):
            foo.param_errors = 5

        params = [2] * len(foo.param_names)
        params[0] = "let me in"
        with pytest.raises(ValueError):
            foo.param_errors = params

        params = [2] * len(foo.param_names)
        params += [5]
        with pytest.raises(ValueError):
            foo.param_errors = params

    def test_func(self):
        foo = self.ff_class()

        for x in (0, 1., np.linspace(10, 30, num=20)):
            assert np.allclose(foo.func(x, *self._test_params),
                               self.func(x, *self._test_params))

        x = [4, 5, 6]
        assert np.allclose(foo.func(x, *self._test_params),
                           self.func(np.array(x), *self._test_params))

        with pytest.raises(ValueError):
            foo.func("hello", *self._test_params)

        with pytest.raises(ValueError):
            params = list(self._test_params)
            params[0] = "hello"
            foo.func(5, *params)


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
    _test_params = (5., 4.)

    @staticmethod
    def func(x, m, b):
        return m * x + b

    def test_basics(self):
        super().test_basics()

        foo = self.ff_class()

        assert foo.param_names == ("m", "b")
        assert foo.latex_str == fr"m \, x + b"
        assert foo.__str__() == f"f(x) = m x + b"
