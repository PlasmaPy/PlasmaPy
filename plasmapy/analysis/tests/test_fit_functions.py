"""
Tests for the fitting function classes defined in `plasmapy.analysis.fit_functions`.
"""
import inspect
import numpy as np
import pytest

from abc import ABC, abstractmethod

import plasmapy.analysis.fit_functions as ffuncs


class TestAbstractFitFunction:
    """
    Tests for fit function class `plasmapy.analysis.fit_functions.AbstractFitFunction`.
    """

    ff_class = ffuncs.AbstractFitFunction

    def test_is_abs(self):
        assert issubclass(self.ff_class, ABC)

    def test_methods(self):
        # required attributes/methods
        for name in ("curve_fit", "func", "func_err", "root_solve", "__call__"):
            assert hasattr(self.ff_class, name)

        # required properties
        for name in (
            "curve_fit_results",
            "latex_str",
            "params",
            "param_errors",
            "param_names",
            "rsq",
        ):
            assert hasattr(self.ff_class, name)
            assert isinstance(getattr(self.ff_class, name), property)

    def test_abstractmethods(self):
        # abstract methods
        for name in ("__str__", "func", "func_err", "latex_str"):
            assert name in self.ff_class.__abstractmethods__


class BaseFFTests(ABC):
    abc = ffuncs.AbstractFitFunction
    _test_params = NotImplemented  # type: tuple
    _test_param_errors = NotImplemented  # type: tuple

    @property
    @abstractmethod
    def ff_class(self):
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def func(x, *args):
        raise NotImplementedError

    @abstractmethod
    def func_err(self, x, params, param_errors, x_err=None):
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
        for name in (
            "curve_fit_results",
            "latex_str",
            "params",
            "param_errors",
            "param_names",
            "rsq",
        ):
            assert hasattr(self.ff_class, name)
            assert isinstance(getattr(self.ff_class, name), property)

        foo = self.ff_class()
        assert foo.__repr__() == f"{foo.__str__()} {foo.__class__}"

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

        for x in (0, 1.0, np.linspace(10, 30, num=20)):
            assert np.allclose(
                foo.func(x, *self._test_params), self.func(x, *self._test_params)
            )

        x = [4, 5, 6]
        assert np.allclose(
            foo.func(x, *self._test_params), self.func(np.array(x), *self._test_params)
        )

        with pytest.raises(ValueError):
            foo.func("hello", *self._test_params)

        with pytest.raises(ValueError):
            params = list(self._test_params)
            params[0] = "hello"
            foo.func(5, *params)

    def test_func_err(self):
        foo = self.ff_class(
            params=self._test_params, param_errors=self._test_param_errors
        )

        for x in (0, 1.0, np.linspace(10, 30, num=20)):
            assert np.allclose(
                foo.func_err(x),
                self.func_err(x, self._test_params, self._test_param_errors),
            )

        x = [4, 5, 6]
        results = foo.func_err(x, x_err=0.1, rety=True)
        assert np.allclose(
            results[0],
            self.func_err(
                np.array(x), self._test_params, self._test_param_errors, x_err=0.1
            ),
        )
        assert np.allclose(results[1], self.func(np.array(x), *self._test_params))

        with pytest.raises(ValueError):
            foo.func_err("hello")

        with pytest.raises(ValueError):
            foo.func_err(5, x_err="goodbye")

        with pytest.raises(ValueError):
            foo.func_err(5, x_err=[0.1, 0.1])

    def test_call(self):
        foo = self.ff_class()
        foo.params = self._test_params
        foo.param_errors = self._test_param_errors

        for x in (0, 1.0, np.linspace(10, 30, num=20)):
            assert np.allclose(foo(x), self.func(x, *self._test_params))

            # also return error
            y, y_err = foo(x, reterr=True)
            assert np.allclose(y, self.func(x, *self._test_params))
            assert np.allclose(
                y_err, self.func_err(x, self._test_params, self._test_param_errors),
            )

        x = [4, 5, 6]
        x_err = 0.05
        assert np.allclose(foo(x), self.func(np.array(x), *self._test_params))
        y, y_err = foo(x, x_err=x_err, reterr=True)
        assert np.allclose(y, self.func(np.array(x), *self._test_params))
        assert np.allclose(
            y_err,
            self.func_err(
                np.array(x), self._test_params, self._test_param_errors, x_err=x_err
            ),
        )

        with pytest.raises(ValueError):
            foo("hello")

        with pytest.raises(ValueError):
            foo(5, x_err=[1, 2], reterr=True)

    @abstractmethod
    def test_root_solve(self):
        raise NotImplementedError

    def test_curve_fit(self):
        foo = self.ff_class()

        xdata = np.linspace(-10, 10)
        ydata = self.func(xdata, *self._test_params)

        assert foo.params is None
        assert foo.param_errors is None
        assert foo.rsq is None
        assert foo.curve_fit_results is None

        foo.curve_fit(xdata, ydata)

        assert foo.curve_fit_results is not None
        assert np.isclose(foo.rsq, 1.0)
        assert np.allclose(
            foo.param_errors, tuple([0] * len(foo.param_names)), atol=1.5e-8,
        )
        assert np.allclose(foo.params, self._test_params)


class TestFFExponential(BaseFFTests):
    """
    Tests for fit function class `plasmapy.analysis.fit_functions.Exponential`.
    """

    ff_class = ffuncs.Exponential
    _test_params = (5.0, 1.0)
    _test_param_errors = (0.1, 0.1)

    @staticmethod
    def func(x, a, alpha):
        return a * np.exp(alpha * x)

    def func_err(self, x, params, param_errors, x_err=None):
        a, alpha = params
        a_err, alpha_err = param_errors
        y = self.func(x, *params)

        a_term = (a_err / a) ** 2
        alpha_term = (x * alpha_err) ** 2

        err = a_term + alpha_term

        if x_err is not None:
            x_term = (alpha * x_err) ** 2
            err += x_term

        err = np.abs(y) * np.sqrt(err)

        return err

    def test_basics(self):
        super().test_basics()

        foo = self.ff_class()

        assert foo.param_names == ("a", "alpha")
        assert foo.latex_str == fr"A \, \exp(\alpha x)"
        assert foo.__str__() == f"f(x) = A exp(alpha x)"

    def test_root_solve(self):
        foo = self.ff_class(params=(1, 1), param_errors=(0, 0))
        root, err = foo.root_solve()
        assert np.isnan(root)
        assert np.isnan(err)


class TestFFExponentialPlusLinear(BaseFFTests):
    """
    Tests for fit function class
    `plasmapy.analysis.fit_functions.ExponentialPlusLinear`.
    """

    ff_class = ffuncs.ExponentialPlusLinear
    _test_params = (2.0, 1.0, 5.0, -10.0)
    _test_param_errors = (0.1, 0.1, 0.1, 0.1)

    @staticmethod
    def func(x, a, alpha, m, b):
        return a * np.exp(alpha * x) + m * x + b

    def func_err(self, x, params, param_errors, x_err=None):
        a, alpha, m, b = params
        a_err, alpha_err, m_err, b_err = param_errors

        exp_y = a * np.exp(alpha * x)

        a_term = (exp_y * a_err / a) ** 2
        alpha_term = (exp_y * x * alpha_err) ** 2
        m_term = (m_err * x) ** 2
        b_term = b_err ** 2

        err = a_term + alpha_term + m_term + b_term

        if x_err is not None:
            x_term = (exp_y * alpha * x_err) ** 2
            x_term += (m * x_err) ** 2
            x_term += 2 * a * alpha * m * np.exp(alpha * x) * (x_err ** 2)

            err += x_term

        err = np.sqrt(err)

        return err

    def test_basics(self):
        super().test_basics()

        foo = self.ff_class()

        assert foo.param_names == ("a", "alpha", "m", "b")
        assert foo.latex_str == fr"A \, \exp(\alpha x) + m x + b"
        assert foo.__str__() == f"f(x) = A exp(alpha x) + m x + b"

    def test_root_solve(self):
        foo = self.ff_class(params=(5.0, 0.5, 1.0, 5.0), param_errors=(0, 0, 0, 0))
        root, err = foo.root_solve(-5)
        assert np.isclose(root, -5.345338)
        assert np.isnan(err)


class TestFFExponentialPlusOffset(BaseFFTests):
    """
    Tests for fit function class
    `plasmapy.analysis.fit_functions.ExponentialPlusOffset`.
    """

    ff_class = ffuncs.ExponentialPlusOffset
    _test_params = (2.0, 1.0, -10.0)
    _test_param_errors = (0.1, 0.1, 0.1)

    @staticmethod
    def func(x, a, alpha, b):
        return a * np.exp(alpha * x) + b

    def func_err(self, x, params, param_errors, x_err=None):
        a, alpha, b = params
        a_err, alpha_err, b_err = param_errors

        exp_y = a * np.exp(alpha * x)

        a_term = (exp_y * a_err / a) ** 2
        alpha_term = (exp_y * x * alpha_err) ** 2
        b_term = b_err ** 2

        err = a_term + alpha_term + b_term

        if x_err is not None:
            x_term = (exp_y * alpha * x_err) ** 2
            err += x_term

        err = np.sqrt(err)

        return err

    def test_basics(self):
        super().test_basics()

        foo = self.ff_class()

        assert foo.param_names == ("a", "alpha", "b")
        assert foo.latex_str == fr"A \, \exp(\alpha x) + b"
        assert foo.__str__() == f"f(x) = A exp(alpha x) + b"

    def test_root_solve(self):
        foo = self.ff_class(params=(3.0, 0.5, -5.0), param_errors=(0, 0, 0))
        root, err = foo.root_solve()
        assert root == np.log(5.0 / 3.0) / 0.5
        assert err == 0

        foo.params = (3.0, 0.5, 5.0)
        root, err = foo.root_solve()
        assert np.isnan(root)
        assert np.isnan(err)


class TestFFLinear(BaseFFTests):
    """
    Tests for fit function class `plasmapy.analysis.fit_functions.Linear`.
    """

    ff_class = ffuncs.Linear
    _test_params = (5.0, 4.0)
    _test_param_errors = (0.1, 0.1)

    @staticmethod
    def func(x, m, b):
        return m * x + b

    def func_err(self, x, params, param_errors, x_err=None):
        m, b = params
        m_err, b_err = param_errors

        m_term = (m_err * x) ** 2
        b_term = b_err ** 2
        err = m_term + b_term

        if x_err is not None:
            x_term = (m * x_err) ** 2
            err += x_term
        err = np.sqrt(err)

        return err

    def test_basics(self):
        super().test_basics()

        foo = self.ff_class()

        assert foo.param_names == ("m", "b")
        assert foo.latex_str == fr"m x + b"
        assert foo.__str__() == f"f(x) = m x + b"

    def test_root_solve(self):
        foo = self.ff_class(params=(1, 1), param_errors=(0, 0))
        assert foo.root_solve() == (-1, 0)

        foo.params = (5.0, 1.3)
        foo.param_errors = (0.1, 0.1)
        root, err = foo.root_solve()
        assert root == -1.3 / 5.0
        assert err == np.abs(root) * np.sqrt((0.1 / 5.0) ** 2 + (0.1 / 1.3) ** 2)
