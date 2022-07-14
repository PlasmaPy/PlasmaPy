"""
Tests for the fitting function classes defined in `plasmapy.analysis.fit_functions`.
"""
import numpy as np
import pytest

from abc import ABC, abstractmethod
from contextlib import nullcontext as does_not_raise

import plasmapy.analysis.fit_functions as ffuncs


class TestAbstractFitFunction:
    """
    Tests for fit function class `plasmapy.analysis.fit_functions.AbstractFitFunction`.

    Notes
    -----
    Since `AbstractFitFunction` can not be instantiated, its complete
    functionality can not be directly tested.  To resolve this, `BaseFFTests`,
    which is the base test class for the fit function classes, will test all
    the functionality within `AbstractFitFunction`.
    """

    ff_class = ffuncs.AbstractFitFunction

    def test_is_abc(self):
        """Test `AbstractFitFunction` is an abstract base class."""
        assert issubclass(self.ff_class, ABC)

    @pytest.mark.parametrize(
        "name, isproperty",
        [
            ("__call__", False),
            ("curve_fit", False),
            ("curve_fit_results", True),
            ("func", False),
            ("func_err", False),
            ("latex_str", True),
            ("param_errors", True),
            ("param_names", True),
            ("params", True),
            ("rsq", True),
            ("root_solve", False),
        ],
    )
    def test_methods(self, name, isproperty):
        """Test for required methods and properties."""
        assert hasattr(self.ff_class, name)

        if isproperty:
            assert isinstance(getattr(self.ff_class, name), property)

    @pytest.mark.parametrize(
        "name",
        ["__str__", "func", "func_err", "latex_str"],
    )
    def test_abstractmethods(self, name):
        """Test for required abstract methods."""
        assert name in self.ff_class.__abstractmethods__


class BaseFFTests(ABC):
    abc = ffuncs.AbstractFitFunction
    _test_params = NotImplemented  # type: tuple
    _test_param_errors = NotImplemented  # type: tuple
    _test_param_names = NotImplemented  # type: tuple
    _test_latex_str = NotImplemented  # type: str
    _test__str__ = NotImplemented  # type: str

    @property
    @abstractmethod
    def ff_class(self):
        """Fit function class to be tested."""
        ...

    @staticmethod
    @abstractmethod
    def func(x, *args):
        """
        Formula/Function that the fit function class is suppose to be modeling.
        This is used to test the fit function `func` method.
        """
        ...

    @abstractmethod
    def func_err(self, x, params, param_errors, x_err=None):
        """
        Function representing the propagation of error associated with fit function
        model. This is used to test the fit function `func_err` method.
        """
        ...

    def test_inheritance(self):
        """Test inheritance from `AbstractFitFunction`."""
        assert issubclass(self.ff_class, self.abc)

    def test_iscallable(self):
        """Test instantiated fit function is callable."""
        assert callable(self.ff_class())

    def test_repr(self):
        """Test __repr__."""
        ff_obj = self.ff_class()
        assert ff_obj.__repr__() == f"{ff_obj.__str__()} {ff_obj.__class__}"

    @pytest.mark.parametrize(
        "name, isproperty",
        [
            ("__call__", False),
            ("_param_names", False),
            ("curve_fit", False),
            ("curve_fit_results", True),
            ("func", False),
            ("func_err", False),
            ("latex_str", True),
            ("param_errors", True),
            ("param_names", True),
            ("params", True),
            ("rsq", True),
            ("root_solve", False),
        ],
    )
    def test_methods(self, name, isproperty):
        """Test attribute/method/property existence."""
        assert hasattr(self.ff_class, name)

        if isproperty:
            assert isinstance(getattr(self.ff_class, name), property)

        if name == "_param_names" and self.ff_class._param_names == NotImplemented:
            pytest.fail(
                f"{self.ff_class} class attribute '_param_names' needs to "
                f" be defined as a tuple of strings representing the names of "
                f"the fit parameters."
            )

    @pytest.mark.parametrize(
        "name, value_ref_name",
        [
            ("param_names", "_test_param_names"),
            ("latex_str", "_test_latex_str"),
            ("__str__", "_test__str__"),
        ],
    )
    def test_abstractmethod_values(self, name, value_ref_name):
        """Test value of all abstract methods, except `func` and `func_err`."""
        ff_obj = self.ff_class()

        value = getattr(ff_obj, name)
        if callable(value):
            value = value()

        exp_value = getattr(self, value_ref_name)
        if exp_value == NotImplemented:
            pytest.fail(
                f"The expected value for abstract method {name} is not "
                f"implemented/defined in the test class attribute {value_ref_name}."
            )

        assert value == exp_value

    @pytest.mark.parametrize(
        "params, param_errors, with_condition",
        [
            (None, None, does_not_raise()),
            ("default", "default", does_not_raise()),
            (5, None, pytest.raises(ValueError)),
            (None, 5, pytest.raises(ValueError)),
            (["wrong"], None, pytest.raises(ValueError)),
            (None, ["wrong"], pytest.raises(ValueError)),
            ("default+", None, pytest.raises(ValueError)),
            (None, "default+", pytest.raises(ValueError)),
        ],
    )
    def test_instantiation(self, params, param_errors, with_condition):
        """Test behavior of fit function class instantiation."""
        if params == "default":
            params = self._test_params
        elif params == "default+":
            params = self._test_params
            params = list(params)
            params.append(5)

        if param_errors == "default":
            param_errors = self._test_param_errors
        elif param_errors == "default+":
            param_errors = self._test_param_errors
            param_errors = list(param_errors)
            param_errors.append(5)

        with with_condition:
            ff_obj = self.ff_class(params=params, param_errors=param_errors)

            assert ff_obj.curve_fit_results is None
            assert ff_obj.rsq is None

            if params is None:
                assert ff_obj.params is None
            else:
                assert ff_obj.params == ff_obj.FitParamTuple(*params)

            if param_errors is None:
                assert ff_obj.param_errors is None
            else:
                assert ff_obj.param_errors == ff_obj.FitParamTuple(*param_errors)

    def test_param_namedtuple(self):
        """
        Test that the namedtuple used for `params` and `param_errors` is
        constructed correctly.
        """
        ff_obj = self.ff_class()
        assert hasattr(ff_obj, "FitParamTuple")
        assert issubclass(ff_obj.FitParamTuple, tuple)
        for name in ff_obj.param_names:
            assert hasattr(ff_obj.FitParamTuple, name)

    def test_param_names(self):
        """Test attribute `param_names` is defined correctly."""
        ff_obj = self.ff_class()
        assert isinstance(ff_obj.param_names, tuple)
        assert len(ff_obj.param_names) != 0
        assert all(isinstance(val, str) for val in ff_obj.param_names)

    @pytest.mark.parametrize(
        "params, extra, with_condition",
        [
            ([2], None, does_not_raise()),
            (5, None, pytest.raises(ValueError)),
            (["wrong"], None, pytest.raises(ValueError)),
            ([3], 10, pytest.raises(ValueError)),
        ],
    )
    def test_params_setting(self, params, extra, with_condition):
        """Tests for property setting of attribute `params`."""
        ff_obj = self.ff_class()

        if isinstance(params, list) and len(params) == 1:
            params = params * len(ff_obj.param_names)
        if extra is not None:
            params.append(extra)

        with with_condition:
            ff_obj.params = params
            assert ff_obj.params == ff_obj.FitParamTuple(*params)

    @pytest.mark.parametrize(
        "param_errors, extra, with_condition",
        [
            ([2], None, does_not_raise()),
            (5, None, pytest.raises(ValueError)),
            (["wrong"], None, pytest.raises(ValueError)),
            ([3], 10, pytest.raises(ValueError)),
        ],
    )
    def test_param_errors_setting(self, param_errors, extra, with_condition):
        """Tests for property setting of attribute `param_errors`."""
        ff_obj = self.ff_class()

        if isinstance(param_errors, list) and len(param_errors) == 1:
            param_errors = param_errors * len(ff_obj.param_names)
        if extra is not None:
            param_errors.append(extra)

        with with_condition:
            ff_obj.param_errors = param_errors
            assert ff_obj.param_errors == ff_obj.FitParamTuple(*param_errors)

    @pytest.mark.parametrize(
        "x, replace_a_param, with_condition",
        [
            (0, None, does_not_raise()),
            (1.0, None, does_not_raise()),
            (np.linspace(10, 30, num=20), None, does_not_raise()),
            ([4, 5, 6], None, does_not_raise()),
            ("hello", None, pytest.raises(TypeError)),
            (5, "hello", pytest.raises(TypeError)),
        ],
    )
    def test_func(self, x, replace_a_param, with_condition):
        """Test the `func` method."""
        ff_obj = self.ff_class()

        params = self._test_params
        if replace_a_param is not None:
            params = list(params)
            params[0] = replace_a_param

        with with_condition:
            y = ff_obj.func(x, *params)

            if isinstance(x, list):
                x = np.array(x)
            y_expected = self.func(x, *params)

            assert np.allclose(y, y_expected)

    @pytest.mark.parametrize(
        "x, kwargs, with_condition",
        [
            (0, {}, does_not_raise()),
            (1.0, {}, does_not_raise()),
            (np.linspace(10, 30, num=20), {}, does_not_raise()),
            ([4, 5, 6], {"x_err": 0.1, "rety": True}, does_not_raise()),
            ("hello", {}, pytest.raises(TypeError)),
            (5, {"x_err": "goodbye"}, pytest.raises(TypeError)),
            (5, {"x_err": [0.1, 0.1]}, pytest.raises(ValueError)),
        ],
    )
    def test_func_err(self, x, kwargs, with_condition):
        """Test the `func_err` method."""
        params = self._test_params
        param_errors = self._test_param_errors
        ff_obj = self.ff_class(params=params, param_errors=param_errors)

        with with_condition:
            results = ff_obj.func_err(x, **kwargs)
            if "rety" in kwargs and kwargs["rety"]:
                y_err, y = results
            else:
                y_err = results
                y = None

            x_err = kwargs["x_err"] if "x_err" in kwargs else None
            if isinstance(x, list):
                x = np.array(x)
            y_err_expected = self.func_err(x, params, param_errors, x_err=x_err)

            assert np.allclose(y_err, y_err_expected)

            if y is not None:
                assert np.allclose(y, self.func(x, *params))

    @pytest.mark.parametrize(
        "x, kwargs, with_condition",
        [
            (0, {}, does_not_raise()),
            (0, {"reterr": True}, does_not_raise()),
            (1.0, {}, does_not_raise()),
            (1.0, {"reterr": True}, does_not_raise()),
            (np.linspace(10, 30, num=20), {}, does_not_raise()),
            (np.linspace(10, 30, num=20), {"reterr": True}, does_not_raise()),
            ([4, 5, 6], {}, does_not_raise()),
            ([4, 5, 6], {"x_err": 0.05, "reterr": True}, does_not_raise()),
            ("hello", {}, pytest.raises(TypeError)),
            (5, {"x_err": [1, 2], "reterr": True}, pytest.raises(ValueError)),
        ],
    )
    def test_call(self, x, kwargs, with_condition):
        """Test __call__ behavior."""
        params = self._test_params
        param_errors = self._test_param_errors
        ff_obj = self.ff_class(params=params, param_errors=param_errors)

        reterr = kwargs["reterr"] if "reterr" in kwargs else False
        x_err = kwargs["x_err"] if "x_err" in kwargs else None
        with with_condition:
            results = ff_obj(x, **kwargs)
            if reterr:
                y = results[0]
                y_err = results[1]
            else:
                y = results

            if isinstance(x, list):
                x = np.array(x)
            y_expected = self.func(x, *params)

            assert np.allclose(y, y_expected)

            if reterr:
                y_err_expected = self.func_err(x, params, param_errors, x_err=x_err)
                assert np.allclose(y_err, y_err_expected)

    @abstractmethod
    def test_root_solve(self):
        ...

    def test_curve_fit(self):
        """Test the `curve_fit` method."""
        ff_obj = self.ff_class()

        xdata = np.linspace(-10, 10)
        ydata = self.func(xdata, *self._test_params)

        assert ff_obj.params is None
        assert ff_obj.param_errors is None
        assert ff_obj.rsq is None
        assert ff_obj.curve_fit_results is None

        ff_obj.curve_fit(xdata, ydata)

        assert ff_obj.curve_fit_results is not None
        assert np.isclose(ff_obj.rsq, 1.0)
        assert np.allclose(
            ff_obj.param_errors,
            tuple([0] * len(ff_obj.param_names)),
            atol=1.5e-8,
        )
        assert np.allclose(ff_obj.params, self._test_params)


class TestFFExponential(BaseFFTests):
    """
    Tests for fit function class `plasmapy.analysis.fit_functions.Exponential`.
    """

    ff_class = ffuncs.Exponential
    _test_params = (5.0, 1.0)
    _test_param_errors = (0.1, 0.1)
    _test_param_names = ("a", "alpha")
    _test_latex_str = r"a \, \exp(\alpha x)"
    _test__str__ = "f(x) = a exp(alpha x)"

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

    def test_root_solve(self):
        ff_obj = self.ff_class(params=(1, 1), param_errors=(0, 0))
        root, err = ff_obj.root_solve()
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
    _test_param_names = ("a", "alpha", "m", "b")
    _test_latex_str = r"a \, \exp(\alpha x) + m x + b"
    _test__str__ = "f(x) = a exp(alpha x) + m x + b"

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
        b_term = b_err**2

        err = a_term + alpha_term + m_term + b_term

        if x_err is not None:
            x_term = (exp_y * alpha * x_err) ** 2
            x_term += (m * x_err) ** 2
            x_term += 2 * a * alpha * m * np.exp(alpha * x) * (x_err**2)

            err += x_term

        err = np.sqrt(err)

        return err

    def test_root_solve(self):
        ff_obj = self.ff_class(params=(5.0, 0.5, 1.0, 5.0), param_errors=(0, 0, 0, 0))
        root, err = ff_obj.root_solve(-5)
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
    _test_param_names = ("a", "alpha", "b")
    _test_latex_str = r"a \, \exp(\alpha x) + b"
    _test__str__ = "f(x) = a exp(alpha x) + b"

    @staticmethod
    def func(x, a, alpha, b):
        return a * np.exp(alpha * x) + b

    def func_err(self, x, params, param_errors, x_err=None):
        a, alpha, b = params
        a_err, alpha_err, b_err = param_errors

        exp_y = a * np.exp(alpha * x)

        a_term = (exp_y * a_err / a) ** 2
        alpha_term = (exp_y * x * alpha_err) ** 2
        b_term = b_err**2

        err = a_term + alpha_term + b_term

        if x_err is not None:
            x_term = (exp_y * alpha * x_err) ** 2
            err += x_term

        err = np.sqrt(err)

        return err

    def test_root_solve(self):
        ff_obj = self.ff_class(params=(3.0, 0.5, -5.0), param_errors=(0, 0, 0))
        root, err = ff_obj.root_solve()
        assert root == np.log(5.0 / 3.0) / 0.5
        assert err == 0

        ff_obj.params = (3.0, 0.5, 5.0)
        with pytest.warns(RuntimeWarning, match="invalid value encountered in log"):
            root, err = ff_obj.root_solve()
            assert np.isnan(root)
            assert np.isnan(err)


class TestFFLinear(BaseFFTests):
    """
    Tests for fit function class `plasmapy.analysis.fit_functions.Linear`.
    """

    ff_class = ffuncs.Linear
    _test_params = (5.0, 4.0)
    _test_param_errors = (0.1, 0.1)
    _test_param_names = ("m", "b")
    _test_latex_str = r"m x + b"
    _test__str__ = "f(x) = m x + b"

    @staticmethod
    def func(x, m, b):
        return m * x + b

    def func_err(self, x, params, param_errors, x_err=None):
        m, b = params
        m_err, b_err = param_errors

        m_term = (m_err * x) ** 2
        b_term = b_err**2
        err = m_term + b_term

        if x_err is not None:
            x_term = (m * x_err) ** 2
            err += x_term
        err = np.sqrt(err)

        return err

    @pytest.mark.parametrize(
        "params, param_errors, root, root_err, conditional",
        [
            ((1, 1), (0, 0), -1, 0, does_not_raise()),
            (
                (5.0, 1.3),
                (0.1, 0.1),
                -1.3 / 5.0,
                np.abs(-1.3 / 5.0) * np.sqrt((0.1 / 5.0) ** 2 + (0.1 / 1.3) ** 2),
                does_not_raise(),
            ),
            ((0.3, 0.0), (0.1, 0.1), 0.0, np.abs(0.1 / 0.3), does_not_raise()),
            ((0.0, 1.0), (0.1, 0.1), np.nan, np.nan, pytest.warns(RuntimeWarning)),
        ],
    )
    def test_root_solve(self, params, param_errors, root, root_err, conditional):
        with conditional:
            ff_obj = self.ff_class(params=params, param_errors=param_errors)
            results = ff_obj.root_solve()

            if np.all(np.isnan([root, root_err])):
                assert np.all(np.isnan(results))
            elif np.isnan(root):
                assert np.isnan(results[0])
                assert np.isclose(results[1], root_err)
            elif np.isnan(root_err):
                assert np.isclose(results[0], root)
                assert np.isnan(results[1])
            else:
                assert np.allclose(results, [root, root_err])
