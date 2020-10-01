"""
`FitFunction` classes designed to assist in curve fitting of swept Langmuir
traces.
"""
__all__ = [
    "AbstractFitFunction",
    "Exponential",
    "ExponentialPlusLinear",
    "ExponentialPlusOffset",
    "Linear",
]

import numpy as np

from abc import ABC, abstractmethod
from collections import namedtuple
from scipy.stats import linregress
from scipy.optimize import curve_fit, fsolve
from typing import Any, Tuple, Union


class AbstractFitFunction(ABC):
    """
    Abstract class for defining fit functions :math:`f(x)` and the tools for
    fitting the function to a set of data.  These were originally designed for
    assisting in fitting curves to swept Langmuir data.
    """

    _parameter_names = NotImplemented  # type: Tuple[str, ...]


    def __init__(self):
        self.FitParamTuple = namedtuple("FitParamTuple", self._parameter_names)
        """
        A named tuple class used for attributes :attr:`parameters` and 
        :attr:`parameters_err`.  The attribute :attr:`parameter_names` defines
        the tuple field names.
        """

        self._parameters = None  # type: Union[None, Tuple[Any, ...]]
        self._parameters_err = None  # type: Union[None, Tuple[Any, ...]]
        self._covariance_matrix = None
        self._rsq = None
        self._curve_fit_results = None

    def __call__(self, x, x_err=None, reterr=False):
        """
        Direct call of the fit function :math:`f(x)``.

        Parameters
        ----------
        x: array_like
            Dependent variables.

        x_err: array_like, optional
            Errors associated with the independent variables `x`.  Must be of
            size one or equal to the size of `x`.

        reterr: bool, optional
            (Default: `False`) If `True`, return an array of uncertainties
            associated with the calculated independent variables

        Returns
        -------
        y: `numpy.ndarray`
            Corresponding dependent variables :math:`y=f(x)` of the independent
            variables :math:`x`.

        y_err: `numpy.ndarray`
            Uncertainties associated with the calculated dependent variables
            :math:`\\delta y`
        """
        if not isinstance(x, np.ndarray):
            x = np.array(x)

        y = self._func(x, *self.parameters)

        if reterr:
            try:
                y_err = self._func_err(x, y, x_err=x_err)
            except NotImplementedError:
                y_err = np.tile(np.nan, x.shape)

            return y, y_err

        return y

    def __repr__(self):
        return f"{self.__str__()} {self.__class__}"

    def __str__(self):
        return f"Unspecified f(x)"

    @staticmethod
    @abstractmethod
    def func(x, *args):
        """
        The fit function.  This signature of the function must first take the
        independent variable followed by the parameters to be fitted as
        separate arguments.

        When sub-classing the definition should look something like::

            def func(x, a, b, c):
                return a * x ** 2 + b * x + c

        Parameters
        ----------
        x: array_like
            Independent variables to be passed to the fit function.

        *args
            The parameters that will be adjusted to make the fit.

        Returns
        -------
        `numpy.ndarray`:
            The calculated dependent variables of the independent variables `x`.
        """
        raise NotImplementedError

    @abstractmethod
    def _func_err(self, x, y, x_err=None):
        """
        Calculate dependent variable uncertainties :math:`\\delta y` for
        dependent variables :math:`y=f(x)`.

        Parameters
        ----------
        x: array_like
            Independent variables to be passed to the fit function.

        y: array_like
            Dependent variables associated with :math:`x`, :math:`f(x)`.

        x_err: array_like, optional
            Errors associated with the independent variables `x`.  Must be of
            size one or equal to the size of `x`.

        Returns
        -------
        `numpy.ndarray`:
            The calculated uncertainties of the dependent variables of the
            independent variables `x`.
        """
        raise NotImplementedError

    @property
    def curve_fit_results(self):
        """
        The results returned by the curve fitting routine used by
        :attr:`curve_fit`.  This is typically from `scipy.stats.linregress` or
        `scipy.optimize.curve_fit`.
        """
        return self._curve_fit_results

    @property
    def parameters(self) -> Union[None, tuple]:
        """The fitted parameters for the fit function."""
        if self._parameters is None:
            return self._parameters
        else:
            return self.FitParamTuple(*self._parameters)

    @parameters.setter
    def parameters(self, val) -> None:
        if isinstance(val, self.FitParamTuple):
            self._parameters = tuple(val)
        elif isinstance(val, (tuple, list)) and len(val) == len(self.parameter_names):
            self._parameters = tuple(val)
        else:
            raise ValueError(f"Got type {type(val)} for 'val', expecting tuple of "
                             f"length {len(self.parameter_names)}.")

    @property
    def parameters_err(self) -> Union[None, tuple]:
        """The associated errors of the fit `parameters`."""
        if self._parameters_err is None:
            return self._parameters_err
        else:
            return self.FitParamTuple(*self._parameters_err)

    @parameters_err.setter
    def parameters_err(self, val) -> None:
        if isinstance(val, self.FitParamTuple):
            self._parameters_err = tuple(val)
        elif isinstance(val, (tuple, list)) and len(val) == len(self.parameter_names):
            self._parameters_err = tuple(val)
        else:
            raise ValueError(f"Got type {type(val)} for 'val', expecting tuple of "
                             f"length {len(self.parameter_names)}.")

    @property
    def parameter_names(self) -> Tuple[str, ...]:
        """Names of the fitted parameters."""
        return self._parameter_names

    @property
    @abstractmethod
    def latex_str(self) -> str:
        """Latex friendly representation of the fit function."""
        raise NotImplementedError

    def root_solve(self, x0, **kwargs):
        """
        Solve for the root of the fit function (i.e. :math:`f(x_r) = 0`).  This
        mehtod used `scipy.optimize.fsolve` to find the function roots.

        Parameters
        ----------
        x0: `~numpy.ndarray`
            The starting estimate for the roots of :math:`f(x_r) = 0`.

        **kwargs
            Any keyword accepted by `scipy.optimize.fsolve`, except for `args`.

        Returns
        -------
        x : `~numpy.ndarray`
            The solution (or the result of the last iteration for an
            unsuccessful call).

        x_err: `~numpy.ndarray`
            The uncertainty associated with the root calculation.  **Currently
            this returns an array of** `numpy.nan` **values equal in shape to**
            `x` **, since there is no determined way to calculate the
            uncertaintyes.**

        Notes
        -----
        If the full output of `scipy.optimize.fsolve` is desired then one can do

            >>> import numpy as np
            >>> import scipy

            >>> class SomeFunc(AbstractFitFunction):
            ...     _parameter_names = ("m", "b")
            ...
            ...     @property
            ...     def latex_str(self) -> str:
            ...         return f"m \\, x + b"
            ...
            ...     @staticmethod
            ...     def func(x, m, b):
            ...         return m * x + b
            ...
            ...     def _func_err(self, x, y, x_err=None):
            ...         m, b = self.parameters
            ...         m_err, b_err = self.parameters_err
            ...
            ...         m_term = x * m_err
            ...         b_term = b_err
            ...         err = m_term ** 2 + b_term ** 2
            ...
            ...         if x_err is not None:
            ...             x_term = m * x_err
            ...             err == x_term ** 2
            ...
            ...         return np.sqrt(err)
            ...
            >>> func = SomeFunc()
            >>> func.parameters = (1., 5.)
            >>> func.parameters_err = (0.0, 0.0)
            >>> roots = scipy.optimize.fsolve(func, -4., full_output=True)
            >>> roots
            (array([-5.]),
             {'nfev': 4,
              'fjac': array([[-1.]]),
              'r': array([-1.]),
              'qtf': array([2.18...e-12]),
              'fvec': array([0.])},
             1,
             'The solution converged.')

        """
        kwargs["args"] = self.parameters
        results = fsolve(self.func, x0, **kwargs)
        if isinstance(results, tuple):
            results = results[0]

        return results, np.tile(np.nan, results.shape)

    @property
    def rsq(self):
        """
        Coefficient of determination (r-squared) value of the fit.

        .. math::

            r^2 &= 1 - \\frac{SS_{res}}{SS_{tot}}

            SS_{res} &= \\sum\\limits_{i} (y_i - f(x_i))^2

            SS_{tot} &= \\sum\\limits_{i} (y_i - \\bar{y})^2

        where :math:`(x_i, y_i)` are the sample data pairs, :math:`f(x_i)` is
        the fitted dependent variable corresponding to :math:`x_i`, and
        :math:`\\bar{y}` is the average of the :math:`y_i` values.

        """
        return self._rsq

    def curve_fit(self, xdata, ydata, **kwargs) -> None:
        """
        Use a non-linear least squares method to fit the fit function to
        (`xdata`, `ydata`), using `scipy.optimize.curve_fit`.  This will set
        the attributes :attr:`parameters`, :attr:`parameters_err`, and
        :attr:`rsq`.

        The results of `scipy.optimize.curve_fit` can be obtained via
        :attr:`curve_fit_results`.

        Parameters
        ----------
        xdata: array_like
            The independent variable where data is measured.  Should be 1D of
            length M.

        ydata: array_like
            The dependent data associated with `xdata`.

        **kwargs
            Any keywords accepted by `scipy.optimize.curve_fit`.

        Raises
        ------
        ValueError
            if either `ydata` or `xdata` contain NaNs, or if incompatible options
            are used.

        RuntimeError
            if the least-squares minimization fails.

        ~scipy.optimize.OptimizeWarning
            if covariance of the parameters can not be estimated.

        """
        popt, pcov = curve_fit(self.func, xdata, ydata, **kwargs)
        self._curve_fit_results = (popt, pcov)
        self.parameters = tuple(popt.tolist())
        self.parameters_err = tuple(np.sqrt(np.diag(pcov)).tolist())

        # calc rsq
        # rsq = 1 - (ss_res / ss_tot)
        residuals = ydata - self.func(xdata, *self.params)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
        self._rsq = 1 - (ss_res / ss_tot)


class Exponential(AbstractFitFunction):
    """
    A sub-class of `AbstractFitFunction` to represent an exponential with an
    offset.

    .. math::

        y &= f(x) = A \\, e^{\\alpha \\, x}

        \\left( \\frac{\\delta y}{|y|} \\right)^2 &=
            \\left( \\frac{\\delta A}{A} \\right)^2
            + (x \\, \\delta \\alpha)^2
            + (\\alpha \\, \\delta x)^2

    where :math:`A` and :math:`\\alpha` are the real constants to be fitted and
    :math:`x` is the independent variable.  :math:`\\delta A`,
    :math:`\\delta \\alpha`, and :math:`\\delta x` are the respective
    uncertainties for :math:`A`, :math:`\\alpha`, and :math:`x`.
    """
    _parameter_names = ("a", "alpha")

    def __str__(self):
        return f"f(x) = A exp(alpha x)"

    @staticmethod
    def func(x, a, alpha):
        return a * np.exp(alpha * x)

    def _func_err(self, x, y, x_err=None):
        a, alpha = self.parameters
        a_err, alpha_err = self.parameters_err

        a_term = (a_err / a) ** 2
        alpha_term = (x * alpha_err) ** 2

        err = a_term + alpha_term

        if x_err is not None:
            x_term = (alpha * x_err) ** 2
            err += x_term

        err = np.abs(y) * np.sqrt(err)

        return err

    @property
    def latex_str(self) -> str:
        return fr"A \, \exp(\alpha \, x)"

    def root_solve(self, *args, **kwargs):
        """
        The root :math:`f(x_r) = 0` for the fit function. **An exponential has no
        real roots.**

        Parameters
        ----------
        *args
            Not needed.  This is to ensure signature comparability with
            `AbstractFitFunction`.

        *kwargs
            Not needed.  This is to ensure signature comparability with
            `AbstractFitFunction`.

        Returns
        -------
        root: float
            The root value for the given fit :attr:`parameters`.

        err: float
            The uncertainty in the calculated root for the given fit
            :attr:`parameters` and :attr:`parameters_err`.
        """

        return np.nan, np.nan


class Linear(AbstractFitFunction):
    """
    A sub-class of `AbstractFitFunction` to represent a linear function.

    .. math::

        y &= f(x) = m \\, x + b

        (\\delta y)^2 &= (x \\, \\delta m)^2 + (m \\, \\delta x)^2 + (\\delta b)^2

    where :math:`m` and :math:`b` are real constants to be fitted and :math:`x` is
    the independent variable.  :math:`\\delta m`, :math:`\\delta b`, and
    :math:`\\delta x` are the respective uncertainties for :math:`m`, :math:`b`,
    and :math:`x`.
    """

    _parameter_names = ("m", "b")

    def __str__(self):
        return f"f(x) = m x + b"

    @staticmethod
    def func(x, m, b):
        """
        The fit function, a linear function.

        .. math::

            f(x) = m \\, x + b

        where :math:`m` and :math:`b` are positive real constants representing the
        slope and intercept, respectively, and :math:`x` is the independent
        variable.

        Parameters
        ----------
        x: array_like
            Independent variable.
        m: float
            value for slope :math:`m`

        b: float
            value for intercept :math:`b`

        Returns
        -------
        y: array_like
            dependent variables corresponding to :math:`x`

        """
        return m * x + b

    def _func_err(self, x, y, x_err=None):
        """
        Calculate dependent variable uncertainties :math:`\\delta y` for
        dependent variables :math:`y=f(x)`.

        .. math::

            (\\delta y)^2 &= (x \\, \\delta m)^2 + (m \\, \\delta x)^2 + (\\delta b)^2

        Parameters
        ----------
        x: array_like
            Independent variables to be passed to the fit function.

        Returns
        -------
        `numpy.ndarray`:
            The calculated uncertainty of the dependent variables of the
            independent variables `x`.
        """
        m, b = self.parameters
        m_err, b_err = self.parameters_err

        m_term = (m_err * x) ** 2
        b_term = b_err ** 2
        err = m_term + b_term

        if x_err is not None:
            x_term = (m * x_err) ** 2
            err += x_term

        return np.sqrt(err)

    @property
    def latex_str(self) -> str:
        return fr"m \, x + b"

    @property
    def rsq(self):
        """
        Coefficient of determination (r-squared) value of the fit.  Calculated
        by `scipy.stats.linregress` from the fit.
        """
        return self._rsq

    def root_solve(self, *args, **kwargs):
        """
        The root :math:`f(x_r) = 0` for the fit function.

        .. math::

            x_r &= \\frac{-b}{m}

            \\delta x_r &= |x_r| \\sqrt{
                \\left( \\frac{\\delta m}{m} \\right)^2
                + \\left( \\frac{\\delta b}{b} \\right)^2
            }

        Parameters
        ----------
        *args
            Not needed.  This is to ensure signature comparability with
            `AbstractFitFunction`.

        *kwargs
            Not needed.  This is to ensure signature comparability with
            `AbstractFitFunction`.

        Returns
        -------
        root: float
            The root value for the given fit :attr:`parameters`.

        err: float
            The uncertainty in the calculated root for the given fit
            :attr:`parameters` and :attr:`parameters_err`.
        """
        m, b = self.parameters
        root = -b / m

        m_err, b_err = self.parameters_err
        err = np.abs(root) * np.sqrt((m_err / m) ** 2 + (b_err / b) ** 2)

        return root, err

    def curve_fit(self, xdata, ydata, **kwargs) -> None:
        """
        Calculate a linear least-squares regression of (`xdata`, `ydata`) using
        `scipy.stats.linregress`.  This will set the attributes
        :attr:`parameters`, :attr:`parameters_err`, and :attr:`rsq`.

        The results of `scipy.stats.linregress` can be obtained via
        :attr:`curve_fit_results`.

        Parameters
        ----------
        xdata: array_like
            The independent variable where data is measured.  Should be 1D of
            length M.

        ydata: array_like
            The dependent data associated with `xdata`.

        **kwargs
            Any keywords accepted by `scipy.stats.linregress.curve_fit`.

        """
        results = linregress(xdata, ydata)
        self._curve_fit_results = results

        m = results[0]
        b = results[1]
        self.parameters = (m, b)

        m_err = results[4]
        b_err = np.sum(xdata ** 2) - ((np.sum(xdata) ** 2) / xdata.size)
        b_err = m_err * np.sqrt(1.0 / b_err)
        self.parameters_err = (m_err, b_err)

        self._rsq = results[2] ** 2


class ExponentialPlusLinear(AbstractFitFunction):
    """
    A sub-class of `AbstractFitFunction` to represent an exponential with an
    linear offset.

    .. math::

        y =& f(x) = A \\, e^{\\alpha \\, x} + m \\, x + b\\\\
        (\\delta y)^2 =&
            \\left( A e^{\\alpha x}\\right)^2 \\left[
                \\left( \\frac{\\delta A}{A} \\right)^2
                + (x \\, \\delta \\alpha)^2
                + (\\alpha \\, \\delta x)^2
            \\right]\\\\
            & + \\left(2 \\, A \\, \\alpha \\, m \\, e^{\\alpha x}\\right)
                (\\delta x)^2\\\\
            & + \\left[(x \\, \\delta m)^2 + (\\delta b)^2 +(m \\, \\delta x)^2\\right]

    where :math:`A`, :math:`\\alpha`, :math:`m`, and :math:`b` are the real
    constants to be fitted and :math:`x` is the independent variable.
    :math:`\\delta A`, :math:`\\delta \\alpha`, :math:`\\delta m`, :math:`\\delta b`,
    and :math:`\\delta x` are the respective uncertainties for :math:`A`,
    :math:`\\alpha`, :math:`m`, and :math:`b`, and :math:`x`.
    """
    _parameter_names = ("a", "alpha", "m", "b")

    def __init__(self):
        super().__init__()
        self._exponential = Exponential()
        self._linear = Linear()

    def __str__(self):
        exp_str = self._exponential.__str__().lstrip("f(x) = ")
        lin_str = self._linear.__str__().lstrip("f(x) = ")
        return f"f(x) = {exp_str} + {lin_str}"

    @property
    def latex_str(self) -> str:
        exp_str = self._exponential.latex_str
        lin_str = self._linear.latex_str
        return fr"{exp_str} + {lin_str}"

    @AbstractFitFunction.parameters.setter
    def parameters(self, val) -> None:
        AbstractFitFunction.parameters.fset(self, val)
        self._exponential.parameters = (self.parameters.a, self.parameters.alpha)
        self._linear.parameters = (self.parameters.m, self.parameters.b)

    @AbstractFitFunction.parameters_err.setter
    def parameters_err(self, val) -> None:
        AbstractFitFunction.parameters_err.fset(self, val)
        self._exponential.parameters_err = (
            self.parameters_err.a,
            self.parameters_err.alpha,
        )
        self._linear.parameters_err = (self.parameters_err.m, self.parameters_err.b)

    def func(self, x, a, alpha, m, b):
        exp_term = self._exponential.func(x, a, alpha)
        lin_term = self._linear.func(x, m, b)
        return exp_term + lin_term

    def _func_err(self, x, y, x_err=None):
        a, alpha, m, b = self.parameters

        exp_y, exp_err = self._exponential(x, x_err=x_err, reterr=True)
        lin_y, lin_err = self._linear(x, x_err=x_err, reterr=True)
        err = exp_err ** 2 + lin_err ** 2

        if x_err is not None:
            blend_err = 2 * a * alpha * m * np.exp(alpha * x) * (x_err ** x)
            err += blend_err

        return np.sqrt(err)


class ExponentialPlusOffset(AbstractFitFunction):
    """
    A sub-class of `AbstractFitFunction` to represent an exponential with a DC
    offset.

    .. math::

        y =& f(x) = A \\, e^{\\alpha \\, x} + m \\, x + b\\\\
        (\\delta y)^2 =&
            \\left( A e^{\\alpha x}\\right)^2 \\left[
                \\left( \\frac{\\delta A}{A} \\right)^2
                + (x \\, \\delta \\alpha)^2
                + (\\alpha \\, \\delta x)^2
            \\right]
            + (\\delta b)^2


    where :math:`A`, :math:`\\alpha`, and :math:`b` are the real constants to
    be fitted and :math:`x` is the independent variable.  :math:`\\delta A`,
    :math:`\\delta \\alpha`, :math:`\\delta b`, and :math:`\\delta x` are the
    respective uncertainties for :math:`A`, :math:`\\alpha`, and :math:`b`, and
    :math:`x`.

    """
    _parameter_names = ("a", "alpha", "b")

    def __init__(self):
        super().__init__()
        self._explin = ExponentialPlusLinear()

    def __str__(self):
        return f"f(x) = A exp(alpha x) + b"

    @property
    def latex_str(self) -> str:
        return fr"A \, \exp(B \, x) + C"

    @AbstractFitFunction.parameters.setter
    def parameters(self, val) -> None:
        AbstractFitFunction.parameters.fset(self, val)
        self._explin.parameters = (
            self.parameters.a,
            self.parameters.alpha,
            0.0,
            self.parameters.b,
        )

    @AbstractFitFunction.parameters_err.setter
    def parameters_err(self, val) -> None:
        AbstractFitFunction.parameters_err.fset(self, val)
        self._explin.parameters_err = (
            self.parameters_err.a,
            self.parameters_err.alpha,
            0.0,
            self.parameters_err.b,
        )

    def func(self, x, a, alpha, b):
        return self._explin.func(x, a, alpha, 0.0, b)


    def _func_err(self, x, y, x_err=None):
        _, err = self._explin(x, x_err=x_err, reterr=True)
        return err

    def root_solve(self, *args, **kwargs):
        """
        The root :math:`f(x_r) = 0` for the fit function.

        .. math::

            x_r &= \\frac{1}{\\alpha} \\ln \\left( \\frac{-b}{A} \\right)

            \\delta x_r &= \\sqrt{
                \\left( \\frac{1}{\\alpha} \\frac{\\delta A}{A} \\right)^2
                + \\left( x_r \\frac{\\delta \\alpha}{\\alpha} \\right)^2
                + \\left( \\frac{1}{\\alpha} \\frac{\\delta b}{b} \\right)^2
            }

        Parameters
        ----------
        *args
            Not needed.  This is to ensure signature comparability with
            `AbstractFitFunction`.

        *kwargs
            Not needed.  This is to ensure signature comparability with
            `AbstractFitFunction`.

        Returns
        -------
        root: float
            The root value for the given fit :attr:`parameters`.

        err: float
            The uncertainty in the calculated root for the given fit
            :attr:`parameters` and :attr:`parameters_err`.
        """
        a, alpha, b = self.parameters
        a_err, b_err, c_err = self.parameters_err

        root = np.log(-b / a) / alpha

        a_term = a_err / (a * alpha)
        b_term = b_err * root / alpha
        c_term = c_err / (alpha * b)
        err = np.sqrt(a_term ** 2 + b_term ** 2 + c_term ** 2)

        return root, err
