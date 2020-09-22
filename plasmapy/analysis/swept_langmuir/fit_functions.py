"""
`FitFunction` classes designed to assist in curve fitting of swept Langmuir
traces.
"""
__all__ = [
    "AbstractFitFunction",
    "ExponentialOffsetFitFunction",
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

    _parameter_names = ()  # type: Tuple[str, ...]
    _parameters = None  # type: Union[None, Tuple[Any, ...]]
    _parameters_err = None  # type: Union[None, Tuple[Any, ...]]
    _covariance_matrix = None
    _rsq = None
    _curve_fit_results = None

    def __init__(self):
        self.FitParamTuple = namedtuple("FitParamTuple", self._parameter_names)
        """
        A named tuple class used for attributes :attr:`parameters` and 
        :attr:`parameters_err`.  The attribute :attr:`parameter_names` defines
        the tuple field names.
        """

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
            (Default: `False`) If `True`, return an array of errors associated
            with the calculated independent variables

        Returns
        -------
        y: `numpy.ndarray`
            Corresponding dependent variables :math:`y=f(x)` of the independent
            variables :math:`x`.

        y_err: `numpy.ndarray`
            Errors associated with the calculated dependent variables
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

    @abstractmethod
    def _func(self, x, *args):
        """
        The fit function.

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
        Calculate dependent variable errors :math:`\\delta y` for dependent
        variables :math:`y=f(x)`.

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
            The calculated errors of the dependent variables of the independent
            variables `x`.
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
            The error associated with the root calculation.  **Currently this
            returns an array of** `numpy.nan` **values equal in shape to**
            `x` **, since there is no determined way to calculate the errors.**

        Notes
        -----
        If the full output of `scipy.optimize.fsolve` is desired then one can do

            >>> func = FitFunction()  # FitFunciton is a subclass of AbstractFitFunction
            >>> roots = scipy.optimize.fsolve(func, 0.)

        """
        kwargs["args"] = self.parameters
        results = fsolve(self._func, x0, **kwargs)
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
        popt, pcov = curve_fit(self._func, xdata, ydata, **kwargs)
        self._curve_fit_results = (popt, pcov)
        self._parameters = tuple(popt.tolist())
        self._parameters_err = tuple(np.sqrt(np.diag(pcov)).tolist())

        # calc rsq
        # rsq = 1 - (ss_res / ss_tot)
        residuals = ydata - self._func(xdata, *self.parameters)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
        self._rsq = 1 - (ss_res / ss_tot)


class ExponentialOffsetFitFunction(AbstractFitFunction):
    """
    A sub-class of `AbstractFitFunction` to represent an exponential with an
    offset.

    .. math::

        y &= f(x) = A \\, \\exp(B \\, x) + C

        (\\delta y)^2 &= (e^{B \\,x} \\delta A)^2
                         + (A \\, B \\, e^{B \\, x} \\delta B)^2
                         + (\\delta C)^2


    where :math:`A`, :math:`B`, and :math:`C` are positive real constants
    and :math:`x` is the independent variable.

    """
    _parameter_names = ("a", "b", "c")

    def __str__(self):
        return f"f(x) = A exp(B x) + C"

    def _func(self, x, a, b, c):
        """
        The fit function, an exponential with an offset.

        .. math::

            f(x) = A \\, \\exp(B \\, x) + C

        where :math:`A`, :math:`B`, and :math:`C` are positive real constants
        and :math:`x` is the independent variable.

        Parameters
        ----------
        x: array_like
            Independent variable.

        a: float
            value for constant :math:`A`

        b: float
            value for constant :math:`B`

        c: float
            value for constant :math:`C`

        Returns
        -------
        y: array_like
            dependent variables corresponding to :math:`x`

        """
        return a * np.exp(b * x) + c

    def _func_err(self, x, y, x_err=None):
        """
        Calculate dependent variable errors :math:`\\delta y` for dependent
        variables :math:`y=f(x)`.

        .. math::

            (\\delta y)^2 = (e^{B \\,x} \\delta A)^2
                             + (A \\, B \\, e^{B \\, x} \\delta B)^2
                             + (\\delta C)^2

        Parameters
        ----------
        x: array_like
            Independent variables to be passed to the fit function.

        Returns
        -------
        `numpy.ndarray`:
            The calculated errors of the dependent variables of the independent
            variables `x`.
        """
        a, b, c = self.parameters
        a_err, b_err, c_err = self.parameters_err
        a_term = (np.exp(b * x) * a_err) ** 2
        b_term = (a * b * np.exp(b * x) * b_err) ** 2
        c_term = c_err ** 2
        return np.sqrt(a_term + b_term + c_term)

    @property
    def latex_str(self) -> str:
        return fr"A \, \exp(B \, x) + C"

    def root_solve(self, *args, **kwargs):
        """
        The root :math:`f(x_r) = 0` for the fit function.

        .. math::

            x_r &= \\frac{1}{B} \\ln \\left( \\frac{-C}{A} \\right)

            \\delta x_r &= \\sqrt{
                \\left( \\frac{\\delta A}{A B} \\right)^2
                + \\left( x_r \\frac{\\delta B}{B} \\right)^2
                + \\left( \\frac{\\delta C}{B C} \\right)^2
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
            The error in the calculated root for the given fit
            :attr:`parameters` and :attr:`parameters_err`.
        """
        a, b, c = self.parameters
        root = np.log(-c / a) / b

        a_err, b_err, c_err = self.parameters_err
        a_term = a_err / (a * b)
        b_term = b_err * root / b
        c_term = c_err / (b * c)
        err = np.sqrt(a_term ** 2 + b_term ** 2 + c_term ** 2)

        return root, err


class Linear(AbstractFitFunction):
    """
    A sub-class of `AbstractFitFunction` to represent a linear function.

    .. math::

        y &= f(x) = m \\, x + b

        (\\delta y)^2 &= (x \\, \\delta m)^2 + (\\delta b)^2

    where :math:`m` and :math:`b` are positive real constants representing the
    slope and intercept, respectively, and :math:`x` is the independent
    variable.

    """

    _parameter_names = ("m", "b")

    def __str__(self):
        return f"f(x) = m x + b"

    def _func(self, x, m, b):
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
        Calculate dependent variable errors :math:`\\delta y` for dependent
        variables :math:`y=f(x)`.

        .. math::

            (\\delta y)^2 &= (x \\, \\delta m)^2 + (\\delta b)^2

        Parameters
        ----------
        x: array_like
            Independent variables to be passed to the fit function.

        Returns
        -------
        `numpy.ndarray`:
            The calculated errors of the dependent variables of the independent
            variables `x`.
        """
        m_term = (self.parameters_err[0] * x) ** 2
        b_term = self.parameters_err[1] ** 2
        return np.sqrt(m_term + b_term)

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
            The error in the calculated root for the given fit
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
        self._parameters = (m, b)

        m_err = results[4]
        b_err = np.sum(xdata ** 2) - ((np.sum(xdata) ** 2) / xdata.size)
        b_err = m_err * np.sqrt(1.0 / b_err)
        self._parameters_err = (m_err, b_err)

        self._rsq = results[2] ** 2
