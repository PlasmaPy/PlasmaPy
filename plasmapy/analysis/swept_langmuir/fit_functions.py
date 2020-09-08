"""
`FitFunction` classes designed to assist in curve fitting of swept Langmuir
traces.
"""
__all__ = [
    "AbstractFitFunction",
    "ExponentialOffsetFitFunction",
    "LinearFitFunction",
]
import numpy as np

from abc import ABC, abstractmethod
from scipy.stats import linregress
from scipy.optimize import curve_fit, fsolve
from typing import Tuple


class AbstractFitFunction(ABC):
    """
    Abstract class for defining fit functions and the tools for fitting the
    function to a set of data.  These were originally designed for assisting in
    fitting curves to swept langmuir data.
    """
    _parameters = None  # type: Tuple
    _parameters_err = None  # type: Tuple
    _covariance_matrix = None
    _rsq = None

    def __call__(self, x):
        """
        The fit function.

        Parameters
        ----------
        x: array_like
            Dependent variables.

        Returns
        -------
        array_like
            Corresponding independent variables of dependent variables `x`.
        """
        return self._func(x, *self.parameters)

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
        array:
            The calculated dependent variables of the independent variables `x`.
        """
        raise NotImplementedError

    @property
    def parameters(self) -> Tuple:
        """The fitted parameters for the fit function."""
        return self._parameters

    @property
    def parameters_err(self) -> Tuple:
        """The associated errors of the fit `parameters`."""
        return self._parameters_err

    @property
    @abstractmethod
    def latex_str(self) -> str:
        """Latex friendly representation of the fit function."""
        raise NotImplementedError

    def root_solve(self, x0, **kwargs):
        """
        Solve for the root of the fit function (i.e. `f(x) = 0`).

        Parameters
        ----------
        x0: `~np.ndarray`
            The starting estimate for the roots of `f(x) = 0`.

        **kwargs
            Any keyword accepted by `scipy.optimize.fsolve`, except for `args`.

        Returns
        -------
        x : ndarray
            The solution (or the result of the last iteration for an
            unsuccessful call).
        infodict : dict
            A dictionary of optional outputs with the keys:

            ``nfev``
                number of function calls
            ``njev``
                number of Jacobian calls
            ``fvec``
                function evaluated at the output
            ``fjac``
                the orthogonal matrix, q, produced by the QR
                factorization of the final approximate Jacobian
                matrix, stored column wise
            ``r``
                upper triangular matrix produced by QR factorization
                of the same matrix
            ``qtf``
                the vector ``(transpose(q) * fvec)``

        ier : int
            An integer flag.  Set to 1 if a solution was found, otherwise refer
            to `mesg` for more information.
        mesg : str
            If no solution is found, `mesg` details the cause of failure.
        """
        kwargs["args"] = self.parameters
        return fsolve(self._func, x0, **kwargs)

    @property
    def rsq(self):
        """Coefficient of determination (r-squared) value of the fit."""
        return self._rsq

    def curve_fit(self, xdata, ydata, **kwargs):
        """
        Use non-linear least squares to fit the fit function to (`xdata`, `ydata`).
        This uses `scipy.optimize.curve_fit`.

        Parameters
        ----------
        xdata: array_like
            The independent variable where data is measured.  Should be 1D of
            length M.

        ydata: array_like
            The dependent data associated with `xdata`.

        **kwargs
            Any keywords accepted by `scipy.optimize.curve_fit`.

        Returns
        -------
        popt : array
            Optimal values for the parameters so that the sum of the squared
            residuals of ``f(xdata, *popt) - ydata`` is minimized.

        pcov : 2-D array
            The estimated covariance of popt. The diagonals provide the variance
            of the parameter estimate. To compute one standard deviation errors
            on the parameters use ``perr = np.sqrt(np.diag(pcov))``.
            How the `sigma` parameter affects the estimated covariance
            depends on `absolute_sigma` argument, as described above.
            If the Jacobian matrix at the solution doesn't have a full rank, then
            'lm' method returns a matrix filled with ``np.inf``, on the other hand
            'trf'  and 'dogbox' methods use Moore-Penrose pseudoinverse to compute
            the covariance matrix.

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
        self._parameters = tuple(popt.tolist())
        self._parameters_err = tuple(np.sqrt(np.diag(pcov)).tolist())

        # calc rsq
        # rsq = 1 - (ss_res / ss_tot)
        residuals = ydata - self._func(xdata, *self.parameters)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
        self._rsq = 1 - (ss_res / ss_tot)

        return popt, pcov


class ExponentialOffsetFitFunction(AbstractFitFunction):
    def __str__(self):
        return f"f(x) = A exp(B x) + C"

    def _func(self, x, a, b, c):
        return a * np.exp(b * x) + c

    @property
    def latex_str(self) -> str:
        return fr"A \, \exp(B \, x) + C"

    def root_solve(self, *args, reterr=True, **kwargs):
        a, b, c = self.parameters
        root = np.log(-c / a) / b

        if reterr:
            a_err, b_err, c_err = self.parameters_err
            a_term = a_err / (a * b)
            b_term = b_err * root / b
            c_term = c_err / (b * c)
            err = np.sqrt(a_term ** 2 + b_term ** 2 + c_term ** 2)

            return root, err
        return root


class LinearFitFunction(AbstractFitFunction):
    def __str__(self):
        return f"f(x) = m x + b"

    def _func(self, x, m, b):
        return m * x + b

    @property
    def latex_str(self) -> str:
        return fr"m \, x + b"

    def root_solve(self, *args, reterr=True, **kwargs):
        m, b = self.parameters
        root = -b / m

        if reterr:
            m_err, b_err = self.parameters_err
            err = np.abs(root) * np.sqrt((m_err / m) ** 2 + (b_err / b) ** 2)

            return root, err
        return root

    def curve_fit(self, xdata, ydata, **kwargs):
        results = linregress(xdata, ydata)

        m = results[0]
        b = results[1]
        self._parameters = (m, b)

        m_err = results[4]
        b_err = np.sum(xdata ** 2) - ((np.sum(xdata) ** 2) / xdata.size)
        b_err = m_err * np.sqrt(1.0 / b_err)
        self._parameters_err = (m_err, b_err)

        self._rsq = results[2] ** 2
        return results
