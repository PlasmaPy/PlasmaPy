"""
Tools for working with dispersion solvers.
"""

import abc


class DispersionInput(abc.ABC):
    """
    Base class that stores the input to a dispersion solver.
    """


class DispersionOutput(abc.ABC):
    """
    Base class that stores the output to a dispersion solver.

    Parameters
    ----------
    input : `DispersionInput`
        The input parameters passed to the dispersion solver that resulted in
        this output.
    k : array
        Wavevectors.
    omega : array
        Angular frequencies. Can be complex, and must be the same size as *k*.
    """
    def __init__(self, input, omega, k):
        self._omega = omega
        self._k = k
        self._input = input

    @property
    def input(self):
        """
        Input parameters to the dispersion relation.
        """
        return self._input

    @property
    def k(self):
        """
        Wavevectors.
        """
        return self._k

    @property
    def omega(self):
        """
        Angular freuencies.
        """
        return self._omega


class DispersionSolver(abc.ABC):
    """
    Base class that stores a dispersion solver. Must be overriden by specific
    dispersion solver implementations.
    """
    @abc.abstractmethod
    def solve(self, input):
        """
        Run the dispersion solver for a given input.

        Parameters
        ----------
        input : DispersionInput

        Returns
        -------
        output : DispersionOutput
        """
        pass
