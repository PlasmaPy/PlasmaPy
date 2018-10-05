"""
Tools for working with dispersion solvers.
"""


class DispersionInput:
    """
    Base class that stores the input to a dispersion solver.
    """


class DispersionOutput:
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
        if omega.size != k.size:
            raise ValueError(
                'Arguments "omega" (size {}) and "k" (size {}) must have the same size.'.format(omega.size, k.size))
        self._omega = omega
        self._k = k
        self._input = input

    @property
    def input(self):
        """
        Input parameters
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
