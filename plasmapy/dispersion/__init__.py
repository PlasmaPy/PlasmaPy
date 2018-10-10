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
        if len(omega.shape) != 1:
            raise ValueError('"omega" must be a 1D array (got {}D)'.format(len(omega.shape)))
        if omega.shape[0] != k.shape[0]:
            raise ValueError(
                'The first dimension of "omega" and "k" must match (got {} and {} respectively)'.format(omega.shape[0], k.shape[0]))
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
