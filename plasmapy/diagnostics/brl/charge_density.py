"""Using the potential, calculate the charge density."""
import numpy as np

from scipy.special import erfc


def eta_1(A, chi, spherical=True):
    r"""Calculate the contribution to the charge density between the region where no particles exist and non-striking particles exist.

    Parameters
    ----------
    A : `numpy.ndarray`
        Lower bound of the integral. Also, `A >= chi`.
    chi : `numpy.ndarray`
        The normalized potential.
    spherical : `bool`, optional
        If `True` the probe will be treated as spherical. If `False` then the probe is cylindrical. Default is `True`.

    Returns
    -------
    `numpy.ndarray`

    Notes
    -----
    This is the term :math:`\eta_1(A)` that appears in various formula in appendix E.
    """
    if spherical:
        # Equation (E.21).
        g = np.pi**0.5 / 2 * np.exp(A - chi) * erfc((A - chi) ** 0.5)
        # Equation (E.23).
        return -np.exp(-A) / np.pi**0.5 * ((A - chi) ** 0.5 + g)
    else:
        # Equation (E.10).
        return 0
