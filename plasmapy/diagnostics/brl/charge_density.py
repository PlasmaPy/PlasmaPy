"""Using the potential, calculate the charge density."""
import numpy as np

from scipy.special import erfc


def _g(x):
    r"""Calculate :math:`g(x)`.

    This is the function defined in (E.21) as
    .. math::

        g(x) = \frac{\sqrt{\pi}}{2} e^{x^2} \erfc(x)

    Parameters
    ----------
    x : `numpy.ndarray`

    Returns
    -------
    `numpy.ndarray`
    """
    return np.pi**0.5 / 2 * np.exp(x**2) * erfc(x)


def eta_1(A, chi, spherical=True):
    r"""Calculate the contribution to the charge density along the :math:`J^2 = 0` line.

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


def eta_2(A, chi, x, spherical=True):
    r"""Calculate the contribution to the charge density between the region where particles are absorbed by the probe and non-striking particles exist.

    Parameters
    ----------
    A : `numpy.ndarray`
        Lower bound of the integral. Also, `A >= kappa` as defined in (E.3).
    chi : `numpy.ndarray`
        The normalized potential.
    chi_p : `float`
        The normalized potential of the probe.
    x : `numpy.ndarray`
        Normalized inverse radius calculated from `~plasmapy.diagnostics.brl.net_spacing.get_x_and_dx_ds`.
    spherical : `bool`, optional
        If `True` the probe will be treated as spherical. If `False` then the probe is cylindrical. Default is `True`.

    Returns
    -------
    `numpy.ndarray`

    Notes
    -----
    This is the term :math:`\eta_2(A)` that appears in various formula in appendix E.
    """
    if spherical:
        # Equation (E.3)
        pass
