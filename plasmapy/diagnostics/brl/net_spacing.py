"""The computational net spacing and it's inverse + derivative."""
import numpy as np


def get_s_points(num_points, s_upper_bound):
    r"""Get the equally spaced points for the variable `s`."""
    return np.linspace(0, s_upper_bound, num_points)


def large_probe_x_and_dx_ds(s_points, normalized_probe_radius):
    r"""The values of `x` and `dx/ds` that correspond to the `s_points` for a probe of large radius.

    Parameters
    ----------
    s_points : `numpy.ndarray`
    normalized_probe_radius : `float`
        The radio of probe radius to debye length, :math:`R_p / \lambda_D`.

    Returns
    -------
    x_points, dx_ds_points : `numpy.ndarray`
        Value at each :math:`s` of :math:`x(s)` and :math:`\frac{dx}{ds}(s)`.

    Notes
    -----
    Laframboise gives no explanation as to why this function is chosen. It can
    be found in appendix I: computer program listing, line 254 (page 3). This
    returns a function of the form
    .. math::

       x(s) = A (1 - s) - B (1 - s)^C

    where :math:`A`, :math:`B`, and :math:`C` are constants that depend on
    the `normalized_probe_radius`.
    """
    shared_constant = 0.5  # SNT in Laframboise
    linear_coefficient = 1 / (1 - shared_constant)  # ZNP in Laframboise
    power_coefficient = linear_coefficient - 1  # AYA in Laframboise
    power = (
        linear_coefficient - 10 / (shared_constant * normalized_probe_radius)
    ) / power_coefficient

    x_points = (
        linear_coefficient * (1 - s_points)
        - power_coefficient * (1 - s_points) ** power
    )
    dx_ds_points = linear_coefficient * -1 + power_coefficient * power * (
        1 - s_points
    ) ** (power - 1)
    return x_points, dx_ds_points


def medium_probe_x_and_dx_ds(s_points):
    r"""The values of `x` and `dx/ds` that correspond to the `s_points` for a probe of medium radius.

    Parameters
    ----------
    s_points : `numpy.ndarray`

    Returns
    -------
    x_points, dx_ds_points : `numpy.ndarray`
        Value at each :math:`s` of :math:`x(s)` and :math:`\frac{dx}{ds}(s)`.

    Notes
    -----
    Laframboise gives no explanation as to why this function is chosen. It can
    be found in appendix I: computer program listing, line 361 (page 3). This
    returns a function of the form
    .. math::

       x(s) = 1 - s.
    """
    x_points = 1 - s_points
    dx_ds_points = -1
    return x_points, dx_ds_points


def small_probe_x_and_dx_ds(s_points):
    r"""The values of `x` and `dx/ds` that correspond to the `s_points` for a probe of small radius.

    Parameters
    ----------
    s_points : `numpy.ndarray`

    Returns
    -------
    x_points, dx_ds_points : `numpy.ndarray`
        Value at each :math:`s` of :math:`x(s)` and :math:`\frac{dx}{ds}(s)`.

    Notes
    -----
    Laframboise gives no explanation as to why this function is chosen. It can
    be found in appendix I: computer program listing, line 360 (page 3). This
    returns a function of the form
    .. math::

       x(s) = e^{-s}
    """
    x_points = np.exp(-s_points)
    dx_ds_points = -np.exp(-s_points)
    return x_points, dx_ds_points
