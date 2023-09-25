"""Evaluation and integration of Poisson's equation."""


def potential_second_derivative_x(eta_net, gamma, x, spherical=True):
    r"""Second derivative of normalized potential in terms of x.

    Parameters
    ----------
    eta_net : `numpy.ndarray`
        Normalized net charge density from equation (9.1).
    gamma : `float`
        :math:`\gamma = R_p^2 / \lambda_D^2` from equation (9.1).
    x : `numpy.ndarray`
        Normalized inverse radius.
    spherical : `bool`, optional
        If `True` the probe will be treated as spherical. If `False` then the probe is cylindrical. Default is `True`.

    Returns
    -------
    `numpy.ndarray`

    Notes
    -----
    This is the term :math:`K_0(s)` as defined in equations (D.1) and (D.12) for the spherical and cylindrical probes respectively. The charge density, `eta_net`, and inverse radius, `x`, should have the same shape and correspond to the same points in space.
    """
    if spherical:
        return -gamma * eta_net / x**4
    else:
        return -gamma * eta_net / x**3
