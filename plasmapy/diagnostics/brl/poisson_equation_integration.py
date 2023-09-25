"""Evaluation and integration of Poisson's equation."""


def evaluate_K0_term(eta_net, gamma, x, spherical=True):
    r"""Second derivative of normalized potential in terms of `x`.

    Parameters
    ----------
    eta_net : `numpy.ndarray`
        Normalized net charge density from equation (9.1).
    gamma : `float`
        :math:`\gamma = R_p^2 / \lambda_D^2` from equation (9.1).
    x : `numpy.ndarray`
        Normalized inverse radius calculated from `~plasmapy.diagnostics.brl.net_spacing.get_x_and_dx_ds`.
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


def evaluate_K1_term(K0_term, x, dx_ds, integration_matrix, spherical=True):
    r"""The :math:`s` dependent term of :math:`d\chi / ds`.

    Parameters
    ----------
    K0_term : `numpy.ndarray`
        The :math:`K_0(s)` term calculated from `~plasmapy.diagnostics.brl.poisson_equation_integration.evaluate_K0_term`.
    x, dx_ds : `numpy.ndarray`
        The the normalized inverse radius and it's `s` derivative calculated from `~plasmapy.diagnostics.brl.net_spacing.get_x_and_dx_ds`.
    integration_matrix : `numpy.ndarray`
        The matrix used for integrating functions from `~plasmapy.diagnostics.brl.integration.construct_integration_matrix`.
    spherical : `bool`, optional
        If `True` the probe will be treated as spherical. If `False` then the probe is cylindrical. Default is `True`.

    Returns
    -------
    `numpy.ndarray`

    Notes
    -----
    This is the term :math:`K_1(s)` as defined in equations (D.4) and (D.14) for the spherical and cylindrical probes respectively.
    """
    if spherical:
        return dx_ds * (integration_matrix @ (K0_term * dx_ds))
    else:
        return (dx_ds / x) * (integration_matrix @ (K0_term * dx_ds))


def evaluate_K2_term(K1_term, chi_0, integration_matrix):
    r"""The integral of :math:`K_1(s)` with the addition of the probe potential.

    Parameters
    ----------
    K1_term : `numpy.ndarray`
        The :math:`K_1(s)` term calculated from `~plasmapy.diagnostics.brl.poisson_equation_integration.evaluate_K1_term`.
    chi_0 : `float`
        The normalized potential of the probe.
    integration_matrix : `numpy.ndarray`
        The matrix used for integrating functions from `~plasmapy.diagnostics.brl.integration.construct_integration_matrix`.

    Returns
    -------
    `numpy.ndarray`

    Notes
    -----
    This is the term :math:`K_2(s)` as defined in equations (D.6). This is the same for both the spherical and cylindrical probes.
    """
    return chi_0 + integration_matrix @ K1_term
