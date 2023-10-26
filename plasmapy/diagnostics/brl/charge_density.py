"""Using the potential, calculate the charge density."""
import numpy as np

from scipy.integrate import quad_vec
from scipy.special import erfc, erfi


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
        # Equation (E.23).
        return -np.exp(-A) / np.pi**0.5 * ((A - chi) ** 0.5 + _g((A - chi) ** 0.5))
    else:
        # Equation (E.10).
        return 0


def eta_2(A, chi, chi_p, x, spherical=True):
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
    kappa = (chi - x**2 * chi_p) / (1 - x**2)
    if spherical:
        # Equation (E.25)
        return (
            -((1 - x**2) ** 0.5)
            * np.exp(-A)
            / np.pi**0.5
            * ((A - kappa) ** 0.5 + _g(A - kappa**0.5))
        )
    else:
        # TODO: Change this to use the approximations given by equations (E.45) - (E.65).

        # Variables defined in (E.45).
        tau = (kappa - chi_p) / 2
        mu = np.max(kappa, chi_p) - tau
        theta = chi - tau
        B = A - tau

        # Integrand from (E.59). (E.59) is a formula that is ready to integrate numerically.
        def integrand(omega):
            return 1 / (
                (-np.log(omega) - theta) * (np.log(omega) ** 2 - mu**2) ** 0.5
            )

        # Evaluate H1 from (E.59).
        H1 = (
            np.exp(-tau)
            * quad_vec(
                lambda omega: (
                    (-np.log(omega) - theta) * ((np.log(omega)) ** 2 - mu**2) ** 0.5
                )
                ** -1,
                0,
                np.exp(-B),
            )[0]
        )

        # Equation (E.44)
        return (
            np.exp(-A)
            / np.pi
            * np.arctan((x**2 / (1 - x**2) * (A - chi_p) / (A - kappa)) ** 0.5)
            + x / (1 - x**2) ** 0.5 * (chi_p - chi) / (2 * np.pi) * H1
        )


def eta_3(A, chi, x, x_B, spherical=True):
    r"""Calculate the contribution to the charge density when there is a finite radius at which the potential is 0.

    This only appears when the temperature of the attracted particles is zero.

    Parameters
    ----------
    A : `numpy.ndarray`
        Lower bound of the integral. Also, `A >= kappa` as defined in (E.3).
    chi : `numpy.ndarray`
        The normalized potential.
    x : `numpy.ndarray`
        Normalized inverse radius calculated from `~plasmapy.diagnostics.brl.net_spacing.get_x_and_dx_ds`.
    x_B : `float`
        Normalized inverse radius at which the potential is 0.
    spherical : `bool`, optional
        If `True` the probe will be treated as spherical. If `False` then the probe is cylindrical. Default is `True`.

    Returns
    -------
    `numpy.ndarray`

    Notes
    -----
    This is the term :math:`\eta_3(A)` that appears in various formula in appendix E.
    """
    # Equation (E.5).
    beta_B = -chi / (x**2 / x_B**2 - 1)
    if spherical:
        # Equation (E.27) where we've simplified the last integral.
        return -(((x**2 / x_B**2 - 1) / np.pi) ** 0.5) * (
            (beta_B - A) ** 0.5 * np.exp(-A)
            - np.exp(-beta_B) * np.pi**0.5 / 2 * erfi((beta_B - A) ** 0.5)
        )
    else:
        # TODO: See if this should be converted to using the approximations given in equations (E.66) - (E.85).
        # Equation (E.11).
        def integrand(beta):
            np.exp(-beta) * np.arcsin((x**2 / x_B**2 * beta / (beta - chi)) ** 0.5)

        return 1 / np.pi * quad_vec(integrand, A, beta_B)[0]


def eta_4(
        s1_indeces, 
        s2_indeces, 
        s_indeces, 
        ds, 
        chi, 
        dchi_ds, 
        x, 
        dx_ds, 
        eta_net, 
        gamma, 
        omega_G, 
        beta_G, 
        spherical=True
    ):
    r"""Calculate the contribution to the charge density across the locus of extrema.

    Parameters
    ----------
    s1_indeces, s2_indeces : `numpy.ndarray[float]`
        Lower and upper indeces of `s` to integrate between where `s1_indeces <= s2_indeces`.
    s_indeces : `numpy.ndarray[int]`
        The `s` index at which `eta_4` is calculated.
    ds : `float`
        Spacing between any two points in `s`.
    chi : `numpy.ndarray`
        The normalized potential.
    x : `numpy.ndarray`
        Normalized inverse radius calculated from `~plasmapy.diagnostics.brl.net_spacing.get_x_and_dx_ds`.
    eta_net : `numpy.ndarray`
        Normalized net charge density from equation (9.1).
    gamma : `float`
        :math:`\gamma = R_p^2 / \lambda_D^2` from equation (9.1).
    omega_G, beta_G : `numpy.ndarray`
        Normalized angular momentum and energy as functions of `s`.
    spherical : `bool`, optional
        If `True` the probe will be treated as spherical. If `False` then the probe is cylindrical. Default is `True`.

    Returns
    -------
    `numpy.ndarray`

    Notes
    -----
    This is the term :math:`\eta_4(\beta_1, \beta_2)` that appears in various 
    formula in appendix E. The integration is defined in the `CAL` function on 
    page 36 of the computer program listing, appendix I.
    """
    # TODO: Figure out how Laframboise takes care of the end point integration.
    num_s_points = chi.size
    s1_indeces_int = np.ceil(s1_indeces)
    s2_indeces_int = np.minimum(np.floor(s2_indeces), num_s_points - 1)
    max_num_s_prime_points = np.max(s2_indeces_int - s1_indeces_int)

    num_points = s1_indeces_int.size
    # Populate the integration matrix with equation (D.22) for each `s_index`.
    integration_matrix = np.zeros((num_points, max_num_s_prime_points))
    index_mask = np.zeros_like(integration_matrix)
    s_index_mask = np.zeros_like(integration_matrix)
    for point_index, (s1_index, s2_index, s_index) in enumerate(zip(s1_indeces_int, s2_indeces_int, s_indeces)):
        index_difference = s2_index - s1_index
        index_mask[point_index, :index_difference + 1] = np.arange(s1_index, s2_index + 1)
        s_index_mask[point_index, :index_difference + 1] = s_index

        if index_difference == 0:
            continue
        elif index_difference == 1:
            integration_matrix[point_index, 0] = 1
        elif index_difference == 2:
            integration_matrix[point_index, :index_difference + 1] = np.array([1, 4, 1]) / 3
        elif index_difference == 3:
            integration_matrix[point_index, :index_difference + 1] = np.array([3, 9, 9, 3]) / 8
        elif index_difference == 4:
            integration_matrix[point_index, :index_difference + 1] = np.array([9, 28, 22, 28, 9]) / 24
        else:
            integration_matrix[point_index, :index_difference + 1] = np.concatenate(
                (np.array([9, 28, 23]), 24 * np.ones(index_difference - 3), np.array([23, 28, 9]))
            ) / 24
    integration_matrix *= ds

    # Equation (E.30).
    # alpha_G = np.zeros_like(integration_matrix)
    # alpha_G[index_mask] = dchi_ds[index_mask] / 2 + gamma * eta_net[index_mask] / (2 * x[index_mask]**3) * dx_ds[index_mask]
    alpha_G = dchi_ds / 2 + gamma * eta_net / (2 * x**3) * dx_ds
    # Equation (E.31).
    # epsilon_G = np.zeros_like(integration_matrix)
    # epsilon_G[index_mask] = alpha_G[index_mask] * np.exp(-beta_G[index_mask])
    epsilon_G = alpha_G * np.exp(-beta_G)
    # Equation (E.32).
    # psi_G = np.zeros_like(integration_matrix)
    # psi_G[index_mask] = (beta_G[index_mask] - chi[s_index_mask] - omega_G[index_mask] * x[s_index_mask]**2)**0.5
    psi_G = (beta_G[:, np.newaxis] - chi[np.newaxis, :] - omega_G[:, np.newaxis] * x[np.newaxis, :]**2)**0.5

    integrand = epsilon_G[:, np.newaxis] * psi_G

    # Equation (E.33).
    # TODO: Fix this integration as it can't use the index_mask or s_index_mask.
    integral = -1 / np.pi**0.5 * np.sum(integration_matrix * integrand, axis=-1)

    # Now include the edges of the integral that we have not yet accounted for.
    # TODO: Use the methods of taylor expanding near the edge points for better inclusion of the edges during integration.

    # The non integer `s1_indeces > 0` because at `s[0] = 1` we are at the probe.
    non_integer_starts = True ^ np.isclose(s1_indeces, s1_indeces_int)
    # If the `integrand` is imaginary at the first index before `s1` then assume that the integration starts where the integrand is 0.
    # Otherwise just do a linear interpolation between `s1 - 1` and `s1`.
    integrand_at_start_index = integrand[s1_indeces_int[non_integer_starts]]
    integrand_before_start_index = integrand[s1_indeces_int[non_integer_starts] - 1]
    integrand_start = np.where(
        np.iscomplex(integrand_before_start_index),
        0,
        (integrand_before_start_index + integrand_at_start_index) / 2
    )
    # Use the trapezoidal rule to add this contribution.
    integral[non_integer_starts] += 0.5 * (integrand_start + integrand_at_start_index) * (s1_indeces_int - s1_indeces)[non_integer_starts] * ds

    non_integer_ends = True ^ np.isclose(s2_indeces, s2_indeces_int)
    integrand_end = np.zeros_like(non_integer_ends)
    # `s2` may be greater than the largest `s` in our grid. Thus we'll linearly interpolate the integrand using the two previous points.
    s2_larger_than_grid = s2_indeces > num_s_points - 1
    integrand_end[s2_larger_than_grid] = (integrand[-1, s2_larger_than_grid] - integrand[-2, s2_larger_than_grid]) / ds * \
        (s2_indeces[s2_larger_than_grid] - ds * (num_s_points - 1)) + integrand[-1, s2_larger_than_grid]

    # If the `integrand` is imaginary at the first index after `s2` then assume that the integration starts where the integrand is 0.
    # Otherwise just do a linear interpolation between `s2` and `s2 + 1`.
    integrand_at_end_index = integrand[s2_indeces_int[non_integer_ends]]
    integrand_after_end_index = integrand[s2_indeces_int[non_integer_ends] + 1]
    integrand_end = np.where(
        np.iscomplex(integrand_after_end_index),
        0,
        (integrand_after_end_index + integrand_at_end_index) / 2
    )
    # Use the trapezoidal rule to add this contribution.
    integral[non_integer_ends] += 0.5 * (integrand_end + integrand_at_end_index) * (s2_indeces - s2_indeces_int)[non_integer_ends] * ds

