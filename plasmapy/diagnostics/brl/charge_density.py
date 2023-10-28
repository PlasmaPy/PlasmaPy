"""Using the potential, calculate the charge density."""
import numpy as np
from diagnostics.brl.integration import integrate
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
        complex_integrand_at_s1,
        complex_integrand_at_s2,
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
    complex_integrand_at_s1, complex_integrand_at_s2 : `numpy.ndarray[bool]`
        Whether the integrand is complex at s1 or s2.
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
    All the terms that enter the integrand (`chi`, `x`, `eta_net`, `omega_G`, 
    and `beta_G`) should be 1D-arrays that contain the value for all points in 
    `s`.

    This is the term :math:`\eta_4(\beta_1, \beta_2)` that appears in various 
    formula in appendix E. The integration is defined in the `CAL` function on 
    page 36 of the computer program listing, appendix I.
    """
    if spherical:
        # Equation (E.30).
        alpha_G = dchi_ds / 2 + gamma * eta_net / (2 * x**3) * dx_ds
        # Equation (E.31).
        epsilon_G = alpha_G * np.exp(-beta_G)
        # Equation (E.32).
        psi_G = (beta_G[np.newaxis, :] - chi[:, np.newaxis] - omega_G[np.newaxis, :] * x[:, np.newaxis]**2)**0.5

        all_integrands = -1 / np.pi**0.5 * epsilon_G[np.newaxis, :] * psi_G # Indexed as `[s, s']` in (E.33).
    else:
        # Equation (E.86).
        alpha_G = dchi_ds + gamma * eta_net / (2 * x**3) * dx_ds
        # Equation (E.31).
        epsilon_G = alpha_G * np.exp(-beta_G)
        # Equation (E.87).
        psi_G = np.arctan((omega_G[np.newaxis, :] * x[:, np.newaxis]**2 / (beta_G[np.newaxis, :] - chi[:, np.newaxis] - omega_G[np.newaxis, :] * x[:, np.newaxis]**2))**0.5)

        all_integrands = 1 / np.pi * epsilon_G[np.newaxis, :] * psi_G
    all_integrands = all_integrands[s_indeces]
    
    return integrate(all_integrands, s1_indeces, s2_indeces, complex_integrand_at_s1, complex_integrand_at_s2, ds)

