r"""
Gyrokinetic (Maxwellian) dispersion relation solver in normalized variables.

Implements a normalized gyrokinetic dispersion relation in terms of:

.. math::

    \bar{\omega} &= \omega / (k_{\parallel} v_A), \\
    k_{\perp} \rho_i, \\
    \beta_i, \\
    \tau &= T_i / T_e,

with an optional ion/electron mass ratio.

This module provides:

- `gyrokinetic_dispersion_residual`: complex residual of the dispersion
  relation
- `solve_gyrokinetic_dispersion`: complex root solve for a single
  ``k_perp_rho_i``
- `solve_gyrokinetic_dispersion_spectrum`: continuation solve over a
  ``k_perp_rho_i`` grid
"""

from __future__ import annotations

__all__ = [
    "gyrokinetic_dispersion_residual",
    "solve_gyrokinetic_dispersion",
    "solve_gyrokinetic_dispersion_spectrum",
]

import numpy as np
from scipy.optimize import root
from scipy.special import ive, wofz


def gyrokinetic_dispersion_residual(
    omega: complex,
    k_perp_rho_i: float,
    beta_i: float,
    tau: float,
    mass_ratio: float = 1836.15267343,
) -> complex:
    r"""
    Return the complex residual of the normalized gyrokinetic dispersion relation.

    Parameters
    ----------
    omega : complex
        Normalized frequency :math:`\bar{\omega} = \omega / (k_\parallel v_A)`.

    k_perp_rho_i : float
        Perpendicular wavenumber normalized to the ion Larmor radius.

    beta_i : float
        Ion plasma beta.

    tau : float
        Ion-to-electron temperature ratio :math:`T_i/T_e`.

    mass_ratio : float, optional
        Ion-to-electron mass ratio (m_i/m_e). Default is 1836.15267343.

    Returns
    -------
    complex
        Complex residual. Roots correspond to solutions of the dispersion relation.
    """
    # Sanity checks
    if beta_i <= 0:
        raise ValueError("beta_i must be > 0.")
    if tau <= 0:
        raise ValueError("tau must be > 0.")
    if k_perp_rho_i < 0:
        raise ValueError("k_perp_rho_i must be >= 0.")
    if mass_ratio <= 0:
        raise ValueError("mass_ratio must be > 0.")

    # Normalized electron Larmor-radius ratio
    rho_e = np.sqrt(mass_ratio * tau)

    alpha_i = 0.5 * k_perp_rho_i**2
    alpha_e = alpha_i / rho_e**2

    # xi = omega_bar / sqrt(beta_i) (ion), then electron scaled by rho_e
    xi_i = omega / np.sqrt(beta_i)
    xi_e = xi_i / rho_e

    # Plasma dispersion function Z(xi)
    zeta_i = 1j * np.sqrt(np.pi) * wofz(xi_i)
    zeta_e = 1j * np.sqrt(np.pi) * wofz(xi_e)

    # Modified Bessel functions (scaled)
    ive0_i = ive(0, alpha_i)
    ive1_i = ive(1, alpha_i)
    ive0_e = ive(0, alpha_e)
    ive1_e = ive(1, alpha_e)

    G0_i = ive0_i * xi_i * zeta_i
    G0_e = ive0_e * xi_e * zeta_e
    G1_i = (ive0_i - ive1_i) * xi_i * zeta_i
    G1_e = (ive0_e - ive1_e) * xi_e * zeta_e

    A = 1 + G0_i + tau * (1 + G0_e)
    B = 1 - ive0_i + tau * (1 - ive0_e)
    C = G1_i - G1_e
    D = G1_i + G1_e / tau
    E = ive0_i - ive1_i - ive0_e + ive1_e

    lhs = (alpha_i * A / omega**2 - A * B + B**2) * (2 * A / beta_i - 2 * A * D + C**2)
    rhs = (A * E + B * C) ** 2

    return complex(lhs - rhs)


def solve_gyrokinetic_dispersion(
    k_perp_rho_i: float,
    beta_i: float,
    tau: float,
    omega_guess: complex = 1.0 + 0.0j,
    mass_ratio: float = 1836.15267343,
    tol: float = 1e-12,
    methods: tuple[str, ...] = ("hybr", "lm"),
) -> complex:
    """
    Solve the normalized gyrokinetic dispersion relation for a single k_perp_rho_i.

    Solves a 2D real system (Re(residual)=0, Im(residual)=0) using scipy.optimize.root
    with optional fallback across methods.

    Returns NaN+1j*NaN if no method succeeds.
    """

    def residual_wrapped(w: np.ndarray) -> np.ndarray:
        omega = w[0] + 1j * w[1]
        res = gyrokinetic_dispersion_residual(
            omega, k_perp_rho_i, beta_i, tau, mass_ratio=mass_ratio
        )
        return np.array([res.real, res.imag], dtype=float)

    x0 = np.array([omega_guess.real, omega_guess.imag])

    for method in methods:
        sol = root(residual_wrapped, x0, method=method, tol=tol)
        if sol.success and np.all(np.isfinite(sol.x)):
             return complex(sol.x[0] + 1j * sol.x[1])

    return np.nan + 1j * np.nan


def solve_gyrokinetic_dispersion_spectrum(
    k_perp_rho_i: np.ndarray,
    beta_i: float,
    tau: float,
    initial_guess: complex = 1.0 + 0.0j,
    mass_ratio: float = 1836.15267343,
    tol: float = 1e-12,
) -> np.ndarray:
    """
    Solve omega_bar(k_perp_rho_i) over an iterable of k values using continuation.

    The previous solution is used as the next initial guess when finite.
    """
    omega_out = []
    guess = initial_guess

    for ik in k_perp_rho_i:
        omega = solve_gyrokinetic_dispersion(ik, beta_i, tau, guess, mass_ratio, tol)
        omega_out.append(omega)

        if np.isfinite(omega.real) and np.isfinite(omega.imag):
            guess = omega

    return np.asarray(omega_out, dtype=np.complex128)
