"""
Implementation of the magnetized transport coefficients presented by
Epperlein and Haines in their 1986 paper "Plasma transport coefficients in a magnetic field by direct numerical solution of
the Fokker–Planck equation" (Phys. Fluids), doi: 10.1063/1.865901
"""

__all__ = [
    "EpperleinHainesPolynomialFit",
    "EpperleinHainesInterpolated",
]

# TODO: allow options for returning transport coefficients in physical units
# this will require more information about the plasma
# look at sec. IIB in Ridgers2008 (pg 3) for relations: these likely hold for EH as well
# EH equivalent is right at the end of section II (pg. 5). Looks the same.

import numpy as np
import os
import warnings

from plasmapy.transport.classical.base import (
    AbstractClassicalTransportCoefficients,
    AbstractInterpolatedCoefficients,
    validate_object,
)

# Get the absolute path to the data files
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
# Load the table of coefficients
c = np.load(os.path.join(data_dir, "epperlein_haines_polynomial_fit_coefficients.npz"))


def _find_nearest(Z):
    """
    Finds the nearest Z-value in the EH tables for a given Z. Prints a
    warning if the Z found is not equal to the Z requested.

    Parameters
    ----------
    Z : float
        An integer charge

    Returns
    -------
    i : int
        The index of the closest Z in the EH tables

    """
    if Z == np.inf:
        return -1

    i = np.argmin(np.abs(c["Z"] - Z))
    if c["Z"][i] != Z:
        warnings.warn(
            f"Value Z = {Z} is not in the Epperlein-Haines table. "
            f"Using the nearest value, Z = {c['Z'][i]}. "
            f"The values in the table are {c['Z']}.",
            RuntimeWarning,
        )
    return i


class EpperleinHainesPolynomialFit(AbstractClassicalTransportCoefficients):
    @property
    def _norm_alpha_para(self):
        i = _find_nearest(self.Z)
        return c["alpha0"][i] * np.ones(self.chi_e.size)

    @property
    def _norm_alpha_perp(self):
        i = _find_nearest(self.Z)
        return 1 - (c["alpha1p"][i] * self.chi_e + c["alpha0p"][i]) / (
            self.chi_e ** 2 + c["a1p"][i] * self.chi_e + c["a0p"][i]
        )

    @property
    def _norm_alpha_cross(self):
        i = _find_nearest(self.Z)
        return (
            self.chi_e
            * (c["alpha1pp"][i] * self.chi_e + c["alpha0pp"][i])
            / (
                self.chi_e ** 3
                + c["a2pp"][i] * self.chi_e ** 2
                + c["a1pp"][i] * self.chi_e
                + c["a0pp"][i]
            )
            ** (8 / 9)
        )

    @property
    @validate_object(properties=["chi_e", "Z"])
    def norm_alpha(self):
        """
        Calculates the normalized alpha coefficients in terms of the
        dimensionless Hall parameter and the ionization fraction.

        Parameters
        ----------
        chi : float (N,)
            The dimensionless hall parameter (ratio of the electron gyrofrequency
            and the electron-ion collision frequency).

        Z : float
            Ionization fraction. The value will be coerced to the nearest value
            from the Epperlein-Haines tables:
            Z = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 20, 30, 60, ∞]

        Returns
        -------
        norm_alpha, `~numpy.ndarray` of `~u.Quantity` instances (3, N)
            The resistivity coefficients:
                [alpha_para, alpha_perp, alpha_cross]

        Notes
        -----
        For details, see "Plasma transport coefficients in a magnetic field by
        direct numerical solution of the Fokker–Planck equation" by
        E. M. Epperlein and M. G. Haines, DOI: `10.1063/1.865901`

        .. _`10.1063/1.865901`: https://aip.scitation.org/doi/10.1063/1.865901

        """
        return np.array(
            [self._norm_alpha_para, self._norm_alpha_perp, self._norm_alpha_cross]
        )

    @property
    def _norm_beta_para(self):
        i = _find_nearest(self.Z)
        return c["beta0"][i] * np.ones(self.chi_e.size)

    @property
    def _norm_beta_perp(self):
        i = _find_nearest(self.Z)
        return (c["beta1p"][i] * self.chi_e + c["beta0p"][i]) / (
            self.chi_e ** 3
            + c["b2p"][i] * self.chi_e ** 2
            + c["b1p"][i] * self.chi_e
            + c["b0p"][i]
        ) ** (8 / 9)

    # TODO: Note that this function doesn't match the EH paper Fig. 1 in the
    # chi -> inf side. The coefficients and polynomial are right...
    # this might be a mistake in the EH tables?
    @property
    def _norm_beta_cross(self):
        i = _find_nearest(self.Z)
        return (
            self.chi_e
            * (c["beta1pp"][i] * self.chi_e + c["beta0pp"][i])
            / (
                self.chi_e ** 3 * c["b2pp"][i] * self.chi_e ** 2
                + c["b1pp"][i] * self.chi_e
                + c["b0pp"][i]
            )
        )

    @property
    @validate_object(properties=["chi_e", "Z"])
    def norm_beta(self):
        """
        Calculates the normalized beta coefficients in terms of the
        dimensionless Hall parameter and the ionization fraction.

        Parameters
        ----------
        chi : float (N,)
            The dimensionless hall parameter (ratio of the electron gyrofrequency
            and the electron-ion collision frequency).

        Z : float
            Ionization fraction. The value will be coerced to the nearest value
            from the Epperlein-Haines tables:
            Z = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 20, 30, 60, ∞]

        Returns
        -------
        norm_beta, `~numpy.ndarray` of `~u.Quantity` instances (3, N)
            The thermoelectric coefficients:
                [beta_para, beta_perp, beta_cross]

        Notes
        -----
        For details, see "Plasma transport coefficients in a magnetic field by
        direct numerical solution of the Fokker–Planck equation" by
        E. M. Epperlein and M. G. Haines, DOI: `10.1063/1.865901`

        .. _`10.1063/1.865901`: https://aip.scitation.org/doi/10.1063/1.865901

        """
        return np.array(
            [self._norm_beta_para, self._norm_beta_perp, self._norm_beta_cross]
        )

    @property
    def _norm_kappa_e_para(self):
        i = _find_nearest(self.Z)
        return c["gamma0"][i] * np.ones(self.chi_e.size)

    @property
    def _norm_kappa_e_perp(self):
        i = _find_nearest(self.Z)
        return (c["gamma1p"][i] * self.chi_e + c["gamma0p"][i]) / (
            self.chi_e ** 3
            + c["c2p"][i] * self.chi_e ** 2
            + c["c1p"][i] * self.chi_e
            + c["c0p"][i]
        )

    @property
    def _norm_kappa_e_cross(self):
        i = _find_nearest(self.Z)
        return (
            self.chi_e
            * (c["gamma1pp"][i] * self.chi_e + c["gamma0pp"][i])
            / (
                self.chi_e ** 3
                + c["c2pp"][i] * self.chi_e ** 2
                + c["c1pp"][i] * self.chi_e
                + c["c0pp"][i]
            )
        )

    @property
    @validate_object(properties=["chi_e", "Z"])
    def norm_kappa_e(self):
        """
        Calculates the normalized kappa_e coefficients in terms of the
        dimensionless Hall parameter and the ionization fraction.

        Parameters
        ----------
        chi : float (N,)
            The dimensionless hall parameter (ratio of the electron gyrofrequency
            and the electron-ion collision frequency).

        Z : float
            Ionization fraction. The value will be coerced to the nearest value
            from the Epperlein-Haines tables:
            Z = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 20, 30, 60, ∞]

        Returns
        -------
        norm_kapp_e, `~numpy.ndarray` of `~u.Quantity` instances (3, N)
            The electron thermal conductivity coefficients:
                [kappa_e_para, kappa_e_perp, kappa_e_cross]

        Notes
        -----
        For details, see "Plasma transport coefficients in a magnetic field by
        direct numerical solution of the Fokker–Planck equation" by
        E. M. Epperlein and M. G. Haines, DOI: `10.1063/1.865901`

        .. _`10.1063/1.865901`: https://aip.scitation.org/doi/10.1063/1.865901

        """
        return np.array(
            [self._norm_kappa_e_para, self._norm_kappa_e_perp, self._norm_kappa_e_cross]
        )


class EpperleinHainesInterpolated(AbstractInterpolatedCoefficients):
    @property
    def _data_file(self):
        return os.path.join(data_dir, "epperlein_haines_data.npz")


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    chi = np.linspace(-2, 2, num=50)
    chi = 10 ** chi

    coef = EpperleinHainesInterpolated(chi_e=chi, Z=5)

    data = coef.norm_alpha[1, :]

    # print(data)

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.plot(chi, data)

    chi = np.linspace(-2, 2, num=100)
    chi = 10 ** chi

    # Instantiate the object
    coef1 = EpperleinHainesPolynomialFit(chi_e=chi, Z=1)
    coef2 = EpperleinHainesPolynomialFit(chi_e=chi, Z=np.inf)

    fig, axarr = plt.subplots(nrows=3, ncols=2, figsize=(10, 10), sharex=True)

    for i in range(3):
        for j in range(2):
            ax = axarr[i][j]
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(1e-2, 1e2)

    axarr[-1][0].set_xlabel("$\\chi$")
    axarr[-1][1].set_xlabel("$\\chi$")

    ax = axarr[0][0]
    ax.set_title("$\\alpha_\\perp$")
    ax.set_ylim(0.2, 1)
    ax.plot(chi, coef1.norm_alpha[1], label="Z=1")
    ax.plot(chi, coef2.norm_alpha[1], label="Z=Infinite")
    ax.legend()

    ax = axarr[1][0]
    ax.set_title("$\\alpha_\\wedge$")
    ax.set_ylim(0.01, 1)
    ax.plot(chi, coef1.norm_alpha[2])
    ax.plot(chi, coef2.norm_alpha[2])

    ax = axarr[2][0]
    ax.set_title("$\\beta_\\perp$")
    ax.set_ylim(0.001, 10)
    ax.plot(chi, coef1.norm_beta[1])
    ax.plot(chi, coef2.norm_beta[1])

    ax = axarr[0][1]
    ax.set_title("$\\beta_\\wedge$")
    ax.set_ylim(0.01, 1)
    ax.plot(chi, coef1.norm_beta[2])
    ax.plot(chi, coef2.norm_beta[2])

    ax = axarr[1][1]
    ax.set_title("$\\kappa_\\perp$")
    ax.set_ylim(0.001, 100)
    ax.plot(chi, coef1.norm_kappa_e[1])
    ax.plot(chi, coef2.norm_kappa_e[1])

    ax = axarr[2][1]
    ax.set_title("$\\kappa_\\wedge$")
    ax.set_ylim(0.01, 10)
    ax.plot(chi, coef1.norm_kappa_e[2])
    ax.plot(chi, coef2.norm_kappa_e[2])

    """
    import astropy.units as u
    chi = np.linspace(-2, 2, num=5)
    chi = 10 ** chi
    coef = EpperleinHainesPolynomialFit(chi_e=chi, Z = 1)

    print(coef.norm_alpha.shape)
    """
