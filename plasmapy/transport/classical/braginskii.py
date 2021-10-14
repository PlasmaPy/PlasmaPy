"""
Implementation of the magnetized transport coefficients presented by
Braginskii in his 1965 review article "Transport Processes in a Plasma"
"""

__all__ = [
    "BraginskiiPolynomialFit",
    "BraginskiiInterpolated",
]

import numpy as np
import os
import warnings

from plasmapy.transport.classical.base import (
    AbstractInterpolatedCoefficients,
    AbstractPolynomialCoefficients,
    validate_object,
)

# Get the absolute path to the data files
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
# Load the table of coefficients
coeff_table = np.load(
    os.path.join(data_dir, "braginskii_polynomial_fit_coefficients.npz")
)

# TODO: Note that the Braginskii Beta coefficient is normalized differently than in EH:
# Compare Braginskii Eq 4.30-4.33 to EH Eq. 7a-7b. Factors of ne and Te?
# The docs need to include a very clear indication that of what these coefficients
# are!
# The EH normalization is at the top of Pg.5 right
# The Braginskii normalization is just implicit in the equtions


class BraginskiiPolynomialFit(AbstractPolynomialCoefficients):
    @property
    def _c(self):
        return coeff_table

    def _Delta(self, i):
        return (
            self.chi_e ** 4
            + self._c["delta1"][i] * self.chi_e ** 2
            + self._c["delta0"][i]
        )

    @property
    def _norm_alpha_para(self):
        i = self._find_nearest_Z(self.Z)
        return self._c["alpha0"][i] * np.ones(self.chi_e.size)

    @property
    def _norm_alpha_perp(self):
        i = self._find_nearest_Z(self.Z)
        return 1 - (
            self._c["alpha1p"][i] * self.chi_e ** 2 + self._c["alpha0p"][i]
        ) / self._Delta(i)

    @property
    def _norm_alpha_cross(self):
        i = self._find_nearest_Z(self.Z)
        return (
            self.chi_e
            * (self._c["alpha1pp"][i] * self.chi_e ** 2 + self._c["alpha0pp"][i])
            / self._Delta(i)
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
            from the Braginskii tables:
            Z = [1, 2, 3, 4, ∞]

        Returns
        -------
        norm_alpha, `~numpy.ndarray` of `~u.Quantity` instances (3, N)
            The resistivity coefficients:
                [alpha_para, alpha_perp, alpha_cross]

        Notes
        -----
        Add Braginskii note
        """
        return np.array(
            [self._norm_alpha_para, self._norm_alpha_perp, self._norm_alpha_cross]
        )

    @property
    def _norm_beta_para(self):
        i = self._find_nearest_Z(self.Z)
        return self._c["beta0"][i] * np.ones(self.chi_e.size)

    @property
    def _norm_beta_perp(self):
        i = self._find_nearest_Z(self.Z)
        return (self._c["beta1p"] * self.chi_e ** 2 + self._c["beta0p"]) / self._Delta(
            i
        )

    @property
    def _norm_beta_cross(self):
        i = self._find_nearest_Z(self.Z)
        return (
            self.chi_e
            * (self._c["beta1pp"] * self.chi_e ** 2 + self._c["beta0pp"])
            / self._Delta(i)
        )

    # TODO: Maybe move these collected norm functions also to the base class?
    # then just re-instantiate here with a modified docstring??

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
            from the Braginskii tables:
            Z = [1, 2, 3, 4, ∞]

        Returns
        -------
        norm_beta, `~numpy.ndarray` of `~u.Quantity` instances (3, N)
            The thermoelectric coefficients:
                [beta_para, beta_perp, beta_cross]

        Notes
        -----
        Braginskii Note

        """
        return np.array(
            [self._norm_beta_para, self._norm_beta_perp, self._norm_beta_cross]
        )

    @property
    def _norm_kappa_e_para(self):
        i = self._find_nearest_Z(self.Z)
        return self._c["gamma0"][i] * np.ones(self.chi_e.size)

    @property
    def _norm_kappa_e_perp(self):
        i = self._find_nearest_Z(self.Z)
        return (
            self._c["gamma1p"][i] * self.chi_e ** 2 + self._c["gamma0p"][i]
        ) / self._Delta(i)

    @property
    def _norm_kappa_e_cross(self):
        i = self._find_nearest_Z(self.Z)
        return (
            self.chi_e
            * (self._c["gamma1pp"][i] * self.chi_e ** 2 + self._c["gamma0pp"][i])
            / self._Delta(i)
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
            from the Braginskii tables:
            Z = [1, 2, 3, 4, ∞]

        Returns
        -------
        norm_kapp_e, `~numpy.ndarray` of `~u.Quantity` instances (3, N)
            The electron thermal conductivity coefficients:
                [kappa_e_para, kappa_e_perp, kappa_e_cross]

        """
        return np.array(
            [self._norm_kappa_e_para, self._norm_kappa_e_perp, self._norm_kappa_e_cross]
        )


class BraginskiiInterpolated(AbstractInterpolatedCoefficients):
    @property
    def _data_file(self):
        return os.path.join(data_dir, "braginskii_data.npz")


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    chi = np.linspace(-2, 2, num=50)
    chi = 10 ** chi

    coef = BraginskiiInterpolated(chi_e=chi, Z=5)

    data = coef.norm_alpha[1, :]

    # print(data)

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.plot(chi, data)

    chi = np.linspace(-2, 2, num=100)
    chi = 10 ** chi

    # Instantiate the object
    coef1 = BraginskiiPolynomialFit(chi_e=chi, Z=1)
    coef2 = BraginskiiPolynomialFit(chi_e=chi, Z=np.inf)

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
