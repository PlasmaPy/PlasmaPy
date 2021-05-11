import matplotlib.pyplot as plt
import numpy as np

from astropy import constants
from dataclasses import dataclass

try:
    from functools import cached_property
except ImportError:
    from cached_property import cached_property

try:
    from scipy.integrate import cumtrapz as cumulative_trapezoid
    from scipy.integrate import trapz as trapezoid
except ImportError:
    from scipy.integrate import cumulative_trapezoid, trapezoid


@dataclass
class FluxSurface:
    """Represents a single flux surface out of a magnetic equilibrium, experimental or otherwise."""

    R: np.ndarray
    Z: np.ndarray
    psi: float
    Brvals: np.ndarray
    Bzvals: np.ndarray
    Bphivals: np.ndarray
    Bprimervals: np.ndarray
    Bprimezvals: np.ndarray
    GradRho2: np.ndarray

    def __post_init__(self):
        self.centroid = np.array([np.mean(self.R), np.mean(self.Z)])
        self.Bvectors = np.stack((self.Brvals, self.Bzvals))
        self.dZ = np.gradient(self.Z)
        self.dR = np.gradient(self.R)
        self.dL = np.sqrt(self.dZ ** 2 + self.dR ** 2)
        self.lp = np.cumsum(self.dL)
        self.Lp = np.sum(self.dL)

        self.toroid_area = self.centroid[0] * 2 * np.pi * self.Lp

        n = np.stack((self.dZ, -self.dR))
        # TODO ej, n według shainga to \vec{B} / |B|. Ale dalej nigdzie go nie używam
        self.n = n / np.linalg.norm(n, axis=0)
        self.nr, self.nz = self.n

        self.Bp = np.linalg.norm(self.Bvectors, axis=0)
        self.B2 = self.Bp ** 2 + self.Bphivals ** 2
        self.Bmag = np.sqrt(self.B2)
        self.Bmax = self.Bmag.max()
        self.Bmin = self.Bmag.min()

        integrand = self.Bmag / self.Bp
        integral = cumulative_trapezoid(integrand, self.lp, initial=0)
        self.gamma = 2 * np.pi / integral[-1]
        integral *= self.gamma
        self.Theta = integral

        F = (
            self.flux_surface_average(self.R * self.Bphivals)
            * 2
            * np.pi
            / constants.mu0
        )  # FS average because this is supposed to be constant on the flux surface
        self.psiprime = trapezoid(self.Bphivals, self.lp)
        self.Fhat = constants.mu0 * F

    def plot(
        self, ax=None, n=False, B=False, legend=True, **kwargs
    ):  # coverage: ignore
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect("equal")
            ax.set_xlabel("Radius $R/R_0$")
            ax.set_ylabel("Height $z/R_0$")
        ax.plot(self.R, self.Z, label=fr"$\psi = ${self.psi:.2f}")
        if B:
            ax.quiver(self.R, self.Z, self.Brvals, self.Bzvals, **kwargs)
        if n:
            ax.quiver(self.R, self.Z, self.nr, self.nz, **kwargs)
        if legend:
            ax.legend()
        return ax

    @cached_property
    def _fsa_denominator(self):
        return trapezoid(1 / self.Bp, self.lp)

    def flux_surface_average(self, quantity):
        return trapezoid(quantity / self.Bp, self.lp) / self._fsa_denominator

    @cached_property
    def theta(self):
        theta = 2 * np.pi * self.dL / self.Lp
        return theta

    @cached_property
    def _h(self):
        h = self.Bmag / self.Bmax
        return h

    @cached_property
    def _hmean(self):
        hmean = self.flux_surface_average(self._h)
        return hmean

    @cached_property
    def _h2mean(self):
        h2mean = self.flux_surface_average(self._h ** 2)
        return h2mean

    @cached_property
    def _f_tu(self):
        hmean, h2mean = self._hmean, self._h2mean
        f_tu = 1 - h2mean / hmean ** 2 * (1 - (1 - hmean) ** 0.5 * (1 + hmean / 2))
        return f_tu

    @cached_property
    def _f_tl(self):
        # Houlberg_1997, equations B5-B7
        h, h2mean = self._h, self._h2mean
        f_tl = 1 - h2mean * self.flux_surface_average(
            h ** -2 * (1 - (1 - h) ** 0.5) * (1 + h / 2)
        )
        return f_tl

    def trapped_fraction(self):
        f_t = 0.75 * self._f_tu + 0.25 * self._f_tl
        return f_t

    @cached_property
    def sqrt_jacobian(self):
        sqrt_jacobian = self.toroid_area / 4 / np.pi**2
        return sqrt_jacobian

    @cached_property
    def BDotNablaThetaFSA(self):
        # Eq B14
        BdotNablaTheta = self.psiprime / (2 * np.pi * self.sqrt_jacobian)
        return self.flux_surface_average(BdotNablaTheta)
