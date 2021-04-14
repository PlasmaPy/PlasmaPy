import matplotlib.pyplot as plt
import numpy as np

from astropy import constants
from dataclasses import dataclass
from scipy import integrate
from typing import Callable

try:
    from functools import cached_property
except ImportError:
    from cached_property import cached_property

try:
    from scipy.integrate import cumtrapz as cumulative_trapezoid
except ImportError:
    from scipy.integrate import cumulative_trapezoid


@dataclass
class FluxSurface:
    R: np.ndarray
    Z: np.ndarray
    psi: float
    Brvals: np.ndarray
    Bzvals: np.ndarray
    Bphivals: np.ndarray
    Bprimervals: np.ndarray
    Bprimezvals: np.ndarray
    GradRho: float = None

    def __post_init__(self):
        self.Bvectors = np.stack((self.Brvals, self.Bzvals))
        self.dZ = np.gradient(self.Z)
        self.dR = np.gradient(self.R)
        self.dL = np.sqrt(self.dZ ** 2 + self.dR ** 2)

        n = np.stack((self.dZ, -self.dR))
        # TODO ej, n według shainga to \vec{B} / |B|. Ale dalej nigdzie go nie używam
        self.n = n / np.linalg.norm(n, axis=0)
        self.nr, self.nz = self.n

        self.lp = np.cumsum(self.dL)
        self.Lp = np.sum(self.dL)
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
        psiprime = 1  # TODO if rho \equiv psi, this should work?
        self.Fhat = constants.mu0 * F

    def plot(self, ax=None, n=False, B=False, legend=True, **kwargs):
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

    def plot_psi(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_xlabel("L")
            ax.set_ylabel(r"$\psi$")
        ax.plot(self.lp, self.psi(self.R, self.Z))
        ax.axhline(self.psi)
        return ax

    def flux_surface_average(self, quantity, simpson=False):
        # TODO test the <B dot ... > = 0 property
        integrand = 1 / self.Bp
        if simpson:
            integrator = integrate.simpson
        else:
            try:
                from scipy.integrate import trapz as trapezoid
            except ImportError:
                from scipy.integrate import trapezoid
            integrator = trapezoid
        return integrator(integrand * quantity, self.lp) / integrator(
            integrand, self.lp
        )

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
        h, hmean, h2mean = self._h, self._hmean, self._h2mean
        f_tl = 1 - h2mean * self.flux_surface_average(
            h ** -2 * (1 - (1 - h) ** 0.5) * (1 + h / 2)
        )
        return f_tl

    def trapped_fraction(self):
        f_t = 0.75 * self._f_tu + 0.25 * self._f_tl
        return f_t

    @cached_property
    def BDotNablaThetaFSA(fs):
        dthetadR = np.gradient(fs.theta, fs.R)
        dthetadZ = np.gradient(fs.theta, fs.Z)
        dot_product = fs.Brvals * dthetadR + fs.Bzvals * dthetadZ
        return fs.flux_surface_average(dot_product)
