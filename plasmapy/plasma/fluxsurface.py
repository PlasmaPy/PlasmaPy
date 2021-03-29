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

    def __post_init__(self):
        self.Bvectors = np.stack((self.Brvals, self.Bzvals))
        self.dZ = np.gradient(self.Z)
        self.dR = np.gradient(self.R)
        self.dL = np.sqrt(self.dZ ** 2 + self.dR ** 2)

        normal_vectors = np.stack((self.dZ, -self.dR))
        # TODO ej, n według shainga to \vec{B} / |B|. Ale dalej nigdzie go nie używam
        self.normal_vectors = normal_vectors / np.linalg.norm(normal_vectors, axis=0)

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
        """Test for documentation of a variable in post_init"""
        self.Theta = integral

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
            nr, nz = self.normal_vectors
            ax.quiver(self.R, self.Z, nr, nz, **kwargs)
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
