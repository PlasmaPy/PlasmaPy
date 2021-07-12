import matplotlib.pyplot as plt
import numpy as np

from astropy import constants
from astropy import units as u
from dataclasses import dataclass, field

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

    R: np.ndarray = field(repr=False)
    Z: np.ndarray = field(repr=False)
    psi: float
    Brvals: np.ndarray = field(repr=False)
    Bzvals: np.ndarray = field(repr=False)
    Bphivals: np.ndarray = field(repr=False)
    Bprimervals: np.ndarray = field(repr=False)
    Bprimezvals: np.ndarray = field(repr=False)
    GradRho2: np.ndarray = field(repr=False)

    def __post_init__(self):
        self.centroid = np.array([np.mean(self.R), np.mean(self.Z)])
        self.Bvectors = np.stack((self.Brvals, self.Bzvals))
        self.dZ = np.gradient(self.Z)
        self.dR = np.gradient(self.R)
        self.dL = np.sqrt(self.dZ ** 2 + self.dR ** 2)
        self.lp = np.cumsum(self.dL)
        self.Lp = np.sum(self.dL)

        self.toroid_area = self.centroid[0] * 2 * np.pi * self.Lp
        R0 = 1
        self.Cp = self.Lp / R0

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
        self, ax=None, *, quantity = None, n=False, B=False, legend=True, colorbar = True, **kwargs
    ):  # coverage: ignore
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect("equal")
            ax.set_xlabel("Radius $R/R_0$")
            ax.set_ylabel("Height $z/R_0$")

        if quantity is None:
            ax.plot(self.R, self.Z, label=fr"$\psi = ${self.psi:.2f}")
        else:
            from matplotlib.collections import LineCollection
            points = np.array([self.R, self.Z]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            # Create a continuous norm to map from data points to colors
            norm = plt.Normalize(quantity.min(), quantity.max())
            lc = LineCollection(segments, cmap='plasma', norm=norm)
            # Set the values used for colormapping
            lc.set_array(quantity)
            lc.set_linewidth(2)
            line = ax.add_collection(lc)
            if colorbar:
                ax.get_figure().colorbar(line, ax=ax)

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
        # lin-liu1995, equation 7
        h, h2mean = self._h, self._h2mean
        innermost_left = np.sqrt(1 - h)
        innermost_right = 1 + 1 / 2 * h
        square_bracket = 1 - innermost_left * innermost_right
        integrand = h ** -2 * square_bracket
        f_tl = 1 - h2mean * self.flux_surface_average(integrand)
        return f_tl

    @cached_property
    def trapped_fraction(self):
        f_t = 0.75 * self._f_tu + 0.25 * self._f_tl
        return f_t

    def trapped_fraction_analytical(self, N_integration_points=1000):
        h = self._h
        Lambda = np.linspace(0, 1, N_integration_points)
        integrand = np.sqrt(1 - Lambda[:, np.newaxis] * h[np.newaxis, :])
        FSA = self.flux_surface_average(integrand)
        integrand2 = Lambda / FSA
        the_integral = trapezoid(integrand2, Lambda)
        estimated_trapped_fraction = 1 - 3 / 4 * self._h2mean * the_integral
        return estimated_trapped_fraction

    @cached_property
    def sqrt_jacobian(self):
        # as per http://fusionwiki.ciemat.es/fusionwiki/index.php?title=Flux_coordinates&oldid=4501
        # TODO double check with PR08
        sqrt_jacobian = self.toroid_area / 4 / np.pi ** 2
        return sqrt_jacobian

    @cached_property
    def BDotNablaThetaFSA(self):
        # Eq B14
        BdotNablaTheta = self.psiprime / (2 * np.pi * self.sqrt_jacobian)
        return self.flux_surface_average(BdotNablaTheta)

    def F_m(self, M: int = 3):
        """Mode weights for the Pfirsch-Schlüter contribution.

        Equation B9 from |Houlberg_1997|.

        Parameters
        ----------
        m : Union[int, np.ndarray]
            m
        flux_surface : FluxSurface
            flux_surface
        """
        m = np.arange(1, M + 1).reshape(-1, 1)
        Theta = self.Theta.reshape(1, -1)
        B20 = (self.Brvals * self.Bprimervals + self.Bzvals * self.Bprimezvals).reshape(1, -1)
        under_average_B16 = np.sin(Theta * m) * B20
        under_average_B15 = under_average_B16 / self.Bmag.reshape(1, -1)
        under_average_B16_cos = np.cos(Theta * m) * B20
        under_average_B15_cos = under_average_B16_cos / self.Bmag.reshape(1, -1)
        B15 = self.flux_surface_average(under_average_B15)
        B16 = self.gamma * self.flux_surface_average(under_average_B16)
        B15_cos = self.flux_surface_average(under_average_B15_cos)
        B16_cos = self.gamma * self.flux_surface_average(under_average_B16_cos)

        B2mean = self.fsa_B2 / u.T**2

        F_m = 2 / B2mean / self.BDotNablaThetaFSA * (B15 * B16 + B15_cos * B16_cos)
        return F_m

    @cached_property
    def fsa_B2(self):
        # flux surface averaged B^2
        return self.flux_surface_average(self.B2) * u.T ** 2

    @cached_property
    def fsa_invB2(self):
        return self.flux_surface_average(1 / self.B2) / u.T ** 2

    @cached_property
    def grbm2(self):
        return self.flux_surface_average(self.GradRho2 / self.B2) / u.T ** 2

    @cached_property
    def _B17(self):
        """Equation B17 from |Houlberg_1997|. Likely bugged!

        Notes
        -----
        Eventually this should allow picking the right `m` in `K` below.

        Parameters
        ----------
        flux_surface :
            flux_surface
        """
        B20 = self.Brvals * self.Bprimervals + self.Bzvals * self.Bprimezvals
        under_average_B17 = (B20 / self.Bmag) ** 2
        return self.flux_surface_average(under_average_B17) / self.flux_surface_average(self.B2)



    @cached_property
    def F_m3(self):
        return self.F_m()
