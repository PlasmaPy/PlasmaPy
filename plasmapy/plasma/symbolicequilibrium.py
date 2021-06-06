import functools
import matplotlib.pyplot as plt
import numpy as np
import plasmaboundaries
import sympy
import warnings

from astropy import constants
from astropy import units as u
from collections import namedtuple
from dataclasses import dataclass
from scipy import interpolate, optimize
from skimage import measure

from plasmapy.plasma.fluxsurface import FluxSurface

grid_and_psi = namedtuple("GridAndPsi", ["R", "Z", "psi"])


@dataclass
class SymbolicEquilibrium:
    """
    Analytical solution of the GS equation (Cerfon, Freidberg),
    based on plasmaboundaries package
    """

    aspect_ratio: float
    A: float
    elongation: float
    triangularity: float
    config: str
    B0: float

    def __post_init__(self):
        if self.triangularity > 0.841:
            warnings.warn(
                f"As per Cerfon (2010), your plasma surfaces may not be convex with a triangularity of {self.triangularity} > 0.841."
            )
        params = dict(
            aspect_ratio=self.aspect_ratio,
            A=self.A,
            elongation=self.elongation,
            triangularity=self.triangularity,
        )

        self.psi = plasmaboundaries.compute_psi(params, config=self.config)
        self.symbols = Rsym, Zsym = sympy.symbols("R Z")
        self.psisym = psisym = self.psi(Rsym, Zsym, pkg="sp")

        psifunc = sympy.lambdify([(Rsym, Zsym)], psisym, modules="numpy")
        minimization = optimize.minimize(psifunc, x0=[1.0, 0.0])
        R0, Z0 = minimization.x
        psi0 = minimization.fun
        Br = self.Br = -psisym.diff(Zsym) / Rsym
        Bz = self.Bz = psisym.diff(Rsym) / Rsym
        B = sympy.sqrt(Br ** 2 + Bz ** 2)
        mu0 = constants.mu0.si.value
        Bdiff_r = B.diff(Rsym)
        Bdiff_z = B.diff(Zsym)

        self.psym = (-(psi0 ** 2) / (mu0 * R0 ** 4) * (1 - self.A) * psisym).simplify()
        self.Bphi2 = (
            R0 ** 2
            / Rsym ** 2
            * (self.B0 ** 2 - 2 * psi0 ** 2 / R0 ** 4 * self.A * psisym)
        )
        self.Bphi = self.Bphi2 ** 0.5
        self.Bphifunc = sympy.lambdify((Rsym, Zsym), self.Bphi)
        self.Brfunc = sympy.lambdify((Rsym, Zsym), Br)
        self.Bzfunc = sympy.lambdify((Rsym, Zsym), Bz)
        self.Brdifffunc = sympy.lambdify((Rsym, Zsym), Bdiff_r)
        self.Bzdifffunc = sympy.lambdify((Rsym, Zsym), Bdiff_z)
        # assert (
        #     (Br * Rsym).diff(Rsym) / Rsym + Bz.diff(Zsym)
        # ).simplify() == 0  # due to toroidal symmetry
        # TODO change to close to 0 evaluated on grid

    # @functools.lru_cache # TODO get this to work somehow
    def get_grid_and_psi(self, rminmaxstep, zminmaxstep):
        rmin, rmax, rstep = rminmaxstep
        zmin, zmax, zstep = zminmaxstep

        r = np.arange(rmin, rmax, step=rstep)
        z = np.arange(zmin, zmax, step=zstep)
        R, Z = np.meshgrid(r, z)
        PSI = self.psi(R, Z)  # compute magnetic flux
        return grid_and_psi(R, Z, PSI)

    def plot(
        self,
        rminmaxstep=(0.6, 1.4, 0.01),
        zminmaxstep=(-0.6, 0.6, 0.01),
        savepath=None,
        vmax=0,
    ):  # coverage: ignore
        R, Z, PSI = self.get_grid_and_psi(rminmaxstep, zminmaxstep)

        levels = np.sort(np.linspace(PSI.min(), 0, num=25))
        fig, ax = plt.subplots()
        CS = ax.contourf(R, Z, PSI, levels=levels, vmax=vmax)
        ax.contour(R, Z, PSI, levels=[0], colors="black")  # display the separatrix

        plt.colorbar(CS, label=r"Magnetic flux $\Psi$")
        ax.set_xlabel("Radius $R/R_0$")
        ax.set_ylabel("Height $z/R_0$")
        ax.set_aspect("equal")
        if savepath is not None:
            ax.savefig(savepath)
        return ax

    def get_flux_surface(
        self,
        psi_value,
        *,
        rminmaxstep=(0.6, 1.4, 0.01),
        zminmaxstep=(-0.6, 0.6, 0.01),
        RZPSI=None,
    ):
        if RZPSI is not None:
            R, Z, PSI = RZPSI
            rmax = R.max()
            rmin = R.min()
            zmax = Z.max()
            zmin = Z.min()
        else:
            rmin, rmax, rstep = rminmaxstep
            zmin, zmax, zstep = zminmaxstep
            R, Z, PSI = self.get_grid_and_psi(rminmaxstep, zminmaxstep)

        contours = measure.find_contours(PSI, psi_value, positive_orientation="high")
        if len(contours) == 0:
            raise ValueError(f"Could not find contour for psi = {psi_value}")
        elif len(contours) > 1:
            # as per skimage docs:
            # (The closed-ness of a contours can be tested by checking whether the beginning point is the same as the end point.)
            for contour in contours:
                contour_is_closed = (contour[0] == contour[-1]).all()
                if contour_is_closed:
                    break
            else:
                raise ValueError(f"Could not find closed contour for psi = {psi_value}")
        else:
            contour = contours[0]
        RcontourArrayUnits, ZcontourArrayUnits = contour[:, 1], contour[:, 0]

        Zcontour = ZcontourArrayUnits / PSI.shape[0] * (zmax - zmin) + zmin
        Rcontour = RcontourArrayUnits / PSI.shape[1] * (rmax - rmin) + rmin

        dZ = np.gradient(Zcontour)
        dR = np.gradient(Rcontour)
        self.dL = np.sqrt(dZ ** 2 + dR ** 2)

        Brvals = self.Brfunc(Rcontour, Zcontour)
        Bzvals = self.Bzfunc(Rcontour, Zcontour)
        Bprimervals = self.Brdifffunc(Rcontour, Zcontour)
        Bprimezvals = self.Bzdifffunc(Rcontour, Zcontour)
        Bphivals = self.Bphifunc(Rcontour, Zcontour)
        ρ = PSI / PSI.min()
        ρprime_r = np.gradient(ρ, R[0], axis=1)
        ρprime_z = np.gradient(ρ, Z[:, 0], axis=0)

        ρprime2 = ρprime_z ** 2 + ρprime_r ** 2
        interpolator = interpolate.RectBivariateSpline(Z[:, 0], R[0], ρprime2)

        interpolated_GradRho2 = interpolator(Zcontour, Rcontour, grid=False)
        fs = FluxSurface(
            Rcontour,
            Zcontour,
            psi_value,
            Brvals,
            Bzvals,
            Bphivals,
            Bprimervals,
            Bprimezvals,
            interpolated_GradRho2,
        )
        return fs

    def rho_to_psi(
        self,
        rho,
        rminmaxstep=(0.6, 1.4, 0.01),
        zminmaxstep=(-0.6, 0.6, 0.01),
    ):
        R, Z, PSI = self.get_grid_and_psi(rminmaxstep, zminmaxstep)
        psi_max = 0
        psi_min = PSI.min()
        psi_span = psi_max - psi_min

        psi_values = np.asarray(rho) * psi_span + psi_min
        return psi_values

    # TODO def psi_to_rho(self,
    #                psi

    def get_multiple_flux_surfaces(
        self,
        *,
        psi_values=None,
        rho_values=None,
        rminmaxstep=(0.6, 1.4, 0.01),
        zminmaxstep=(-0.6, 0.6, 0.01),
        permissive=True,
    ):
        if psi_values is None and rho_values is None:
            raise ValueError("psi_values and rho_values cannot both be None")
        elif psi_values is not None and rho_values is not None:
            raise ValueError("psi_values and rho_values cannot both be specified")
        rmin, rmax, rstep = rminmaxstep
        zmin, zmax, zstep = zminmaxstep
        R, Z, PSI = self.get_grid_and_psi(rminmaxstep, zminmaxstep)

        if rho_values is not None:
            psi_values = self.rho_to_psi(rho_values)

        for psi in psi_values:
            try:
                surface = self.get_flux_surface(
                    psi,
                    RZPSI=(R, Z, PSI),
                )
            except ValueError as e:
                if permissive:
                    print(psi, e)
                    yield (psi, None)
                else:
                    raise

            yield (psi, surface)


if __name__ == "__main__":
    params = plasmaboundaries.ITER
    eq = SymbolicEquilibrium(**params, config="non-null")
    eq.plot()
    fs = eq.get_flux_surface(0)
    fs.plot(B=True, n=True)
