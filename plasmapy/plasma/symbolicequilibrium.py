import matplotlib.pyplot as plt
import numpy as np
import plasmaboundaries
import sympy
import warnings

from astropy import constants
from astropy import units as u
from dataclasses import dataclass
from scipy import optimize
from skimage import measure

from plasmapy.plasma.fluxsurface import FluxSurface


@dataclass
class SymbolicEquilibrium:
    aspect_ratio: float
    A: float
    elongation: float
    triangularity: float
    config: str
    B0: float

    def __post_init__(self):
        if self.triangularity > 0.841:
            warnings.warn(
                "You might get a non-convex plasma with a triangularity above 0.841."
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

        psym = (-(psi0 ** 2) / (mu0 * R0 ** 4) * (1 - self.A) * psisym).simplify()
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

    def plot(
        self,
        rminmaxstep=(0.6, 1.4, 0.01),
        zminmaxstep=(-0.6, 0.6, 0.01),
        savepath=None,
        vmax=0,
    ):
        rmin, rmax, rstep = rminmaxstep
        zmin, zmax, zstep = zminmaxstep

        r = np.arange(rmin, rmax, step=rstep)
        z = np.arange(zmin, zmax, step=zstep)
        R, Z = np.meshgrid(r, z)
        PSI = self.psi(R, Z)  # compute magnetic flux

        levels = np.sort(np.linspace(PSI.min(), 0, num=25))
        fig, ax = plt.subplots()
        CS = ax.contourf(R, Z, PSI, levels=levels, vmax=vmax)
        ax.contour(R, Z, PSI, levels=[0], colors="black")  # display the separatrix

        plt.colorbar(CS, label="Magnetic flux $\Psi$")
        ax.set_xlabel("Radius $R/R_0$")
        ax.set_ylabel("Height $z/R_0$")
        ax.set_aspect("equal")
        if savepath is not None:
            ax.savefig(savepath)
        return ax

    def get_flux_surface(
        self, psi_value, rminmaxstep=(0.6, 1.4, 0.01), zminmaxstep=(-0.6, 0.6, 0.01),
    ):
        rmin, rmax, rstep = rminmaxstep
        zmin, zmax, zstep = zminmaxstep

        r = np.arange(rmin, rmax, step=rstep)
        z = np.arange(zmin, zmax, step=zstep)
        R, Z = np.meshgrid(r, z)
        PSI = self.psi(R, Z)  # compute magnetic flux
        contours = measure.find_contours(PSI, psi_value, positive_orientation="high")
        if len(contours) != 1:
            if len(contours) > 1:
                for contour in contours:
                    RcontourArrayUnits, ZcontourArrayUnits = (
                        contour[:, 1],
                        contour[:, 0],
                    )

                    Zcontour = ZcontourArrayUnits / PSI.shape[0] * (zmax - zmin) + zmin
                    Rcontour = RcontourArrayUnits / PSI.shape[1] * (rmax - rmin) + rmin
                    plt.plot(Rcontour, Zcontour)

                rmin, rmax, rstep = rminmaxstep
                zmin, zmax, zstep = zminmaxstep

                r = np.arange(rmin, rmax, step=rstep)
                z = np.arange(zmin, zmax, step=zstep)
                R, Z = np.meshgrid(r, z)
                PSI = self.psi(R, Z)  # compute magnetic flux

                levels = np.linspace(PSI.min(), 0, num=25)
                CS = plt.contourf(R, Z, PSI, levels=levels, vmax=0)
                plt.contour(
                    R, Z, PSI, levels=[0], colors="black"
                )  # display the separatrix

                plt.colorbar(CS, label="Magnetic flux $\Psi$")
                plt.xlabel("Radius $R/R_0$")
                plt.ylabel("Height $z/R_0$")
                plt.gca().set_aspect("equal")
                plt.show()
            raise ValueError(
                f"Could not find contour for psi = {psi_value} (len(contours)={len(contours)})"
            )

        contour = contours[0]
        RcontourArrayUnits, ZcontourArrayUnits = contour[:, 1], contour[:, 0]

        Zcontour = ZcontourArrayUnits / PSI.shape[0] * (zmax - zmin) + zmin
        Rcontour = RcontourArrayUnits / PSI.shape[1] * (rmax - rmin) + rmin

        dZ = np.gradient(Zcontour)
        dR = np.gradient(Rcontour)
        dL = np.sqrt(dZ ** 2 + dR ** 2)

        Brvals = self.Brfunc(Rcontour, Zcontour)
        Bzvals = self.Bzfunc(Rcontour, Zcontour)
        Bprimervals = self.Brdifffunc(Rcontour, Zcontour)
        Bprimezvals = self.Bzdifffunc(Rcontour, Zcontour)
        Bphivals = self.Bphifunc(Rcontour, Zcontour)
        fs = FluxSurface(
            Rcontour,
            Zcontour,
            psi_value,
            Brvals,
            Bzvals,
            Bphivals,
            Bprimervals,
            Bprimezvals,
        )
        return fs


if __name__ == "__main__":
    params = plasmaboundaries.ITER
    eq = SymbolicEquilibrium(**params, config="non-null")
    eq.plot()
    fs = eq.get_flux_surface(0)
    fs.plot(B=True, n=True)
