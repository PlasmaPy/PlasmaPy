"""
Defines SyntheticProtonRadiograph, a subclass of the LineIntegratedDiagnostic class
"""

__all__ = [
        "SyntheticProtonRadiograph",
        ]

import numpy as np
import astropy.units as u
import astropy.constants as constants

from plasmapy.diagnostics.line_integrated_diagnostic import LineIntegratedDiagnostic


class SyntheticProtonRadiograph(LineIntegratedDiagnostic):

    def integrand(self, pts):

        nx,ny,nz,ndim = pts.shape

        # Interpolate field quantites of interest
        pts = np.reshape(pts, (nx*ny*nz, ndim))
        Ex, Ey, Ez, Bx, By, Bz= self.grid.volume_averaged_interpolator(pts, 'E_x', 'E_y', 'E_z','B_x', 'B_y', 'B_z' )

        # Create E and B arrays
        E = np.zeros([nx,ny,nz,3])*u.V/u.m
        E[...,0] = np.reshape(Ex, (nx,ny,nz))
        E[...,1] = np.reshape(Ey, (nx,ny,nz))
        E[...,2] = np.reshape(Ez, (nx,ny,nz))

        B = np.zeros([nx,ny,nz,3])*u.T
        B[...,0] = np.reshape(Bx, (nx,ny,nz))
        B[...,1] = np.reshape(By, (nx,ny,nz))
        B[...,2] = np.reshape(Bz, (nx,ny,nz))


        # Compute components of E and B in the horizontal and
        # vertical directions (relative to the detector plane coordinates)
        Ehax = np.dot(E, self.det_hax)
        Ehax = np.reshape(Ehax, (nx,ny,nz))

        Evax = np.dot(E, self.det_vax)
        Evax  = np.reshape(Evax , (nx,ny,nz))

        Bhax = np.dot(B, self.det_hax)
        Bhax  = np.reshape(Bhax , (nx,ny,nz))

        Bvax = np.dot(B, self.det_vax)
        Bvax  = np.reshape(Bvax , (nx,ny,nz))

        return Ehax, Evax, Bhax, Bvax




    def linear_radiograph(self,size=None, bins=None, collimated=None, num=None,
                          proton_energy = 15*u.MeV):



        hax, vax, Ehax, Evax, Bhax, Bvax = self.line_integral(size=size,
                                                                  bins=bins, collimated=collimated,
                                                                  num=num)

        c1 = constants.e.si/(2*proton_energy)
        c2 = -constants.e.si/np.sqrt(2*constants.m_p.si*proton_energy)

        ahax = (c1*Ehax + c2*Bhax).to(u.dimensionless_unscaled)
        avax = (c1*Evax - c2*Bvax).to(u.dimensionless_unscaled)


        # Apply Eq. 7 from Kugland2012 to calculate the Jacobian
        dh0 = np.mean(np.gradient(hax))
        dv0 = np.mean(np.gradient(vax))

        l1 = np.linalg.norm(self.source)*u.m
        l2 = np.linalg.norm(self.detector)*u.m


        dhdh0 = np.gradient(ahax, axis=0)/dh0
        dvdv0 = np.gradient(avax, axis=1)/dv0
        dhdv0 = np.gradient(ahax, axis=1)/dv0
        dvdh0= np.gradient(avax, axis=0)/dh0

        t1 = self.mag + dhdh0*l2
        t2 = self.mag + dvdv0*l2

        t3 = l2**2*dhdv0*dvdh0
        det = np.abs(t1*t2 - t3)

        # Apply Eq. 6 from Kugland2012 to calculate the intensity

        # Small parameter
        epsilon = 1e-20
        intensity = (1+epsilon)/(det + epsilon) # Possibly *M**2 depending how you define I0?

        return hax, vax, intensity





