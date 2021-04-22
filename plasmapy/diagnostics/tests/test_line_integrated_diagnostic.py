
import numpy as np
import astropy.units as u

from plasmapy.diagnostics.line_integrated_diagnostic import LineIntegrateScalarQuantities
from plasmapy.plasma.grids import CartesianGrid

import matplotlib.pyplot as plt

"""
class DensityIntegrator(LineIntegratedDiagnostic):

    def integrand(self, pts):
        rho = self.grid.volume_averaged_interpolator(pts, 'rho')
        return rho
"""



def test_scalar_integrator():
    # Make a little grid
    ax = np.linspace(-1,1, num=100)*u.mm
    xarr, yarr, zarr = np.meshgrid(ax,ax,ax, indexing='ij')
    
    radius = np.sqrt(xarr**2 + yarr**2)
    
    field = np.where(radius < 0.5*u.mm, 1, 0) *  1*u.kg / u.m ** 3

    grid = CartesianGrid(xarr, yarr, zarr)
    grid.add_quantities(rho=field)

    source = (0*u.mm,  -5*u.mm, 0*u.mm)
    detector = ( 0*u.mm, 5*u.mm,  0*u.mm)


    obj = LineIntegrateScalarQuantities(grid, source, detector, 'rho')
    
    print(obj.source)

    size = np.array([[-2,2],[-2,2]])*u.mm
    hax, vax, integral = obj.line_integral(size=size, bins=[100,100],
                          collimated=False, num=50)
    

    plt.pcolormesh(hax.to(u.mm).value, vax.to(u.mm).value, 
                   integral.T, shading='auto')
    plt.show()


if __name__ == '__main__':
    test_scalar_integrator()