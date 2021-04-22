
import numpy as np
import astropy.units as u

from plasmapy.diagnostics.line_integrated_diagnostic import LineIntegratedDiagnostic
from plasmapy.plasma.grids import CartesianGrid

import matplotlib.pyplot as plt

class DensityIntegrator(LineIntegratedDiagnostic):
    def _integrand(self, pts):
        return self.grid['rho'][pts]




def test_scalar_integrator():
    # Make a little grid
    ax = np.linspace(-1,1, num=50)*u.mm
    xarr, yarr, zarr = np.meshgrid(ax,ax,ax, indexing='ij')
    
    radius = np.sqrt(xarr**2 + yarr**2)
    
    field = np.where(radius < 0.5*u.mm, (1/radius.value +0.05)*u.kg / u.m ** 3, 0*u.kg / u.m ** 3)
    
    
    grid = CartesianGrid(xarr, yarr, zarr)
    grid.add_quantities(rho=field)
    
    print(grid.pts0.unit)


    source = (0*u.mm,  -3*u.mm, 0*u.mm)
    detector = ( 0*u.mm, 5*u.mm,  0*u.mm)


    obj = DensityIntegrator(grid, source, detector)

    size = np.array([[-1,1],[-1,1]])*u.mm
    hax, vax, integral, _ = obj.evaluate_integral(size=size, bins=[20,20],
                          collimated=True, num=10)


if __name__ == '__main__':
    test_scalar_integrator()