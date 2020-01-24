from plasmapy.classes.plasma_base import GenericPlasma
import astropy.units as u
import numpy as np

E_unit = u.V / u.m
class Coils(GenericPlasma):
    """
    Work-in-progress class for vacuum magnetic fields due to MagnetoStatics`.

    This is primarily helpful for `plasmapy.simulation.ParticleTracker`.
    """

    def __init__(self, *magnetostatics):
        self.magnetostatics = magnetostatics

    def _interpolate_E(self, r: u.m):
        return np.zeros(r.shape)
    
    def interpolate_E(self, r: u.m):
        return u.Quantity(self._interpolate_E(r.si.value), E_unit)
    
    def _interpolate_B(self, r):
        B = np.zeros(r.shape)
        for ms in self.magnetostatics:
            field = ms._magnetic_field(r, ms.pt, ms.dl, ms.current, ms.w)
            B += field
        return B

    def interpolate_B(self, r: u.m):
        return u.Quantity(self._interpolate_B(r.si.value), u.T)
        
    @classmethod
    def is_datasource_for(cls, **kwargs):
        match = 'interpolate_B' in kwargs.keys()
        return match

    def visualize(self, figure = None):   # coverage: ignore
        from mayavi import mlab
        if figure is None:
            fig = mlab.figure()
        else:
            fig = figure

        for ms in self.magnetostatics:
            ms.visualize(fig)

        return fig
