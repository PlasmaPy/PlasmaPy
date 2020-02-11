from plasmapy.classes.plasma_base import GenericPlasma
import astropy.units as u
import numpy as np
from plasmapy.formulary import magnetostatics

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
        match = "interpolate_B" in kwargs.keys()
        return match

    def visualize(self, figure=None):  # coverage: ignore
        import pyvista as pv

        if figure is None:
            fig = pv.Plotter(notebook=True)
        else:
            fig = figure

        for ms in self.magnetostatics:
            ms.visualize(fig)

        if figure is None:
            fig.show()

        return fig

    @classmethod
    def toykamak(
        cls,
        minor_radius=0.3 * u.m,
        radius=1 * u.m,
        main_current=15 * u.MA,
        coil_currents=8 * [10 * u.MA],
    ):
        """
        Creates a set of coils for the Toykamak model.
        """
        n_coils = len(coil_currents)
        currents = u.Quantity(coil_currents)

        coil_angles = np.linspace(0, 2 * np.pi, n_coils, endpoint=False)

        coils = []
        for coil_angle, current in zip(coil_angles, currents):
            x = radius * np.cos(coil_angle)
            y = radius * np.sin(coil_angle)
            normal_angle = np.pi / 2 + coil_angle
            normal = u.Quantity([np.cos(normal_angle), np.sin(normal_angle), 0])
            center = u.Quantity([x, y, 0 * u.m])
            coil = magnetostatics.CircularWire(normal, center, minor_radius, current)
            coils.append(coil)

        if main_current.si.value != 0.0:
            plasma_wire = magnetostatics.CircularWire(
                [0, 0, 1], u.Quantity((0, 0, 0), u.m), radius, main_current
            )
            coils.append(plasma_wire)

        c = cls(*coils)
        return c
