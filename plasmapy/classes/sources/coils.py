from typing import Iterable

import astropy.units as u
import numpy as np

from plasmapy.classes.plasma_base import GenericPlasma
from plasmapy.formulary import magnetostatics

E_unit = u.V / u.m


class Coils(GenericPlasma):
    def __init__(self, *magnetostatics: magnetostatics.CircularWire):
        """
        Work-in-progress class for vacuum magnetic fields due to MagnetoStatics`.

        This is primarily helpful for `plasmapy.simulation.ParticleTracker`.

        Note: currently this accepts only circular coils as magnetostatics.

        Parameters
        ----------
        magnetostatics : plasmapy.formulary.magnetostatics.CircularWire
            any number of CircularWire
        """
        self.magnetostatics = magnetostatics

    @staticmethod
    def _interpolate_E(r: np.ndarray):
        return np.zeros(r.shape)

    def interpolate_E(self, r: u.m):
        """
        Returns the electric field, with units.

        Note that MagnetoStatics provide no electric field, so this is effectively zero.

        Parameters
        ----------
        r : u.m
            Position; can be a (3,) array or a (N, 3) array.

        Returns
        -------
        electric field array of the same shape as the input.
        """
        return u.Quantity(self._interpolate_E(r.si.value), E_unit)

    def _interpolate_B(self, r: np.ndarray):
        B = np.zeros(r.shape, dtype=float)
        for ms in self.magnetostatics:
            field = ms._magnetic_field(r, ms.pt, ms.dl, ms.current, ms.w)
            B += field
        return B

    def interpolate_B(self, r: u.m):
        """
        Returns the magnetic field from all the MagnetoStatics, with units.

        Parameters
        ----------
        r : u.m
            Position; can be a (3,) array or a (N, 3) array.

        Returns
        -------
        Magnetic field array of the same shape as the input.
        """
        return u.Quantity(self._interpolate_B(r.si.value), u.T, copy=False)

    def __repr__(self):
        names = []
        for ms in self.magnetostatics:
            name = f"{ms.__class__.__name__}({ms.current * u.A})"
            names.append(name)
        from collections import Counter

        items = []
        for key, value in Counter(names).items():
            items.append(f"{value} x {key}")
        return f"{self.__class__.__name__}({', '.join(items)})"

    @classmethod
    def is_datasource_for(cls, **kwargs):
        match = "interpolate_B" in kwargs.keys()
        return match

    def visualize(self, figure=None):  # coverage: ignore
        """
        Visualizes the set of coils in 3D using PyVista.

        Parameters
        ----------
        figure : pyvista.BasePlotter, optional
            If not provided, a new plotter is created.

        Returns
        -------
        the same `pyvista.BasePlotter` instance as (if) provided;
        this allows method chaining.
        """
        import pyvista as pv

        if figure is None:
            fig = pv.Plotter()
        else:
            fig = figure

        for ms in self.magnetostatics:
            ms.visualize(fig)

        return fig

    @classmethod
    def toykamak(
        cls,
        minor_radius: u.m = 0.3 * u.m,
        radius: u.m = 1 * u.m,
        main_current: u.A = 15 * u.MA,
        coil_currents=None,
    ):
        """
        Creates a set of coils for the Toykamak model.

        Parameters
        ----------
        minor_radius : u.m
            radius of the poloidal coils.
        radius : u.m
            major radius of the single coil representing the plasma
        main_current : u.A
            current through the single plasma-representing coil
        coil_currents :  Iterable[u.A]
            list of currents through the coils. By default, 8 15-MA coils are created.

        Returns
        -------
        instance of `Coils`
        """
        if coil_currents is None:
            coil_currents = 8 * [10 * u.MA]
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
