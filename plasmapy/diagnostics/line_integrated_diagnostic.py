

import astropy.constants as const
import astropy.units as u
import numpy as np
import scipy.interpolate as interp
import warnings


class LineIntegratedDiagnostic:


    def __init__(self, grid : u.m
                 parameters : dict
                 source: u.m,
                 detector: u.m,
                 geometry="cartesian",
                 verbose=True,
                 ):

        self.grid = grid
        self.param = parameters
        self.verbose = verbose


        # Convert geometrical inputs between coordinates systems
        if geometry == "cartesian":
            x, y, z = source
            self.source = np.zeros(3) * u.m
            self.source[0] = x.to(u.m)
            self.source[1] = y.to(u.m)
            self.source[2] = z.to(u.m)

            x, y, z = detector
            self.detector = np.zeros(3) * u.m
            self.detector[0] = x.to(u.m)
            self.detector[1] = y.to(u.m)
            self.detector[2] = z.to(u.m)

        elif geometry == "cylindrical":
            r, t, z = source
            r = r.to(u.m)
            t = t.to(u.rad).value
            z = z.to(u.m)
            self.source = np.zeros(3) * u.m
            self.source[0] = r * np.cos(t)
            self.source[1] = r * np.sin(t)
            self.source[2] = z

            r, t, z = detector
            r = r.to(u.m)
            t = t.to(u.rad).value
            z = z.to(u.m)
            self.detector = np.zeros(3) * u.m
            self.detector[0] = r * np.cos(t)
            self.detector[1] = r * np.sin(t)
            self.detector[2] = z

        elif geometry == "spherical":
            r, t, p = source
            r = r.to(u.m)
            t = t.to(u.rad).value
            p = p.to(u.rad).value
            self.source = np.zeros(3) * u.m
            self.source[0] = r * np.sin(t) * np.cos(p)
            self.source[1] = r * np.sin(t) * np.sin(p)
            self.source[2] = r * np.cos(t)

            r, t, p = detector
            r = r.to(u.m)
            t = t.to(u.rad).value
            p = p.to(u.rad).value
            self.detector = np.zeros(3) * u.m
            self.detector[0] = r * np.sin(t) * np.cos(p)
            self.detector[1] = r * np.sin(t) * np.sin(p)
            self.detector[2] = r * np.cos(t)

        self._log("Source: " + str(self.source.to(u.mm)))
        self._log("Detector: " + str(self.detector.to(u.mm)))


    def _log(self, msg):
        if self.verbose:
            print(msg)


    def evaluate(self, size=None, bins=None, collimated=False,
                  ds=None):
        """
        Analagous to run

        1) Create detector grid from size and bin keywords

        2) For each cell of detector grid, create an array of points of
        separation ds from there to the source point (or, for collimated, the source plane).

        3) Evaluate the integrand at each position in the 3D grid

        4) Integrate along the line-integrated dimension

        """
        pass


    def _integrand(self, pos):
        """
        Returns the integrand at a particular position.

        This method is over-written by subclasses to reflect actual
        diagnostic physics, and can pull plasma parameter information from
        the parameters dict provided.
        """
        pass






