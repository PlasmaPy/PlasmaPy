

import astropy.constants as const
import astropy.units as u
import numpy as np
import scipy.interpolate as interp
import warnings



def _calculate_plane_points(point, size, bins):
        """
        Given a position vector, create an array of points in the plane that
        contains the point and is normal to its position vector. The number
        of points is determine by bins, while the size in space is determined
        by the size argument.
        """




class LineIntegratedDiagnostic:


    def __init__(self, grid : u.m,
                 parameters : dict,
                 source: u.m,
                 detector: u.m,
                 geometry = 'cartesian',
                 verbose=True,
                 ):

        self.grid = grid
        self.param = parameters
        self.verbose = verbose

        # TODO: auto-detect geometry based on units of source and/or detector

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

        # Calculate normal vectors (facing towards the grid origin) for both
        # the source and detector planes
        self.src_n = self.source.value / np.linalg.norm(
            self.source.value
        )
        self.det_n = -self.detector.value / np.linalg.norm(
            self.detector.value
        )
        # Vector directly from source to detector
        self.src_det_vec = self.detector - self.source

        # Experiment axis is the unit vector from the source to the detector
        self.src_det_n = self.src_det_vec.value / np.linalg.norm(
            self.src_det_vec.value
        )



    def _log(self, msg):
        if self.verbose:
            print(msg)


    # ************************************************************************
    # Functions used by both the integration and particle-tracing method
    # ************************************************************************



    # ************************************************************************
    # Integration Method (linear approximation)
    # ************************************************************************


    def evaluate(self, size=None, bins=None, collimated=True,
                  num=10):
        """
        Analagous to run

        1) Create detector grid from size and bin keywords

        2) For each cell of detector grid, create an array of points of
        separation ds from there to the source point (or, for collimated, the source plane).

        3) Evaluate the integrand at each position in the 3D grid

        4) Integrate along the line-integrated dimension

        """

        # Create 2D grids of detector points
        # Define plane  horizontal axis as being perpendicular to both the
        # position vector and the z-axis. In the case where the pos vector
        # is aligned with the z-axis, just automatically chose the horizontal
        # axis to be the x-axis.
        if np.allclose(np.abs(self.det_n), np.array([0, 0, 1])):
            nx = np.array([1, 0, 0])
        else:
            nx = np.cross(np.array([0, 0, 1]), self.det_n)
        nx = nx / np.linalg.norm(nx)

        # Define the detector vertical axis as being orthogonal to the
        # detector axis and the horizontal axis
        ny = np.cross(nx, self.det_n)
        ny = ny / np.linalg.norm(ny)

        xax = np.linspace(size[0][0], size[0][1], num=bins[0])
        yax = np.linspace(size[1][0], size[1][1], num=bins[1])
        x_offset, y_offset = np.meshgrid(xax, yax, indexing='ij')


        det_pts = np.outer(x_offset,nx) + np.outer(y_offset,ny) + self.detector

        det_pts = np.reshape(det_pts, [bins[0], bins[1],3])


        # Create 2D grids of source points
        if collimated:
            src_pts = det_pts - self.src_det_vec
        else:
            # If not collimated, assume a point source
            src_pts = np.outer( np.ones([bins[0], bins[1]]), self.source)


        # Determine where the grid begins and ends as fractions of the
        # source-to-detector vector
        source_to_det = np.linalg.norm(self.src_det_vec).to(u.mm)
        source_to_grid = np.min(np.linalg.norm(self.grid -
                                               self.source, axis=3))
        grid_to_det = np.min(np.linalg.norm(self.grid -
                                               self.detector, axis=3))

        # Paramter for parameteric equation
        start = source_to_grid/source_to_det
        stop = 1 - grid_to_det/source_to_det
        i = np.linspace(start, stop, num=num)
        ds = (stop-start)*source_to_det/num

        # Create an array of points along lines between points in src_pts
        # and det_pts of the equation m*i + b where
        # m = det_pts - src_pts
        # Integrating along the lines corresponds to summation along ax=2
        mi = np.outer(det_pts - src_pts, i)
        mi = np.reshape(mi, [bins[0], bins[1], 3, num] )

        b = np.outer(src_pts, np.ones(num))
        b = np.reshape(b, [bins[0], bins[1], 3, num])

        pts = mi+b
        pts = np.moveaxis(pts, 2, 3)

        # Evaluate the integrand at each position
        integrand = self._integrand(pts)

        # Integrate
        integral = np.trapz(integrand, axis=2)*ds

        return xax, yax, integral





    def _create_points(self, point, size, bins):
        """
        Create an array of points
        """


        return points


    def _integrand(self, pts):
        """
        Returns the integrand at a particular position.

        This method is over-written by subclasses to reflect actual
        diagnostic physics, and can pull plasma parameter information from
        the parameters dict provided.
        """
        return pts[:,:,:,1]



    # ************************************************************************
    # Particle Tracing Method (for non-linear problems)
    # ************************************************************************




if __name__ == '__main__':

    # Make a little grid
    ax = np.linspace(-1,1,3)*u.mm
    xarr, yarr, zarr = np.meshgrid(ax,ax,ax, indexing='ij')
    grid = np.array([xarr, yarr, zarr])*u.mm
    grid = np.moveaxis(grid, 0, -1)

    parameters = {}
    source = (-3*u.mm, 0*u.mm, 0*u.mm)
    detector = (5*u.mm, 0*u.mm, 0*u.mm)


    obj = LineIntegratedDiagnostic(grid, parameters, source, detector)

    obj.evaluate(size=np.array([[-1,1],[-1,1]])*u.mm, bins=[3,3])




