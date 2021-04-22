"""
Defines the LineIntegratedDiagnostic base class
"""

__all__ = [
        "LineIntegratedDiagnostic",
        ]

import astropy.units as u
import numpy as np


def _rot_a_to_b(a, b):
    r"""
    Calculates the 3D rotation matrix that will rotate vector a to be aligned
    with vector b.
    """
    # Normalize both vectors
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)

    # Manually handle the case where a and b point in opposite directions
    if np.dot(a, b) == -1:
        return -np.identity(3)

    axb = np.cross(a, b)
    c = np.dot(a, b)
    vskew = np.array(
        [[0, -axb[2], axb[1]], [axb[2], 0, -axb[0]], [-axb[1], axb[0], 0]]
    ).T  # Transpose to get right orientation

    return np.identity(3) + vskew + np.dot(vskew, vskew) / (1 + c)



class LineIntegratedDiagnostic:


    def __init__(self, grid : u.m,
                 source,
                 detector,
                 verbose=True,
                 ):

        # self.grid is the grid object
        self.grid = grid
        # self.grid_arr is the grid positions in si. This is created here
        # so that it isn't continously called later
        self.grid_arr = grid.grid

        self.verbose = verbose

        # Auto-detect source and detector geometry based on units
        geo_units = [x.unit for x in source]
        if geo_units[2].is_equivalent(u.rad):
            geometry = 'spherical'
        elif geo_units[1].is_equivalent(u.rad):
            geometry = 'cylindrical'
        else:
            geometry = 'cartesian'

        # Convert geometrical inputs between coordinates systems
        if geometry == "cartesian":
            x, y, z = source

            self.source = np.zeros(3)
            self.source[0] = x.to(u.m).value
            self.source[1] = y.to(u.m).value
            self.source[2] = z.to(u.m).value

            x, y, z = detector
            self.detector = np.zeros(3)
            self.detector[0] = x.to(u.m).value
            self.detector[1] = y.to(u.m).value
            self.detector[2] = z.to(u.m).value

        elif geometry == "cylindrical":
            r, t, z = source
            r = r.to(u.m)
            t = t.to(u.rad).value
            z = z.to(u.m)
            self.source = np.zeros(3)
            self.source[0] = (r * np.cos(t)).to(u.m).value
            self.source[1] = (r * np.sin(t)).to(u.m).value
            self.source[2] = z.to(u.m).value

            r, t, z = detector
            r = r.to(u.m)
            t = t.to(u.rad).value
            z = z.to(u.m)
            self.detector = np.zeros(3)
            self.detector[0] = (r * np.cos(t)).to(u.m).value
            self.detector[1] = (r * np.sin(t)).to(u.m).value
            self.detector[2] = z.to(u.m).value

        elif geometry == "spherical":
            r, t, p = source
            r = r.to(u.m)
            t = t.to(u.rad).value
            p = p.to(u.rad).value
            self.source = np.zeros(3)
            self.source[0] = (r * np.sin(t) * np.cos(p)).to(u.m).value
            self.source[1] = (r * np.sin(t) * np.sin(p)).to(u.m).value
            self.source[2] = (r * np.cos(t)).to(u.m).value

            r, t, p = detector
            r = r.to(u.m)
            t = t.to(u.rad).value
            p = p.to(u.rad).value
            self.detector = np.zeros(3)
            self.detector[0] = (r * np.sin(t) * np.cos(p)).to(u.m).value
            self.detector[1] = (r * np.sin(t) * np.sin(p)).to(u.m).value
            self.detector[2] = (r * np.cos(t)).to(u.m).value

        self._log(f"Source: {self.source} m")
        self._log(f"Detector: {self.detector} m")

        # Calculate normal vectors (facing towards the grid origin) for both
        # the source and detector planes
        self.src_n = self.source / np.linalg.norm(
            self.source
        )
        self.det_n = -self.detector / np.linalg.norm(
            self.detector
        )
        # Vector directly from source to detector
        self.src_det_vec = self.detector - self.source

        # Experiment axis is the unit vector from the source to the detector
        self.src_det_n = self.src_det_vec / np.linalg.norm(
            self.src_det_vec
        )

        self.mag = 1 + np.linalg.norm(self.detector)/np.linalg.norm(self.source)

        # Create unit vectors that define the detector plane
        self._create_detector_plane()



    def _log(self, msg):
        if self.verbose:
            print(msg)

    # ************************************************************************
    # Methods used by both integration and particle tracing methods
    # ************************************************************************

    def _create_detector_plane(self):
        r"""
        Defines the horizontal and vertical axes of the detector plane. The
        horizontal axis is defined as being perpendicular to both the
        source-detector axis and the z-axis. In the case where the pos vector
        is aligned with the z-axis, this is automatically chosen to be the
        x-axis. THe vertical axis is then chosen to be orthogonal and
        right-handed with respect to the horizontal axis and the
        source-detector axis.

        """
        # Create 2D grids of detector points
        # Define plane  horizontal axis
        if np.allclose(np.abs(self.det_n), np.array([0, 0, 1])):
            nx = np.array([1, 0, 0])
        else:
            nx = np.cross(np.array([0, 0, 1]), self.det_n)
        nx = nx / np.linalg.norm(nx)
        self.det_hax = nx # Unit vector for hax, detector horizontal axis

        # Define the detector vertical axis as being orthogonal to the
        # detector axis and the horizontal axis
        ny = np.cross(nx, self.det_n)
        ny = -ny / np.linalg.norm(ny)
        self.det_vax = ny # Unit vector for vax, detector vertical axis




    # ************************************************************************
    # Integration Method (linear approximation)
    # ************************************************************************


    def line_integral(self, size=np.array([[-1,1],[-1,1]])*u.cm, bins=[50,50], collimated=True,
                    num=100):
        """
        Analagous to run

        1) Create detector grid from size and bin keywords

        2) For each cell of detector grid, create an array of points of
        separation ds from there to the source point (or, for collimated, the source plane).

        3) Evaluate the integrand at each position in the 3D grid

        4) Integrate along the line-integrated dimension

        Parameters
        ----------
        size : `u.Quantity` array of shape [2,2]
            The bounds of the detector region. The default is [[-1,1],[-1,1]] cm.
        bins : integer ndarray array of shape [2,2]
            Number of bins in each direction of the detector region. The
            default is [50,50].
        collimated : Boolean, optional
            If True, the source will be assumed to be collimated. If False,
            a point source will be used. The default is True (collimated).
        num : int, optional
            Number of integration points along the line (within the grid region).
            The default is 100.

        Returns
        -------
        xax : TYPE
            DESCRIPTION.
        yax : TYPE
            DESCRIPTION.
        integral : TYPE
            DESCRIPTION.

        """

        # Create arrays of detector grid points
        xax = np.linspace(size[0][0].to(u.m).value, size[0][1].to(u.m).value, num=int(bins[0]))
        yax = np.linspace(size[1][0].to(u.m).value, size[1][1].to(u.m).value, num=int(bins[1]))
        x_offset, y_offset = np.meshgrid(xax, yax, indexing='ij')

        # Shift those points in space to be in the detector plane
        det_pts = np.outer(x_offset,self.det_hax) + np.outer(y_offset,self.det_vax) + self.detector
        det_pts = np.reshape(det_pts, [bins[0], bins[1],3])


        # Create 2D grids of source points
        if collimated:
            src_pts = det_pts - self.src_det_vec
        else:
            # If not collimated, assume a point source
            src_pts = np.outer( np.ones([bins[0], bins[1]]), self.source)
            src_pts = np.reshape(src_pts, (bins[0], bins[1], 3))


        # Determine where the grid begins and ends as fractions of the
        # source-to-detector vector
        source_to_det = np.linalg.norm(self.src_det_vec)
        source_to_grid = np.min(np.linalg.norm(self.grid_arr -
                                               self.source, axis=3))
        grid_to_det = np.min(np.linalg.norm(self.grid_arr -
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

        # Evaluate the integrands
        integrands = self.integrand(pts)

        if not isinstance(integrands, tuple):
            integrands = [integrands,]

        # Integrate
        integral = []
        for integrand in integrands:
            integral.append(np.trapz(integrand, axis=2)*(ds*u.m))

        return (xax*u.m, yax*u.m, *integral)


    def integrand(self, pts):
        """
        Returns the integrand at a particular position.

        This method is over-written by subclasses to reflect actual
        diagnostic physics, and can pull plasma parameter information from
        the parameters dict provided.

        Parameters
        ----------

        pts: np.array (bins[0], bins[1], num ,3)
            Positions (in si units) at which the integrand will be
            evaluated. The third axis is the one over which the integrands
            will subsequently be integrated

        """
        raise NotImplementedError("The integrand method must be implemented"
                             " for this diagnostic.")



    




if __name__ == '__main__':


    """
    class TestIntegrator(LineIntegratedDiagnostic):
        def _integrand(self):
            return self.integration_pts[...,0], self.integration_pts[...,1]

    # Make a little grid
    ax = np.linspace(-1,1,3)*u.mm
    xarr, yarr, zarr = np.meshgrid(ax,ax,ax, indexing='ij')
    grid = CartesianGrid(xarr, yarr, zarr)

    source = (0*u.mm,  -3*u.mm, 0*u.mm)
    detector = ( 0*u.mm, 5*u.mm,  0*u.mm)


    obj = TestIntegrator(grid, source, detector)

    hax, vax, integral, _ = obj.evaluate_integral(size=np.array([[-1,1],[-1,1]])*u.mm, bins=[2,2],
                          collimated=False, num=10)

    """






