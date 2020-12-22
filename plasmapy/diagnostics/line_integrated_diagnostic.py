"""
Defines the LineIntegratedDiagnostic base class
"""

__all__ = [
        "LineIntegratedDiagnostic",
        ]

import astropy.units as u
import numpy as np

from plasmapy.plasma.grids import CartesianGrid


class LineIntegratedDiagnostic:


    def __init__(self, grid : u.m,
                 source,
                 detector,
                 geometry = 'cartesian',
                 verbose=True,
                 ):

        self.grid = grid
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


    def evaluate_integral(self, size=np.array([[-1,1],[-1,1]])*u.cm, bins=[50,50], collimated=True,
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
        size : u.quantity array of shape [2,2]
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
        self.hax_n = nx # Unit vector for hax, detector horizontal axis

        # Define the detector vertical axis as being orthogonal to the
        # detector axis and the horizontal axis
        ny = np.cross(nx, self.det_n)
        ny = -ny / np.linalg.norm(ny)
        self.vax_n = ny # Unit vector for vax, detector vertical axis

        xax = np.linspace(size[0][0], size[0][1], num=int(bins[0]))
        yax = np.linspace(size[1][0], size[1][1], num=int(bins[1]))
        x_offset, y_offset = np.meshgrid(xax, yax, indexing='ij')


        det_pts = np.outer(x_offset,nx) + np.outer(y_offset,ny) + self.detector

        det_pts = np.reshape(det_pts, [bins[0], bins[1],3])


        # Create 2D grids of source points
        if collimated:
            src_pts = det_pts - self.src_det_vec
        else:
            # If not collimated, assume a point source
            src_pts = np.outer( np.ones([bins[0], bins[1]]), self.source)
            src_pts = np.reshape(src_pts, (bins[0], bins[1], 3))


        # Calculate the unit vector of integration at each point
        # for use in integrands later
        lines = det_pts - src_pts
        lines = np.moveaxis(lines, -1,0)
        unit_v = lines/np.linalg.norm(lines, axis=0)
        unit_v = np.outer(unit_v, np.ones(num))
        unit_v = np.reshape(unit_v, (3, bins[0], bins[1], num))
        unit_v = np.moveaxis(unit_v, 0, 3)
        self.unit_v = unit_v

        # Determine where the grid begins and ends as fractions of the
        # source-to-detector vector
        source_to_det = np.linalg.norm(self.src_det_vec).to(u.mm)
        source_to_grid = np.min(np.linalg.norm(self.grid.grid*self.grid.unit -
                                               self.source, axis=3))
        grid_to_det = np.min(np.linalg.norm(self.grid.grid*self.grid.unit -
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
        self.integration_pts = pts

        integrands = self._integrand()

        if not isinstance(integrands, tuple):
            integrands = [integrands,]

        # Integrate
        integral = []
        for integrand in integrands:
            integral.append(np.trapz(integrand, axis=2)*ds)

        return (xax, yax, *integral)


    def _integrand(self):
        """
        Returns the integrand at a particular position.

        This method is over-written by subclasses to reflect actual
        diagnostic physics, and can pull plasma parameter information from
        the parameters dict provided.
        """
        raise NotImplementedError("The integrand method must be implemented"
                             " for this diagnostic.")



    # ************************************************************************
    # Particle Tracing Method (for non-linear problems)
    # ************************************************************************


    def _generate_particles(self):
        """
        Generate particles

        """


    def _advance_to_grid(self):
        """
        Advance particles to close to the start of the grid
        """

    def run(self):
        """
        Run the particle tracer
        """


        # Create an instance of plasmapy's ParticleTracker using the
        # particles previously generated


        # Run the particle tracker


    def _advance_to_detector(self):
        """
        Advance particles to the detector plane and discard any particles that
        will never reach the detector plane (eg. moving away)
        """


    def create_image(self):
        """
        Create a histogram image of particle locations in the image plane
        """




if __name__ == '__main__':


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




