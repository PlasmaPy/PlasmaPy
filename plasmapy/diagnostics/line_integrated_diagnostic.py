"""
Defines the LineIntegratedDiagnostic base class and subclasses
"""

__all__ = [
    "LineIntegratedDiagnostic",
    "LineIntegrateScalarQuantities",
    "Interferometer",
]

import astropy.constants as const
import astropy.units as u
import numpy as np

from typing import Union

from plasmapy.diagnostics.proton_radiography import _coerce_to_cartesian_si


class LineIntegratedDiagnostic:
    """
    An abstract line integrated diagnostic
    """

    def __init__(
        self,
        grid: u.m,
        source,
        detector,
        verbose=True,
    ):

        # self.grid is the grid object
        self.grid = grid
        # self.grid_arr is the grid positions in si. This is created here
        # so that it isn't continously called later
        self.grid_arr = grid.grid.to(u.m).value

        self.verbose = verbose

        self.source = _coerce_to_cartesian_si(source)
        self.detector = _coerce_to_cartesian_si(detector)

        self._log(f"Source: {self.source} m")
        self._log(f"Detector: {self.detector} m")

        # Calculate normal vectors (facing towards the grid origin) for both
        # the source and detector planes
        self.src_n = self.source / np.linalg.norm(self.source)
        self.det_n = -self.detector / np.linalg.norm(self.detector)
        # Vector directly from source to detector
        self.src_det_vec = self.detector - self.source

        # Experiment axis is the unit vector from the source to the detector
        self.src_det_n = self.src_det_vec / np.linalg.norm(self.src_det_vec)

        self.mag = 1 + np.linalg.norm(self.detector) / np.linalg.norm(self.source)

        # Create unit vectors that define the detector plane
        self._create_detector_plane()

    def _log(self, msg):
        if self.verbose:
            print(msg)

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
        self.det_hax = nx  # Unit vector for hax, detector horizontal axis

        # Define the detector vertical axis as being orthogonal to the
        # detector axis and the horizontal axis
        ny = np.cross(nx, self.det_n)
        ny = -ny / np.linalg.norm(ny)
        self.det_vax = ny  # Unit vector for vax, detector vertical axis

    def line_integral(
        self,
        size=np.array([[-1, 1], [-1, 1]]) * u.cm,
        bins=[50, 50],
        collimated=True,
        num=100,
    ):
        """
        Calculates the line-integral tof the integrand function through the
        provided grid. This is accomplished through the following steps:

        1) Create an array of points in the detector plane using the detector
        location and the size and bins keywords.

        2) For each cell of detector grid, create an array of points of
        separation ds from there to the source point (or, when collimated=True,
                                                      the source plane).

        3) Evaluate the integrand function at each point.

        4) Integrate along the line-integrated dimension to obtain the
        line-integrated quantity in the detector plane.

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
        xax : `~u.Quantity` array (Nh,)
            The horizontal axis of the detector plane

        yax : `~u.Quantity` array of shape (Nv,)
            The vertical axis of the detector plane

        integral : `~u.Quantity` array or list of arrays of shape (Nh, Ny)
           The line-integrated values in the detector plane.

        """

        # Create arrays of detector grid points
        xax = np.linspace(
            size[0][0].to(u.m).value, size[0][1].to(u.m).value, num=int(bins[0])
        )
        yax = np.linspace(
            size[1][0].to(u.m).value, size[1][1].to(u.m).value, num=int(bins[1])
        )
        x_offset, y_offset = np.meshgrid(xax, yax, indexing="ij")

        # Shift those points in space to be in the detector plane
        det_pts = (
            np.outer(x_offset, self.det_hax)
            + np.outer(y_offset, self.det_vax)
            + self.detector
        )
        det_pts = np.reshape(det_pts, [bins[0], bins[1], 3])

        # Create 2D grids of source points
        if collimated:
            src_pts = det_pts - self.src_det_vec
        else:
            # If not collimated, assume a point source
            src_pts = np.outer(np.ones([bins[0], bins[1]]), self.source)
            src_pts = np.reshape(src_pts, (bins[0], bins[1], 3))

        # Determine where the grid begins and ends as fractions of the
        # source-to-detector vector
        source_to_det = np.linalg.norm(self.src_det_vec)
        source_to_grid = np.min(np.linalg.norm(self.grid_arr - self.source, axis=3))
        grid_to_det = np.min(np.linalg.norm(self.grid_arr - self.detector, axis=3))

        # Paramter for parameteric equation
        start = source_to_grid / source_to_det
        stop = 1 - grid_to_det / source_to_det
        i = np.linspace(start, stop, num=num)
        ds = (stop - start) * source_to_det / num

        # Create an array of points along lines between points in src_pts
        # and det_pts of the equation m*i + b where
        # m = det_pts - src_pts
        # Integrating along the lines corresponds to summation along ax=2
        mi = np.outer(det_pts - src_pts, i)
        mi = np.reshape(mi, [bins[0], bins[1], 3, num])

        b = np.outer(src_pts, np.ones(num))
        b = np.reshape(b, [bins[0], bins[1], 3, num])

        pts = (mi + b) * u.m
        pts = np.moveaxis(pts, 2, 3)

        # Reshape the pts array from grid shape (nx, ny, nz, 3) to a list
        # of points (nx*ny*nz, 3) as required by the grids interpolators
        nx, ny, nz, ndim = pts.shape
        pts = np.reshape(pts, (nx * ny * nz, ndim))

        # Evaluate the integrands
        integrands = self.integrand(pts)

        # If a single integrand is returned, put it in a list
        if not isinstance(integrands, tuple):
            integrands = (integrands,)

        # Integrate
        integral = []
        for integrand in integrands:
            # Reshape the integrands from (nx*ny*nz) to (nx, ny, nz)
            integrand = np.reshape(integrand, (nx, ny, nz))
            # Integrate
            integral.append(np.trapz(integrand, axis=2) * (ds * u.m))

        return (xax * u.m, yax * u.m, *integral)

    def integrand(self, pts):
        """
        Returns the integrand at a particular position.

        This method is over-written by subclasses to reflect actual
        diagnostic physics, and can pull plasma parameter information from
        the parameters dict provided.

        Parameters
        ----------

        pts: `~u.Quantity` (nx*ny*nz, 3)
            Positions at which the integrand will be
            evaluated.

        Returns
        -------

        arr: `~u.Quantity`(nx*ny*nz) or list of same
            Integrand value at each of the points provided. Some integrand
            functions may return multiple arrays of interpolated values as
            a list, each of which will then be integrated separately.
        """
        raise NotImplementedError(
            "The integrand method must be implemented" " for this diagnostic."
        )


class LineIntegrateScalarQuantities(LineIntegratedDiagnostic):
    """
    Line-integrates a scalar quantity
    """

    def __init__(
        self,
        grid: u.m,
        source,
        detector,
        quantities: Union[str, list, tuple],
        verbose=True,
    ):

        # Validate the quantities input
        if isinstance(quantities, str):
            quantities = [
                quantities,
            ]

        for quantity in quantities:
            if quantity not in grid.quantities:
                raise ValueError(
                    f"quantity {quantity} is not defined on the " "provided grid."
                )

        self.quantities = quantities

        # Continue with the rest of the parent class init
        super().__init__(grid, source, detector, verbose=verbose)

    def integrand(self, pts):
        """
        Returns the scalar value of the quantity or quantities specified
        at each point.

        Parameters
        ----------

        pts: `~u.Quantity` (nx*ny*nz, 3)
            Positions at which the integrand will be
            evaluated.

        Returns
        -------

        arr: `~u.Quantity`(nx*ny*nz) or list of same
            Integrand value at each of the points provided. Some integrand
            functions may return multiple arrays of interpolated values as
            a list, each of which will then be integrated separately.

        """
        arr = self.grid.volume_averaged_interpolator(pts, *self.quantities)
        return arr


class Interferometer(LineIntegrateScalarQuantities):
    def __init__(self, grid, source, detector, verbose=False):
        super().__init__(grid, source, detector, quantities="n_e", verbose=verbose)

    def interferogram(
        self,
        probe_freq: u.Hz,
        size=np.array([[-1, 1], [-1, 1]]) * u.cm,
        bins=[50, 50],
        collimated=True,
        num=100,
        unwrapped=True,
    ):

        # TODO: implement an actual critical density function for PlasmaPy
        # Critical density in cm^-3
        n_c = (
            (const.eps0.si * const.m_e / const.e.si ** 2)
            * (2 * np.pi) ** 2
            * probe_freq ** 2
        )
        n_c = n_c.to(u.cm ** -3)

        hax, vax, int_ne = self.line_integral(
            size=size, bins=bins, collimated=collimated, num=num
        )

        phase_shift = (-np.pi * probe_freq / (const.c.si * n_c)) * int_ne

        phase_shift = phase_shift.to(u.dimensionless_unscaled).value

        if not unwrapped:
            # Thanks to helpful stack exchange answer for this compact expression
            # https://stackoverflow.com/a/15927914
            phase_shift = (phase_shift + np.pi) % (2 * np.pi) - np.pi

        return hax, vax, phase_shift
