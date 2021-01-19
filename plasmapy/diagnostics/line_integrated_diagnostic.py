"""
Defines the LineIntegratedDiagnostic base class
"""

__all__ = [
        "LineIntegratedDiagnostic",
        ]

import astropy.units as u
import astropy.constants as const

import numpy as np
import warnings

from plasmapy.plasma.grids import CartesianGrid
from plasmapy.simulation.particle_integrators import boris_push


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

        self.grid = grid
        self.verbose = verbose

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
    def _max_theta_grid(self):
        r"""
        Using the grid and the source position, compute the maximum particle
        theta that will impact the grid. This value can be used to determine
        which particles are worth tracking.
        """
        ind = 0
        theta = np.zeros([8])
        for x in [0, -1]:
            for y in [0, -1]:
                for z in [0, -1]:
                    # Souce to grid corner vector
                    vec = self.grid.grid[x, y, z, :] - self.source

                    # Calculate angle between vec and the source-to-detector
                    # axis, which is the central axis of the proton beam
                    theta[ind] = np.arccos(
                        np.dot(vec.value, self.source_to_detector.value)
                        / np.linalg.norm(vec.value)
                        / np.linalg.norm(self.source_to_detector.value)
                    )
                    ind += 1
        return np.max(theta)


    def _generate_particles(self):
        r"""
        Generates the angular distributions about the Z-axis, then
        rotates those distributions to align with the source-to-detector axis.

        By default, protons are generated over almost the entire pi/2. However,
        if the detector is far from the source, many of these particles will
        never be observed. The max_theta keyword allows these extraneous
        particles to be neglected to focus computational resources on the
        particles who will actually hit the detector.

        """

        self._log("Creating Particles")


        self.charge = const.e.si
        self.mass = const.m_p.si

        # Calculate the velocity corresponding to the proton energy
        #self.v_nonrel = np.sqrt(2 * self.proton_energy / const.m_p.si).to(u.m / u.s)
        ER = (self.proton_energy/(self.mass*const.c.si**2)).to(u.dimensionless_unscaled)
        self.v0 = const.si.c*np.sqrt(1 - 1/(ER+1)**2)


        max_theta = self.max_theta.to(u.rad).value
        # Create a probability vector for the theta distribution
        # Theta must follow a sine distribution in order for the proton
        # flux per solid angle to be uniform.
        arg = np.linspace(0, max_theta, num=int(1e5))
        prob = np.sin(arg)
        prob *= 1 / np.sum(prob)

        # Randomly choose theta's weighted with the sine probabilities
        theta = np.random.choice(arg, size=self.nparticles, replace=True, p=prob)

        # Also generate a uniform phi distribution
        phi = np.random.uniform(size=self.nparticles) * 2 * np.pi

        # Determine the angle aboe which particles will not hit the grid
        # these particles can be ignored until the end of the simulation,
        # then immediately advanced to the detector grid with their original
        # velocities
        max_theta_grid = self._max_theta_grid()

        # This array holds the indices of all particles that WILL hit the grid
        self.grid_ind = np.where(theta < max_theta_grid)[0]
        self.nparticles_grid = len(self.grid_ind)

        # Construct the velocity distribution around the z-axis
        self.v = np.zeros([self.nparticles, 3]) * u.m / u.s
        self.v[:, 0] = self.v0 * np.sin(theta) * np.cos(phi)
        self.v[:, 1] = self.v0 * np.sin(theta) * np.sin(phi)
        self.v[:, 2] = self.v0 * np.cos(theta)

        # Calculate the rotation matrix
        a = np.array([0, 0, 1])
        b = self.detector - self.source
        rot = _rot_a_to_b(a, b)

        # Apply rotation matrix to calculated velocity distribution
        self.v = np.matmul(self.v, rot)

        # Place particles at the source
        self.x = np.outer(np.ones(self.nparticles), self.source)

        # Create flags for tracking when particles during the simulation
        # on_grid -> zero if the particle is off grid, 1
        self.on_grid = np.zeros([self.nparticles_grid])
        # Entered grid -> non-zero if particle EVER entered the grid
        self.entered_grid = np.zeros([self.nparticles_grid])

    def _adaptive_dt(self, Ex, Ey, Ez, Bx, By, Bz):
        r"""
        Calculate the appropraite dt based on a number of considerations
        including the local grid resolution (ds) and the gyroperiod of the
        particles in the current fields.
        """
        # If dt was explicitly set, ignore this fcn
        if self.dt is not None:
            return self.dt

        # Compute the timestep indicated by the grid resolution
        # min is taken for irregular grids, in which different particles
        # have different local grid resolutions
        gridstep = 0.5*(np.min(self.ds) / self.v0).to(u.s)

        #print(f"gridstep = {gridstep.to(u.s):.1e}")

        # If not, compute a number of possible timesteps
        # Compute the cyclotron gyroperiod
        Bmag = np.max( np.sqrt(Bx**2 + By**2 + Bz**2))

        if Bmag == 0:
            gyroperiod = np.inf * u.s
        else:
            gyroperiod = (2 * np.pi * const.m_p.si / (const.e.si * np.max(Bmag))).to(
                u.s
            )

        # TODO: introduce a minimum timestep based on electric fields too!

        # TODO: Clean up the way timestep ranges are specified in the
        # user API.


        # Create an array of all the possible time steps we computed
        candidates = np.array([gyroperiod.value, gridstep.value]) * u.s

        if all(candidates > self.dt_range[1]):
            return self.dt_range[1]

        # Unless it interferes with the range, always choose the smallest
        # time step
        dt = np.min(candidates)

        if dt > self.dt_range[0]:
            return dt
        else:
            return self.dt_range[0]


    def _adaptive_ds(self):
        r"""
        Compute the local grid resolution for each particle (used for determining
        a maximum timestep based on grid spacing).
        """

        # TODO: clean this up? Probably will only accept uniform grids
        # right now because calculating ds for non-uniform grids is a
        # problem in itself.


        if self.grid.regular_grid:
            # If self.ds is a scalar (as setup by the init function)
            # extend it to be the length of npartiles_grid
            if self.ds.size == 1:
                self.ds = self.ds * np.ones(self.nparticles_grid)
            else:
                pass
        else:
            raise NotImplementedError("Adaptive timestep is not yet supportd for "
                                      "non-uniform grids, because the adaptive "
                                      "grid resolution is not supported.")



    def _advance_to_grid(self):
        r"""
        Advances all particles to the timestep when the first particle should
        be entering the grid. Doing in this in one step (rather than pushing
        the particles through zero fields) saves computation time.
        """
        # Distance from the source to the nearest gridpoint
        dist = np.min(np.linalg.norm(self.grid.grid - self.source, axis=3))

        # Time for fastest possible particle to reach the grid.
        t = (dist / self.v0).to(u.s)

        self.x = self.x + self.v * t



    def _advance_to_detector(self):
        r"""
        Advances all particles to the detector plane. This method will be
        called after all particles have cleared the grid.

        This step applies to all particles, including those that never touched
        the grid.
        """
        dist_remaining = np.dot(self.x, self.det_n) + np.linalg.norm(self.detector)

        v_towards_det = np.dot(self.v, -self.det_n)

        # Time remaining for each particle to reach detector plane
        t = dist_remaining / v_towards_det

        # If particles have not yet reached the detector plane and are moving
        # away from it, they will never reach the detector.
        # So, we can remove them from the arrays
        condition = np.logical_and(v_towards_det < 0, dist_remaining > 0)
        ind = np.nonzero(np.where(condition, 0, 1))[0]
        self.x = self.x[ind, :]
        self.v = self.v[ind, :]
        self.nparticles_grid = self.x.shape[0]
        t = t[ind]

        self.x += self.v * np.outer(t, np.ones(3))

        # Check that all points are now in the detector plane
        # (Eq. of a plane is nhat*x + d = 0)
        plane_eq = np.dot(self.x, self.det_n) + np.linalg.norm(self.detector)
        assert np.allclose(plane_eq, np.zeros(self.nparticles_grid), atol=1e-6)


    def _push(self):
        r"""
        Advance particles using an implementation of the time-centered
        Boris algorithm
        """

        # Calculate the local grid resolution for each particle
        self._adaptive_ds()

        pos = self.x[self.grid_ind,:]

        # Update the list of particles on and off the grid
        self.on_grid = grid.on_grid(pos)
        # entered_grid is zero at the end if a particle has never
        # entered the grid
        self.entered_grid += self.on_grid


        # Estimate the E and B fields for each particle
        Ex, Ey, Ez, Bx, By, Bz = self.grid.volume_averaged_interpolator(pos, "E_x", "E_y", "E_z",
                                                                        "B_x", "B_y", "B_y")

        E = np.array([Ex, Ey, Ez])
        E = np.moveaxis(E, 0, -1)
        B = np.array([Bx, By, Bz])
        B = np.moveaxis(B, 0, -1)
        print(E.shape)

        # Calculate the adaptive timestep from the fields currently experienced
        # by the particles
        # If user sets dt explicitly, that's handled in _adpative_dt
        dt = self._adaptive_dt(Ex, Ey, Ez, Bx, By, Bz)


        # TODO: Test v/c and implement relativistic Boris push when required
        #vc = np.max(v)/const.c.si

        boris_push(self.x.si.value,
                   self.v.si.value,
                   B.si.value,
                   E.si.value,
                   self.q.si.value,
                   self.m.si.value,
                   dt.si.value)

    def _stop_condition(self):
        r"""
        The stop condition is that most of the particles have entered the grid
        and almost all have now left it.
        """
        # Count the number of particles who have entered, which is the
        # number of non-zero entries in entered_grid
        n_entered = np.nonzero(self.entered_grid)[0].size

        # How many of the particles have entered the grid
        entered = np.sum(n_entered) / self.nparticles_grid

        # Of the particles that have entered the grid, how many are currently
        # on the grid?
        # if/else avoids dividing by zero
        if np.sum(n_entered) > 0:
            still_on = np.sum(self.on_grid) / np.sum(n_entered)
        else:
            still_on = 0.0

        if entered > 0.1 and still_on < 0.001:
            self._log(
                f"Stop condition reached: {entered*100:.0f}% entered "
                f"({n_entered})"
                f", {still_on*100:.0f}% are still on the grid"
            )

            # Warn user if < 10% of the particles ended up on the grid
            if n_entered < 0.1 * self.nparticles:
                warnings.warn(
                    f"Only {100*n_entered/self.nparticles:.2f}% of "
                    "particles entered the field grid: consider "
                    "decreasing the max_theta to increase this "
                    "number.",
                    RuntimeWarning,
                )

            return True
        else:
            return False



    def run(
        self,
        nparticles,
        max_theta=0.9 * np.pi / 2 * u.rad,
        dt=None,
        dt_range=np.array([0, np.infty]) * u.s,
        field_weighting="nearest neighbor",
    ):
        r"""
        Runs a particle-tracing simulation.
        Timesteps are adaptively calculated based on the
        local grid resolution of the particles and the electric and magnetic
        fields they are experiencing. Both regular (uniform) and irregular
        grids are supported, although the former is faster. After all particles
        have left the simulated field volume, they are advanced to the
        detector plane where they can be used to construct a synthetic
        diagnostic image.

        Parameters
        ----------
        nparticles : integer
            The number of particles to include in the simulation. The default
            is 1e5.

        max_theta : `~astropy.units.Quantity`, optional
            The largest velocity vector angle (measured from the
            source-to-detector axis) for which particles should be generated.
            Decreasing this angle can eliminate particles that would never
            reach the detector region of interest. The default is 0.9*pi/2.
            Units must be convertable to radians.

        dt : `~astropy.units.Quantity`, optional
            An explicitly set timestep in units convertable to seconds.
            Setting this optional keyword overrules the adaptive time step
            capability and forces the use of this timestep throughout.

        dt_range : `~astropy.units.Quantity`, array shape (2,), optional
            A range into which the adaptive dt will be coerced.
            The default is np.array([0, np.infty])*u.s.

        field_weighting : str
            String that selects the field weighting algorithm used to determine
            what fields are felt by the particles. Options are:

            * 'nearest neighbor': Particles are assigned the fields on
                the grid vertex closest to them.

            * 'volume averaged' : The fields experienced by a field are a
                volume-average of the eight grid points surrounding them.

        Returns
        -------
        None.

        """

        # Load inputs
        self.nparticles = int(nparticles)
        self.dt = dt
        self.dt_range = dt_range
        self.field_weighting = field_weighting

        # TODO: Error check that grid contains E and B variables required

        # Advance the particles to the near the start of the simulation
        # volume
        self._advance_to_grid()

        # Push the particles until the stop condition is satisfied
        # (no more particles on the simulation grid)
        while not self._stop_condition():
            if self.verbose:
                fract = 100 * np.sum(self.on_grid) / self.nparticles_grid
                self._log(f"{fract:.1f}% on grid ({np.sum(self.on_grid)})")
            self._push()


        # Advance the particles to the image plane
        # At this stage, remove any particles that are not going to ever
        # hit the grid.
        self._advance_to_detector()

        self._log("Run completed")


    def particle_image(
        self, size=None, bins=None, null=False, optical_density=False
    ):
        r"""
        Calculate a particle count histogram in the
        image plane. The horizontal axis in the detector plane is defined to
        be perpendicular to both the source-to-detector vector and the z-axis
        (unless the source-to-detector axis is parallel to the z axis, in which
        case the horizontal axis is the x-axis). The vertical axis is defined
        to be orthgonal to both the source-to-detector vector and the
        horizontal axis.

        Parameters
        ----------
        size : `~astropy.units.Quantity`, shape (2,2), optional
            The size of the detector array, specified as the minimum
            and maximum values included in both the horizontal and vertical
            directions in the detector plane coordinates. Shape is
            [[hmin,hmax], [vmin, vmax]]. Units must be convertable to meters.

        bins : array of integers, shape (2)
            The number of bins in each direction in the format [hbins, vbins].
            The default is [250,250].

        null: bool
            If True, returns the intensity in the image plane in the absence
            of simulated fields.


        optical_density: bool
            If True, return the optical density rather than the intensity

            .. math::
                OD = -log_{10}(Intensity/I_0)

            where I_O is the intensity on the detector plane in the absence of
            simulated fields. Default is False.

        Returns
        -------
        hax : `~astropy.units.Quantity` array shape (hbins,)
            The horizontal axis of the synthetic radiograph in meters.

        vax : `~astropy.units.Quantity` array shape (vbins, )
            The vertical axis of the synthetic radiograph in meters.

        intensity : ndarray, shape (hbins, vbins)
            The number of protons counted in each bin of the histogram.
        """

        # Note that, at the end of the simulation, all particles were moved
        # into the image plane.

        # Define detector horizontal axis as being perpendicular to both the
        # detector axis and the z-axis. In the case where the detector axis
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
        ny = -ny / np.linalg.norm(ny)

        # If null is True, use the predicted positions in the absence of
        # simulated fields
        if null:
            r = self.r0
        else:
            r = self.r

        # Determine locations of points in the detector plane using unit
        # vectors
        xloc = np.dot(r - self.detector, nx)
        yloc = np.dot(r - self.detector, ny)

        if size is None:
            # If a detector size is not given, choose lengths based on the
            # dimensions of the grid
            w = self.mag * np.max(
                [
                    np.max(np.abs(self.grid.xarr.value)),
                    np.max(np.abs(self.grid.yarr.value)),
                    np.max(np.abs(self.grid.zarr.value)),
                ]
            )

            size = np.array([[-w, w], [-w, w]]) * self.grid.unit

        # If #bins is not set, make a guess
        if bins is None:
            bins = [200, 200]

        # Generate the histogram
        intensity, h, v = np.histogram2d(
            xloc.si.value, yloc.si.value, range=size.si.value, bins=bins
        )

        # Throw a warning if < 50% of the particles are included on the
        # histogram
        percentage = 100 * np.sum(intensity) / self.nparticles
        if percentage < 50:
            warnings.warn(
                f"Only {percentage:.2f}% of the particles are shown "
                " on this synthetic radiograph. Consider increasing "
                " the size to include more.",
                RuntimeWarning,
            )

        if optical_density:
            # Generate the null radiograph
            x, y, I0 = self.particle_image(size=size, bins=bins, null=True)

            # Calcualte I0 as the mean of the non-zero values in the null
            # histogram. Zeros are just outside of the illuminate area.
            I0 = np.mean(I0[I0 != 0])

            # Overwrite any zeros in intensity to avoid log10(0)
            intensity[intensity == 0] = 1

            # Calculate the optical_density
            intensity = -np.log10(intensity / I0)

        return h * u.m, v * u.m, intensity




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




