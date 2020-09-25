"""
Routines for the analysis of proton radiographs. These routines can be broadly
classified as either creating synthetic radiographs from prescribed fields or
methods of 'inverting' experimentally created radiographs to reconstruct the
original fields (under some set of assumptions).
"""

__all__ = [
    "SyntheticProtonRadiograph",
]

import astropy.constants as const
import astropy.units as u
import numpy as np
import scipy.interpolate as interp
import warnings

from abc import ABC, abstractmethod
from scipy.special import erf as erf

from plasmapy.plasma import fields as fields


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


class SyntheticProtonRadiograph:
    r"""
    Represents a proton radiography experiment with simulated or
    calculated E and B fields given at positions defined by a grid of spatial
    coordinates. The proton source and detector plane are defined by vectors
    from the origin of the field grid, and the energy of the protons
    can be set.

    Parameters
    ----------
    grid : `~astropy.units.Quantity`, shape (nx,ny,nz,3)
        An array giving the positions of each grid point. Units must be
        convertable to meters.

    E : `~astropy.units.Quantity`, shape (nx,ny,nz,3)
        The vector electric field at each gridpoint. Units must be
        convertable to V/m.

    B : `~astropy.units.Quantity`, shape (nx,ny,nz,3)
        The vector magnetic field at each grid point. Units must be
        convertable to Tesla.

    source : `~astropy.units.Quantity`, shape (3)
        A vector pointing from the origin of the field grid to the location
        of the proton point source. This vector will be interpreted as
        being in either cartesian, cylindrical, or spherical coordinates
        based on the geometry keyword. The units of the vector must be
        compatible with the geometry chosen:
        * Cartesian (x,y,z) : (meters, meters, meters)
        * cylindrical (r, theta, z) : (meters, radians, meters)
        * spherical (r, theta, phi) : (meters, radians, radians)
        In spherical coordinates theta is the polar angle.

    detector : `~astropy.units.Quantity`, shape (3)
        A vector pointing from the origin of the field grid to the center
        of the detector plane. The vector from the source point to this
        point defines the normal vector of the detector plane. This vector
        can also be specified in cartesian, cylindrical, or spherical
        coordinates by setting the geometry keyword and using the
        appropriate units.

    proton_energy : `~astropy.units.Quantity`, optional
        The energy of the protons, convertable to eV. The default is
        14 MeV.

    geometry : string, optional
        A keyword that allows the source and detector vectors to be
        specified in different coordinate systems. Valid values are
        'cartesian', 'cylindrical', and 'spherical'. The default
        value is 'cartesian'.

    verbose : bool, optional
        If true, updates on the status of the program will be printed
        into the command line while running.
    """

    def __init__(
        self,
        grid: u.m,
        E: u.V / u.m,
        B: u.T,
        source: u.m,
        detector: u.m,
        proton_energy=14 * u.MeV,
        geometry="cartesian",
        verbose=True,
    ):
        r"""
        Initalize the simPrad object, carry out coordinate transformations,
        and compute several quantities that will be used elsewhere.
        """

        # Convert grid to PosGrid object
        if isinstance(grid, fields.PosGrid):
            self.grid = grid
        else:
            self.grid = fields.PosGrid(grid=grid)

        self.E = E
        self.B = B
        self.proton_energy = proton_energy
        self.verbose = verbose

        # Validate input arrays
        arr = {"grid": self.grid.grid, "E": E, "B": B}
        for x in arr.keys():
            if not np.isfinite(arr[x]).all():
                raise ValueError(
                    f"Input arrays must be finite: {x} contains "
                    "either NaN or infinite values."
                )

        self.charge = const.e.si
        self.mass = const.m_p.si
        # Calculate the velocity corresponding to the proton energy
        self.v0 = np.sqrt(2 * self.proton_energy / const.m_p.si).to(u.m / u.s)

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

        # Check that source-detector vector actually passes through the grid
        if not self.grid.vector_intersects(self.source, self.detector):
            raise ValueError(
                "The vector from the source to the detector "
                "does not pass through the grid! "
                f"Source: {self.source}, "
                f"Detector: {self.detector}"
            )

        # Calculate some parameters involving the source and detector locations
        self.det_n = -self.detector.si.value / np.linalg.norm(
            self.detector.si.value
        )  # Plane normal vec

        # Vector directly from source to detector
        self.source_to_detector = self.detector - self.source
        # Experiment axis is the unit vector from the source to the detector
        self.exp_ax = self.source_to_detector.si.value / np.linalg.norm(
            self.source_to_detector.si.value
        )

        # Compute the magnification
        self.mag = 1 + (
            np.linalg.norm(self.detector.si.value)
            / np.linalg.norm(self.source.si.value)
        )

    def _log(self, msg):
        if self.verbose:
            print(msg)

    def _init_interpolator(self):
        r""""
        Auto-detects whether the given grid is uniform or irregular and
        creates an appropriate interpolator. Also calculates the grid
        resolution for use in choosing an appropriate timestep.
        """
        # Create a grid of indices for use in interpolation
        nx, ny, nz, x = self.grid.shape
        indgrid = np.indices([nx, ny, nz])
        indgrid = np.moveaxis(indgrid, 0, -1)

        if self.grid.regular_grid:
            self._log("Auto-detected a regularly spaced grid.")
        else:
            self._log("Auto-detected an irregularly spaced grid.")

        if self.grid.regular_grid:
            # Create axes under the regular grid assumption
            dvec = (
                np.array([self.grid.dx.value, self.grid.dy.value, self.grid.dz.value])
                * self.grid.dx.unit
            )

            # Estimate the grid-point spacing along the source_to_detector vector
            # Will be expanded to a constant array of length nparticles_grid
            # when _adaptive_ds() is called.
            self.ds = np.linalg.norm(np.dot(dvec, self.exp_ax))

            # Initialize the interpolator
            pts = (
                self.grid.xaxis.si.value,
                self.grid.yaxis.si.value,
                self.grid.zaxis.si.value,
            )
            self._log("Creating regular grid interpolator")
            self.interpolator = interp.RegularGridInterpolator(
                pts, indgrid, method="nearest", bounds_error=False, fill_value=-1
            )

        else:
            # Flat arrays of points for irregular grid interpolation fcn
            pts = np.zeros([nx * ny * nz, 3])
            pts[:, 0] = self.grid.xarr.flatten().si.value
            pts[:, 1] = self.grid.yarr.flatten().si.value
            pts[:, 2] = self.grid.xarr.flatten().si.value

            # Flatten the index grid for irregular grid interpolation fcn
            indgrid2 = np.zeros([nx * ny * nz, 3])
            indgrid2[:, 0] = indgrid[:, :, :, 0].flatten()
            indgrid2[:, 1] = indgrid[:, :, :, 1].flatten()
            indgrid2[:, 2] = indgrid[:, :, :, 2].flatten()

            # Initialize the interpolator
            self._log("Creating irregular grid interpolator")
            self.interpolator = interp.NearestNDInterpolator(pts, indgrid2)

            # If dt is not explicitly set, create an array of the
            # distance to the nearest neighbor of each grid
            if self.dt is None:
                # TODO
                # self.ds is used for determining when particles are on-grid
                # Somewhat ambiguous how to chose a single value for this for an
                # irregular grid: this may not be the best solution.
                self.ds = np.median(self.grid.nearest_neighbor)

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

    def _generate_particles(self, max_theta=0.9 * np.pi / 2 * u.rad):
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
        max_theta = max_theta.to(u.rad).value

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
        self.gi = np.where(theta < max_theta_grid)[0]
        self.nparticles_grid = len(self.gi)

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
        self.r = np.outer(np.ones(self.nparticles), self.source)

        # Create flags for tracking when particles during the simulation
        # on_grid -> zero if the particle is off grid, 1
        self.on_grid = np.zeros([self.nparticles_grid])
        # Entered grid -> non-zero if particle EVER entered the grid
        self.entered_grid = np.zeros([self.nparticles_grid])

    def _generate_null(self):
        r"""
        Calculate the distribution of particles on the detector in the absence
        of any simulated fields.
        """
        # Calculate the unit vector from the source to the detector
        dist = np.linalg.norm(self.source_to_detector)
        uvec = self.source_to_detector.to(u.m).value / dist.to(u.m).value

        # Calculate the remaining distance each particle needs to travel
        # along that unit vector
        remaining = np.dot(self.source, uvec)

        # Calculate the time remaining to reach that plane and push
        t = ((dist - remaining) / np.dot(self.v, uvec)).to(u.s)

        # Calculate the particle positions for that case
        self.r0 = self.source + self.v * np.outer(t, np.ones(3))

    def _advance_to_grid(self):
        r"""
        Advances all particles to the timestep when the first particle should
        be entering the grid (to save time)
        """
        # Distance from the source to the nearest gridpoint
        dist = np.min(np.linalg.norm(self.grid.grid - self.source, axis=3))

        # Time for fastest possible particle to reach the grid.
        t = (dist / self.v0).to(u.s)

        self.r = self.r + self.v * t

    def _advance_to_detector(self):
        r"""
        Once all particles have cleared the grid, advance them to the detector
        plane.

        This step applies to all particles, not just those that will hit the
        grid.
        """
        dist_remaining = np.dot(self.r, self.det_n) + np.linalg.norm(self.detector)

        v_towards_det = np.dot(self.v, -self.det_n)

        # Time remaining for each particle to reach detector plane
        t = dist_remaining / v_towards_det

        # If particles have not yet reached the detector plane and are moving
        # away from it, they will never reach the detector.
        condition = np.logical_and(v_towards_det < 0, dist_remaining > 0)
        ind = np.nonzero(np.where(condition, 0, 1))[0]
        self.r = self.r[ind, :]
        self.v = self.v[ind, :]
        self.nparticles_grid = self.r.shape[0]
        t = t[ind]

        self.r += self.v * np.outer(t, np.ones(3))

        # Check that all points are now in the detector plane
        # (Eq. of a plane is nhat*x + d = 0)
        plane_eq = np.dot(self.r, self.det_n) + np.linalg.norm(self.detector)
        assert np.allclose(plane_eq, np.zeros(self.nparticles_grid), atol=1e-6)

    def _place_particles(self):
        r"""
        For each particle, find the indicies of the nearest field assuming that
        the fields are placed on a regular grid
        """
        # Interpolate the grid indices that best match each particle position
        i = self.interpolator(self.r[self.gi, :].si.value)

        # Store the grid positions
        self.xi = i[:, 0].astype(np.int32)
        self.yi = i[:, 1].astype(np.int32)
        self.zi = i[:, 2].astype(np.int32)

    def _estimate_fields(self):
        """
        Return the field experienced by each particle.

        If field_weighting = 'volume averaged', use a volume-averaged field
        weighting scheme. If not, use a neareset-neighbor scheme.
        """

        if self.field_weighting == "nearest neighbor":
            E = self.E[self.xi, self.yi, self.zi, :]
            B = self.B[self.xi, self.yi, self.zi, :]

            # Set fields for off-grid particles
            E = E * np.outer(self.on_grid, np.ones(3))
            B = B * np.outer(self.on_grid, np.ones(3))

        elif self.field_weighting == "volume averaged":
            if not self.grid.regular_grid:
                raise ValueError(
                    "Volume weighting fields is currently "
                    " only supported on regular grids."
                )

            # Load the positions of each particle, and of the grid point
            # interpolated for it
            rpos = self.r[self.gi, :]
            gpos = self.grid.grid[self.xi, self.yi, self.zi, :]

            nx, ny, nz, nax = self.grid.shape

            # Determine the points bounding the grid cell containing the
            # particle
            x0 = np.where(rpos[:, 0] > gpos[:, 0], self.xi, self.xi - 1)
            x1 = x0 + 1
            y0 = np.where(rpos[:, 1] > gpos[:, 1], self.yi, self.yi - 1)
            y1 = y0 + 1
            z0 = np.where(rpos[:, 2] > gpos[:, 2], self.zi, self.zi - 1)
            z1 = z0 + 1

            # Calculate the cell volume
            # TODO: if this could be done for an arbitrary mesh of 8 points,
            # this algorithm could be used on irregualr grids. Make it so?
            cell_vol = self.grid.dx * self.grid.dy * self.grid.dz

            E = np.zeros([len(self.xi), 3]) * self.E.unit
            B = np.zeros([len(self.xi), 3]) * self.B.unit

            for x in [x0, x1]:
                for y in [y0, y1]:
                    for z in [z0, z1]:

                        # Determine if gridpoint is within bounds
                        valid = (
                            (x >= 0)
                            & (x < nx)
                            & (y >= 0)
                            & (y < ny)
                            & (z >= 0)
                            & (z < nz)
                        )
                        out = np.where(valid == False)

                        # Distance from grid vertex to particle position
                        d = np.abs(self.grid.grid[x, y, z, :] - rpos)

                        # Fraction of cell volume that is closest to the
                        # current point
                        weight = d[:, 0] * d[:, 1] * d[:, 2] / cell_vol
                        weight[out] = 0
                        weight = np.outer(weight, np.ones(3))

                        E += weight * self.E[x, y, z, :]
                        B += weight * self.B[x, y, z, :]

        return E, B

    def _adaptive_dt(self, B):
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
        gridstep = (np.min(self.ds) / self.v0).to(u.s)

        # If not, compute a number of possible timesteps
        # Compute the cyclotron gyroperiod
        Bmag = np.max(np.linalg.norm(B, axis=1))  # B is [nparticles,3] here
        if Bmag == 0:
            gyroperiod = np.inf * u.s
        else:
            gyroperiod = (2 * np.pi * const.m_p.si / (const.e.si * np.max(Bmag))).to(
                u.s
            )

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
        Compute the local grid resolution for each particle (for determining
        if the particle should be influenced by the grid or not). For regular
        grids this is a constant array, but for irregular grids it changes
        as the particles move
        """

        if self.grid.regular_grid:
            # If self.ds is a scalar (as setup by the init function)
            # extend it to be the length of npartiles_grid
            if self.ds.size == 1:
                self.ds = self.ds * np.ones(self.nparticles_grid)
            else:
                pass
        else:
            # Calculate ds for each particle as the nearest-neighbor
            # distance of its assigned gridpoint
            self.ds = self.grid.nearest_neighbor[self.xi, self.yi, self.zi]

    def _push(self):
        r"""
        Advance particles using an implementation of the time-centered
        Boris algorithm
        """
        # Calculate the indices of the field grid points nearest to each particle
        # Note that this is the most time-intensive part of each push
        self._place_particles()

        # Calculate the local grid resolution for each particle
        self._adaptive_ds()

        # Update the list of particles on and off the grid
        dist = np.linalg.norm(
            self.r[self.gi, :] - self.grid.grid[self.xi, self.yi, self.zi, :], axis=1
        )
        self.on_grid = np.where(dist < self.ds, 1, 0)
        self.entered_grid += self.on_grid

        # Estimate the E and B fields for each particle
        E, B = self._estimate_fields()

        # Calculate the adaptive timestep from the fields currently experienced
        # by the particles
        dt = self._adaptive_dt(B)
        dt2 = dt * self.charge / self.mass / 2

        # Push only particles on a grid trajectory
        v = self.v[self.gi, :]

        # Execute the Boris push algorithm
        vminus = v + E * dt2
        t = -B * dt2
        s = 2 * t / (1 + (t * t).sum(axis=1, keepdims=True))
        vprime = vminus + np.cross(vminus.si.value, t) * u.m / u.s
        vplus = vminus + np.cross(vprime.si.value, s) * u.m / u.s
        vnew = vplus + E * dt2

        # Update the velocities of the particles that are being pushed
        self.v[self.gi, :] = vnew
        # Update the positions
        self.r[self.gi, :] += self.v[self.gi, :] * dt

    def run(
        self,
        nparticles,
        max_theta=0.9 * np.pi / 2 * u.rad,
        dt=None,
        dt_range=np.array([0, np.infty]) * u.s,
        field_weighting="nearest neighbor",
    ):
        r"""
        Runs a particle-tracing simulation using the geometry defined in the
        SimPrad object. Timesteps are adaptively calculated based on the
        local grid resolution of the particles and the electric and magnetic
        fields they are experiencing. Both regular (uniform) and irregular
        grids are supported, although the former is faster. After all particles
        have left the simulated field volume, they are advanced to the
        detector plane where they can be used to construct a synthetic
        proton radiograph.

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

        # Initialize the interpolators and grid resolution matrices
        self._init_interpolator()

        # Initialize variables and create the particle distribution
        self._generate_particles(max_theta=max_theta)

        # Generate a null distribution (where the particles would go without
        # simulated fields)
        self._generate_null()

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

    def calc_ke(self, total=True):
        r"""
        Calculate the total kinetic energy of some or all particles. This calculation
        is currently done on velocity time steps (half-integer time steps)
        but it's good enough for ensuring energy is conserved.'
        """
        ke = 0.5 * self.mass * np.sum(self.v ** 2, axis=1)

        if total:
            return np.sum(ke).to(u.J)
        else:
            return ke.to(u.J)

    def synthetic_radiograph(
        self, size=None, bins=None, null=False, optical_density=False
    ):
        r"""
        Calculate a "synthetic radiograph" (particle count histogram in the
        image plane). The horizontal axis in the detector plane is defined to
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
        ny = ny / np.linalg.norm(ny)

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
            x, y, I0 = self.synthetic_radiograph(size=size, bins=bins, null=True)

            # Calcualte I0 as the mean of the non-zero values in the null
            # histogram. Zeros are just outside of the illuminate area.
            I0 = np.mean(I0[I0 != 0])

            # Overwrite any zeros in intensity to avoid log10(0)
            intensity[intensity == 0] = 1

            # Calculate the optical_density
            intensity = -np.log10(intensity / I0)

        return h * u.m, v * u.m, intensity
