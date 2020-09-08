"""
Routines for the analysis of proton radiographs. These routines can be broadly
classified as either creating synthetic radiographs from prescribed fields or
methods of 'inverting' experimentally created radiographs to reconstruct the
original fields (under some set of assumptions).
"""

__all__ = [
    "test_fields",
    "SimPrad",
]

import astropy.constants as const
import astropy.units as u
import numpy as np
import scipy.interpolate as interp
import time


def test_fields(
    grid=None,
    mode="electrostatic gaussian sphere",
    regular_grid=True,
    num=(100, 100, 100),
    length=1 * u.mm,
):
    """
    This function generates test fields based on analytical models for
    testing or demonstrating the proton radiography module.

    Parameters
    ----------
    grid : `~astropy.units.Quantity` array shape [nx,ny,nz,3], optional
        A grid of positions on which to calculate fields. If no grid is
        specified, one will be created.

    mode : str, optional
        The field model to load. The default represents the electrostatic
        E-field created by a spherical Gaussian potential.

    regular_grid : bool, optional
        If a grid is being generated, setting this keyword to False will
        generate a grid with random spacing between points. If True
        (the default) the grid spacing will be uniform.

    num : int or tuple of 3 ints, optional
        If a grid is being constructed, this variable sets the number of
        elements in each dimension. If set to a single integer, all three
        dimensions are taken to have the same length. The default is 100.

    length : `~astropy.units.Quantity`, optional
        The length of each dimension, which can be given in several ways:
            - If a single value is given, then the length of each dimension
            will be set to (-length, length)
            - If an array of shape (3) is given, each dimension is set to be
            symmetric using each value, eg. xdim = [-length[0], length[0]]...
            - If an array of shape (2,3) is given, then L[:,i] is the min
            and max of the ith dimension.

    Returns
    -------
    grid : `~astropy.units.Quantity` array, shape [nx,ny,nz,3]
        The grid of positions corresponding to the field grids.

    E : `~astropy.units.Quantity` array, shape [nx,ny,nz,3]
        Electric field array in units of V/m

    B : `~astropy.units.Quantity` array, shape [nx,ny,nz,3]
        Magnetic field array in units of Tesla

    """

    # TODO: Not all fields are compatible with non-uniform grids: notably the
    # electric gaussian sphere relies on a gradient...

    # Process different options for inputs
    if isinstance(num, int):
        num = (num, num, num)

    L = np.zeros([2, 3]) * length.unit
    if isinstance(L, u.Quantity):
        L[:, 0] = np.array([-length.value, length.value]) * length.unit
        L[:, 1] = np.array([-length.value, length.value]) * length.unit
        L[:, 2] = np.array([-length.value, length.value]) * length.unit
    elif L.shape == (3):
        L[:, 0] = np.array([-length[0].value, length[0].value]) * length.unit
        L[:, 1] = np.array([-length[1].value, length[1].value]) * length.unit
        L[:, 2] = np.array([-length[2].value, length[2].value]) * length.unit
    else:
        L = length

    # If no grid is specified, create a grid
    if grid is None:
        grid = np.zeros([num[0], num[1], num[2], 3]) * u.cm

        # If regulr_grid keyword is False, create a non-uniform grid with
        # random point spacing in all directions
        if not regular_grid:
            for d in [0, 1, 2]:
                ax = np.random.uniform(low=L[0, d].value, high=L[1, d].value, size=num)
                ax = np.sort(ax, axis=d)
                grid[..., d] = ax * L.unit

        # If regular_grid keyword is set, create a uniformly spaced grid
        else:
            xaxis = np.linspace(L[0, 0], L[1, 0], num=num[0])
            yaxis = np.linspace(L[0, 1], L[1, 1], num=num[1])
            zaxis = np.linspace(L[0, 2], L[1, 2], num=num[2])
            xarr, yarr, zarr = np.meshgrid(xaxis, yaxis, zaxis, indexing="ij")
            grid[..., 0] = xarr
            grid[..., 1] = yarr
            grid[..., 2] = zarr

    # Create axes and arrays of each coordinate from the grid
    xaxis = grid[:, 0, 0, 0]
    yaxis = grid[0, :, 0, 1]
    zaxis = grid[0, 0, :, 2]
    xarr = grid[..., 0]
    yarr = grid[..., 1]
    zarr = grid[..., 2]

    # Calculate length of each axis, for use in scaling field features
    length = np.abs(L[0, :] - L[1, :])

    # Calculate radius arrays in spherical and cylindrical coordinates
    # for use in generating different field structures
    radius = np.sqrt(xarr ** 2 + yarr ** 2 + zarr ** 2)
    pradius = np.sqrt(xarr ** 2 + yarr ** 2)

    # Create blank arrays for fields
    E = np.zeros(grid.shape) * u.V / u.m
    B = np.zeros(grid.shape) * u.T

    if mode == "no fields":
        pass

    elif mode == "electrostatic gaussian sphere":
        print("Generating Electrostatic Gaussian Sphere")

        a = length[0] / 3
        potential = -np.exp(-(xarr ** 2 + yarr ** 2 + zarr ** 2) / a ** 2) * u.V

        Ex, Ey, Ez = np.gradient(potential, xaxis, yaxis, zaxis)

        Ex = np.where(radius < 0.5 * length[0], Ex, 0)
        Ey = np.where(radius < 0.5 * length[0], Ey, 0)
        Ez = np.where(radius < 0.5 * length[0], Ez, 0)

        E[:, :, :, 0], E[:, :, :, 1], E[:, :, :, 2] = Ex, Ey, Ez

        # Normalize to a desired maximum field
        E = (E / np.max(E)).to(u.dimensionless_unscaled)
        E = E * (5e9 * u.V / u.m)

    elif mode == "axial magnetic field":
        print("Generating Axial Magnetic Field")

        a = length[0] / 4
        B[:, :, :, 2] = np.where(pradius < a, 400 * u.T, 0 * u.T)

    return grid, E, B


def _rot_a_to_b(a, b):
    """
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


def _nearest_neighbor(array):
    """
    Given a 3D array of cartesian postions of shape [nx,ny,nz,3],
    return a Quantity array of the distance to the nearest neighbor from each
    point of shape [nx,ny,nz]. This output is used to define the local
    scalar grid resolution on an irregular grid.
    """
    nx, ny, nz, x = array.shape

    dist = np.zeros([nx, ny, nz, 26]) * array.unit

    ind = 0
    for x in [-1, 0, 1]:
        for y in [-1, 0, 1]:
            for z in [-1, 0, 1]:

                # Skip the zero shift case.
                if x == 0 and y == 0 and z == 0:
                    continue

                # Shift the array
                shifted_arr = np.roll(array, shift=(x, y, z), axis=(0, 1, 2))
                # Compute the distance between the shifted array points and the
                # previous array points
                dist[..., ind] = np.linalg.norm(array - shifted_arr, axis=3)
                # Increment counter
                ind += 1

    # Return the minimum distance between all 27 points
    return np.min(dist, axis=3)


# TODO: error checking
# ERRORS
# - Validate that no input vaules contain nans or anything like that

# WARNINGS
# - Throw a warning if source-detector vector does not pass through grid volume
# - Throw a warning if highly non-linear deflections (path crossings) are expected
# - Throw a warning if <10% of the particles ended up on the grid by the end of the run
# - Throw a warning if <50% of the particles ended up on the histogram


class SimPrad:
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
                cartesian (x,y,z) : (meters, meters, meters)
                cylindrical (r, theta, z) : (meters, radians, meters)
                spherical (r, theta, phi) : (meters, radians, radians)
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
        """
        Initalize the simPrad object, carry out coordinate transformations,
        and compute several quantities that will be used elsewhere.
        """

        self.grid = grid
        self.E = E
        self.B = B
        self.proton_energy = proton_energy
        self.verbose = verbose

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

    def _init_particle_sim(self):
        """"
        Auto-detects whether the given grid is uniform or irregular and
        creates an appropriate interpolator. Also calculates the grid
        resolution for use in choosing an appropriate timestep.
        """
        # Create a grid of indices for use in interpolation
        nx, ny, nz, x = self.grid.shape
        indgrid = np.indices([nx, ny, nz])
        indgrid = np.moveaxis(indgrid, 0, -1)

        # Auto-detect grid regularity
        variance = np.zeros([3])
        dx = np.gradient(self.grid[..., 0], axis=0)
        variance[0] = np.std(dx) / np.mean(dx)
        dy = np.gradient(self.grid[..., 1], axis=1)
        variance[1] = np.std(dy) / np.mean(dy)
        dz = np.gradient(self.grid[..., 2], axis=2)
        variance[2] = np.std(dz) / np.mean(dz)

        if np.allclose(variance, 0.0, atol=1e-6):
            self._log("Auto-detected a regularly spaced grid.")
            self.regular_grid = True
        else:
            self._log("Auto-detected an irregularly spaced grid.")
            self.regular_grid = False

        if self.regular_grid:
            # Create axes under the regular grid assumption
            xaxis = self.grid[:, 0, 0, 0]
            yaxis = self.grid[0, :, 0, 1]
            zaxis = self.grid[0, 0, :, 2]
            dx = np.mean(np.gradient(xaxis))
            dy = np.mean(np.gradient(yaxis))
            dz = np.mean(np.gradient(zaxis))
            dvec = np.array([dx.value, dy.value, dz.value]) * dx.unit

            # Estimate the grid-point spacing along the source_to_detector vector
            # Will be expanded to a constant array of length nparticles_grid
            # when _adaptive_ds() is called.
            self.ds = np.linalg.norm(np.dot(dvec, self.exp_ax))

            # Initialize the interpolator
            pts = (xaxis.si.value, yaxis.si.value, zaxis.si.value)
            self._log("Creating regular grid interpolator")
            self.interpolator = interp.RegularGridInterpolator(
                pts, indgrid, method="nearest", bounds_error=False, fill_value=-1
            )

        else:
            # Flat arrays of points for irregular grid interpolation fcn
            pts = np.zeros([nx * ny * nz, 3])
            pts[:, 0] = self.grid[:, :, :, 0].flatten().si.value
            pts[:, 1] = self.grid[:, :, :, 1].flatten().si.value
            pts[:, 2] = self.grid[:, :, :, 2].flatten().si.value

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
                self._log("Creating nearest-neighbor grid")
                self.nearest_neighbor = _nearest_neighbor(self.grid)

            # TODO
            # self.ds is used for determining when particles are on-grid
            # Somewhat ambiguous how to chose a single value for this for an
            # irregular grid: this may not be the best solution.
            self.ds = np.median(self.nearest_neighbor)

    def _max_theta_grid(self):
        """
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
                    vec = self.grid[x, y, z, :] - self.source

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
        """
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

        # Calculate where particles would hit the detector in the absence
        # of simulated  fields
        self._generate_null()

        # Create flags for tracking when particles during the simulation
        # on_grid -> zero if the particle is off grid, 1
        self.on_grid = np.zeros([self.nparticles_grid])
        # Entered grid -> non-zero if particle EVER entered the grid
        self.entered_grid = np.zeros([self.nparticles_grid])

    def _advance_to_grid(self):
        """
        Advances all particles to the timestep when the first particle should
        be entering the grid (to save time)
        """
        # Distance from the source to the nearest gridpoint
        dist = np.min(np.linalg.norm(self.grid - self.source, axis=3))

        # Time for fastest possible particle to reach the grid.
        t = (dist / self.v0).to(u.s)

        self.r = self.r + self.v * t

    def _advance_to_detector(self):
        """
        Once all particles have cleared the grid, advance them to the detector
        plane.

        This step applies to all particles, not just those that will hit the
        grid.
        """
        dist_remaining = np.dot(self.r, self.det_n) + np.linalg.norm(self.detector)

        v_towards_det = np.dot(self.v, -self.det_n)

        # If particles have not yet reached the detector plane and are moving
        # away from it, they will never reach the detector.
        condition = np.logical_and(v_towards_det < 0, dist_remaining > 0)
        ind = np.nonzero(np.where(condition, 0, 1))[0]
        self.r = self.r[ind, :]
        self.v = self.v[ind, :]
        self.nparticles_grid = self.r.shape[0]

        # Time remaining for each particle to reach detector plane
        t = dist_remaining / v_towards_det

        self.r += self.v * np.outer(t, np.ones(3))

        # Check that all points are now in the detector plane
        # (Eq. of a plane is nhat*x + d = 0)
        plane_eq = np.dot(self.r, self.det_n) + np.linalg.norm(self.detector)
        assert np.allclose(plane_eq, np.zeros(self.nparticles_grid), atol=1e-6)

    def _generate_null(self):
        """
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

    def _place_particles(self):
        """
        For each particle, find the indicies of the nearest field assuming that
        the fields are placed on a regular grid
        """
        # Interpolate the grid indices that best match each particle position
        i = self.interpolator(self.r[self.gi, :].si.value)

        # Store the grid positions
        self.xi = i[:, 0].astype(np.int32)
        self.yi = i[:, 1].astype(np.int32)
        self.zi = i[:, 2].astype(np.int32)

    def _adaptive_dt(self, B):
        """
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
        Bmag = np.linalg.norm(B, axis=1)  # B is [nparticles,3] here
        gyroperiod = (2 * np.pi * const.m_p.si / (const.e.si * np.max(Bmag))).to(u.s)

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
        """
        Compute the local grid resolution for each particle (for determining
        if the particle should be influenced by the grid or not). For regular
        grids this is a constant array, but for irregular grids it changes
        as the particles move
        """

        if self.regular_grid:
            # If self.ds is a scalar (as setup by the init function)
            # extend it to be the length of npartiles_grid
            if self.ds.size == 1:
                self.ds = self.ds * np.ones(self.nparticles_grid)
            else:
                pass
        else:
            # Calculate ds for each particle as the nearest-neighbor
            # distance of its assigned gridpoint
            self.ds = self.nearest_neighbor[self.xi, self.yi, self.zi]

    def _push(self):
        """
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
            self.r[self.gi, :] - self.grid[self.xi, self.yi, self.zi, :], axis=1
        )
        self.on_grid = np.where(dist < self.ds, 1, 0)
        self.entered_grid += self.on_grid

        # Retreive the fields
        E = self.E[self.xi, self.yi, self.zi, :] * np.outer(self.on_grid, np.ones(3))
        B = self.B[self.xi, self.yi, self.zi, :] * np.outer(self.on_grid, np.ones(3))

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
    ):
        """
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

        Returns
        -------
        None.

        """
        # Load inputs
        self.nparticles = int(nparticles)
        self.dt = dt
        self.dt_range = dt_range

        # Initialize the interpolators and grid resolution matrices
        self._init_particle_sim()

        # Initialize variables and create the particle distribution
        self._generate_particles(max_theta=max_theta)

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
        """
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
            return True
        else:
            return False

    def calc_ke(self, total=True):
        """
        Calculate the total kinetic energy of some or all particles. This calculation
        is currently done on velocity time steps (half-integer time steps)
        but it's good enough for ensuring energy is conserved.'
        """
        ke = 0.5 * self.mass * np.sum(self.v ** 2, axis=1)

        if total:
            return np.sum(ke).to(u.J)
        else:
            return ke.to(u.J)

    def synthetic_radiograph(self, size=None, bins=None):
        """
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
            The number of bins in each direction in the format
            [hbins, vbins]. The default is [250,250].

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

        # Determine locations of points in the detector plane using unit
        # vectors
        xloc = np.dot(self.r - self.detector, nx)
        yloc = np.dot(self.r - self.detector, ny)

        if size is None:
            # If a detector size is not given, choose lengths based on the
            # dimensions of the grid
            w = self.mag * np.max(
                [
                    np.max(np.abs(self.grid[..., 0].si.value)),
                    np.max(np.abs(self.grid[..., 1].si.value)),
                    np.max(np.abs(self.grid[..., 2].si.value)),
                ]
            )

            size = np.array([[-w, w], [-w, w]])
        else:
            size = size.si.value

        # If #bins is not set, make a guess
        if bins is None:
            bins = [200, 200]

        # Generate the histogram
        intensity, h, v = np.histogram2d(
            xloc.si.value, yloc.si.value, range=size, bins=bins
        )

        return h * u.m, v * u.m, intensity
