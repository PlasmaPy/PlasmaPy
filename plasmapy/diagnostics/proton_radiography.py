"""
Routines for the analysis of proton radiographs. These routines can be broadly
classified as either creating synthetic radiographs from prescribed fields or
methods of 'inverting' experimentally created radiographs to reconstruct the
original fields (under some set of assumptions).
"""

__all__ = [
    "test_fields",
    "SimPrad",
    "calc_ke",
    "run"
    "synthetic_radiograph",
]

import astropy.constants as const
import astropy.units as u

import numpy as np

import time

import scipy.interpolate as interp


def test_fields(grid=None, mode='electrostatic gaussian sphere'):

    # Create or load a grid
    if grid is None:
        L = 2.5*u.mm
        xaxis = np.linspace(-L, L, num=250)
        yaxis = np.linspace(-L, L, num=250)
        zaxis = np.linspace(-L, L, num=250)
        xarr, yarr, zarr = np.meshgrid(xaxis, yaxis, zaxis, indexing='ij')
        grid = np.zeros([xaxis.size, yaxis.size, zaxis.size, 3])*u.cm
        grid[...,0] = xarr
        grid[...,1] = yarr
        grid[...,2] = zarr

    else:
        xaxis = grid[:,0,0,0]
        yaxis = grid[0,:,0,1]
        zaxis = grid[0,0,:,2]

        xarr = grid[...,0]
        yarr = grid[...,1]
        zarr = grid[...,2]

        L = np.max(xaxis) - np.min(xaxis)

    # Create arrays for fields
    E = np.zeros(grid.shape)*u.V/u.m
    B = np.zeros(grid.shape)*u.T


    if mode == 'no fields':
        pass

    # Electrostatic cylinder creates a radial electric field with the axis
    # pointing along the z axis
    elif mode == 'electrostatic gaussian sphere':
        print("Generating Electrostatic Gaussian Sphere")

        alen = L/4
        blen = L/2
        potential = -np.exp(-(xarr**2 + yarr**2)/alen**2 -
                               zarr**2/blen**2)

        Ex, Ey, Ez = np.gradient(potential, xaxis, yaxis, zaxis)
        E[...,0], E[...,1], E[...,2] = Ex, Ey, Ez

        # Normalize to a desired maximum field
        E = E/np.max(E)
        E = E*1e10*u.V/u.m

    return grid, E, B


def _rot_a_to_b(a,b):
    """
    Calculates the 3D rotation matrix that will rotate vector a to be aligned
    with vector b.

    Algorithm based on this discussion:
    https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
    """
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)

    # TODO: manually handle the case where b = -zhat: that breaks this fcn
    axb = np.cross(a,b)
    c = np.dot(a,b)
    vskew = np.array([[0, -axb[2], axb[1]],
             [axb[2], 0, -axb[0]],
             [-axb[1], axb[0], 0]]).T # Transpose to get right orientation

    return np.identity(3) + vskew + np.dot(vskew, vskew)/(1 + c)



class SimPrad():
    """

    source -> 3-element ndarray vector from the origin to the source
    of the protons.

    detector -> 3-element ndarray specifying center of the collection plane
    (the point an undeflected proton would take straight from the source).

    geometry -> Choose cartesian, spherical, or cylindrical coordinates for the
    source and detector locations (grid is ALWAYS in cartesian and the
    detector plane is also always flat).

    """

    def __init__(self, grid: u.m,
                 E: u.V/u.m,
                 B: u.T,
                 source: u.m,
                 detector: u.m,
                 proton_energy = 14*u.MeV,
                 geometry = 'cartesian',
                 solid = None,
                 verbose = True):

        self.verbose = verbose
        self.grid = grid
        self.E = E
        self.B = B
        self.proton_energy = proton_energy

        self.charge = const.e.si
        self.mass = const.m_p.si
        # Calculate the velocity corresponding to the proton energy
        self.v0 = np.sqrt(2*self.proton_energy/const.m_p.si).to(u.m/u.s)

        # Extract axes from the given grid
        self.xaxis = self.grid[:,0,0,0]
        self.yaxis = self.grid[0,:,0,1]
        self.zaxis = self.grid[0,0,:,2]
        self.dx = np.mean(np.gradient(self.xaxis))
        self.dy = np.mean(np.gradient(self.yaxis))
        self.dz = np.mean(np.gradient(self.zaxis))

        # Compute maximum separation between gridpoints
        # this is used for calculating the timestep
        self.ds = np.sqrt(3/2)*max(self.dx, self.dy, self.dz)

        # Convert geometrical inputs between coordinates systems
        if geometry == 'cartesian':
            self.source = source
            self.detector = detector

        elif geometry == 'spherical':
            r, t, p = source
            t = np.deg2rad(t)
            p = np.deg2rad(p)
            self.source = np.zeros(3)*r.unit
            self.source[0] = r*np.sin(t)*np.cos(p)
            self.source[1] = r*np.sin(t)*np.sin(p)
            self.source[2] = r*np.cos(t)

            r, t, p = detector
            t = np.deg2rad(t)
            p = np.deg2rad(p)
            self.detector = np.zeros(3)*r.unit
            self.detector[0] = r*np.sin(t)*np.cos(p)
            self.detector[1] = r*np.sin(t)*np.sin(p)
            self.detector[2] = r*np.cos(t)

        # TODO: Implement cylindrical coordinates for completeness

        # Calculate some parameters involving the source and detector locations
        # Detector plane normal vector
        self.det_n = self.detector/np.linalg.norm(self.detector)

        # Vector directly from source to detector
        self.source_to_detector = self.detector - self.source

        # Compute the magnification
        self.mag = 1 + (np.linalg.norm(self.detector)
                        /np.linalg.norm(self.source))

    def _log(self, msg):
        if self.verbose:
            print(msg)

    #*************************************************************************
    # Simulated Particle Tracing Code (Psim) functions
    #*************************************************************************

    def _init_particle_sim(self):
        """"
        Initialize variables necessary for running the particle prad sim
        """

        # This grid of indices is used when interpolating particles onto the
        # fields. Generated here to avoid the cost of re-generating it with
        # each push
        self.ind_grid = np.indices([self.xaxis.size,
                                   self.yaxis.size,
                                   self.zaxis.size])
        self.ind_grid = np.moveaxis(self.ind_grid, 0, -1)

        # These comparison axes will be used for finding the best index matches
        # when placing particles on the field grid
        self.xaxis_compare = np.outer(self.xaxis, np.ones(self.nparticles))
        self.yaxis_compare = np.outer(self.yaxis, np.ones(self.nparticles))
        self.zaxis_compare = np.outer(self.zaxis, np.ones(self.nparticles))


    def _generate_particles(self):
        """
        Generates the angular distributions about the Z-axis, then
        shifts those distributions to align with the source-to-origin axis
        """

        self._log("Creating Particles")
        # Calculate angle range for particles to possibly hit the sim volume
        # there's no point including particles at larger angles, which will
        # never interact with the simulated fields.
        max_angle = 2.5(*np.arctan(np.max(self.grid)
                                   /np.linalg.norm(self.source)).value)

        # Create a probability vector
        arg = np.linspace(0, max_angle, num=int(1e4))
        prob = np.sin(arg)
        prob *= 1/np.sum(prob)

        # Randomly choose theta's weighted with the sine probabilities
        theta = np.random.choice(arg, size=self.nparticles,
                                 replace=True, p=prob)

        # Uniform Phi distribution
        phi = np.random.uniform(size=self.nparticles)*2*np.pi

        # Construct the velocity distribution around the z-axis
        self.v = np.zeros([self.nparticles, 3]) * u.m/u.s
        self.v[:,0] = self.v0*np.sin(theta)*np.cos(phi)
        self.v[:,1] = self.v0*np.sin(theta)*np.sin(phi)
        self.v[:,2] = self.v0*np.cos(theta)

        # Calculate the rotation matrix
        a = np.array([0,0,1])
        b = self.detector - self.source
        rot = _rot_a_to_b(a,b)

        # Apply rotation matrix to calculated velocity distribution
        self.v = np.matmul(self.v, rot)

        # Place particles at the source
        self.r = np.outer(np.ones(self.nparticles), self.source)

        # Calculate where particles would hit the detector in the absence
        # of simulated  fields
        self._generate_null()

        # Create flags for tracking when particles during the simulation
        # on_grid -> zero if the particle is off grid, 1
        self.on_grid = np.zeros([self.nparticles])
        # Entered grid -> non-zero if particle EVER entered the grid
        self.entered_grid = np.zeros([self.nparticles])


    def _advance_to_grid(self):
        """
        Advances all particles to the timestep when the first particle should
        be entering the grid (to save time)
        """
        # Distance from the source to the nearest gridpoint
        dist = np.min(np.linalg.norm(self.grid - self.source, axis=3))

        # Source-to-origin center unit vector
        uvec = -self.source/ np.linalg.norm(self.source)

        # Time for fastest possible particle to reach the grid.
        t = (dist*uvec/self.v0).to(u.s)

        self.r = self.r + self.v*t

    def _advance_to_detector(self):
        """
        Once all particles have cleared the grid, advance them to the detector
        plane.
        """
        # Calculate the remaining distance each particle needs to travel to
        # reach the plane
        # see https://mathworld.wolfram.com/Point-PlaneDistance.html
        remaining = np.dot(self.r - self.detector, self.det_n)

        # Component of the particle velocity in the direction of the detector
        # plane
        v_in_dir = np.dot(self.v, self.det_n)

        # Calculate the time remaining to reach that plane and push
        # t can be negative if particle has already passed the detector plane
        t = (remaining/v_in_dir).to(u.s)

        # TODO: Negative times should only really be allowed here if particles
        # have passed the detector plane (otherwise they correspond) to
        # particles going the wrong way. Maybe solve by checking that
        # v dot uvec is positive?

        self.r = self.r + self.v*np.outer(t,np.ones(3))

    def _generate_null(self):
        """
        Calculate the distribution of particles on the detector in the absence
        of any simulated fields.
        """
        # Calculate the unit vector from the source to the detector
        dist = np.linalg.norm(self.source_to_detector)
        uvec = self.source_to_detector.to(u.m).value/dist.to(u.m).value

        # Calculate the remaining distance each particle needs to travel
        # along that unit vector
        remaining = np.dot(self.source, uvec)

        # Calculate the time remaining to reach that plane and push
        t = ((dist - remaining)/np.dot(self.v, uvec)).to(u.s)

        # Calculate the particle positions for that case
        self.r0 = self.source + self.v*np.outer(t,np.ones(3))


    def _grid_particles(self):
        """
        For each particle, find the indicies of the nearest field grid point
        """
        # Interpolate the grid indices that best match each particle position
        pts = (self.xaxis.si.value, self.yaxis.si.value, self.zaxis.si.value)
        ind_interp = interp.RegularGridInterpolator(pts,
                                           self.ind_grid, method='nearest',
                                           bounds_error = False,
                                           fill_value = -1)
        i = ind_interp(self.r.si.value)

        # Store the grid positions
        self.xi = i[:,0].astype(np.int32)
        self.yi = i[:,1].astype(np.int32)
        self.zi = i[:,2].astype(np.int32)

        # Update the list of particles on and off the grid
        dist = np.linalg.norm(self.r - grid[self.xi,self.yi,self.zi, :], axis=1)
        self.on_grid = np.where(dist < self.ds, 1, 0)
        self.entered_grid += self.on_grid


    def _adaptive_dt(self, B):
        """
        Calculate the appropraite dt based on a number of considerations
        """
        # Compute the cyclotron gyroperiod
        gyroperiod = (2*np.pi*const.m_p.si/(const.e.si*np.max(np.abs(B)))).to(u.s)

        # Compute the time required for a proton to cross one grid cell
        gridstep = (min(self.dx, self.dy, self.dz)/self.v0).to(u.s)

        return min(gyroperiod, gridstep)*self.dt_mult


    def _push(self):
        """
        Advance particles using an implementation of the time-centered
        Boris algorithm
        """
        # Calculate the indices of the field grid points nearest to each particle
        self._grid_particles()

        E = self.E[self.xi, self.yi, self.zi, :]*np.outer(self.on_grid, np.ones(3))
        B = self.B[self.xi, self.yi, self.zi, :]*np.outer(self.on_grid, np.ones(3))

        # Calculate the adaptive timestep from the fields currently experienced
        # by the particles
        dt = self._adaptive_dt(B)
        dt2 = dt*self.charge/self.mass/2

        # Execute the Boris push algorithm
        vminus = self.v + E*dt2
        t = -B*dt2
        s = 2 * t / (1 + (t * t).sum(axis=1, keepdims=True))
        vprime = vminus + np.cross(vminus.si.value, t) * u.m/u.s
        vplus = vminus + np.cross(vprime.si.value, s) * u.m/u.s
        vnew = vplus + E*dt2
        self.v = vnew

        self.dr = dt*self.v

        self.r += self.dr


    def run(self, nparticles=1e5, dt_mult=1):
        """
        Setup and run a particle-tracing simulated radiograph
        """
        # Load inputs
        self.nparticles = int(nparticles)
        self.dt_mult = dt_mult

        # Initialize some particle-sim-specific variables
        self._init_particle_sim()

        # Initialize variables and create the particle distribution
        self._generate_particles()

        # Advance the particles to the near the start of the simulation
        # volume
        self._advance_to_grid()

        # Push the particles until the stop condition is satisfied
        # (no more particles on the simulation grid)
        while not self._stop_condition():
            if self.verbose:
                fract = 100*np.sum(self.on_grid)/self.nparticles
                self._log(f"{fract:.1f}% on grid")
            self._push()

        # Advance the particles to the image plane
        self._advance_to_detector()


    def _stop_condition(self):
        """
        The stop condition is that most of the particles have entered the grid
        and almost all have now left it.
        """
        return (np.sum(self.entered_grid)/self.nparticles > .5 and
                np.sum(self.on_grid)/self.nparticles < 0.01)


    def calc_ke(self, total=True):
        """
        Calculate the total kinetic energy of some or all particles. This calculation
        is currently done on velocity time steps (half-integer time steps)
        but it's good enough for ensuring energy is conserved.'
        """
        ke = 0.5*self.mass*np.sum(self.v**2, axis=1)

        if total:
            return np.sum(ke).to(u.J)
        else:
            return ke.to(u.J)

    def synthetic_radiograph(self, size=None, bins=None):
        """
        Calculate a "synthetic radiograph" (particle count histogram in the
        image plane)
        """
        # Note that, at the end of the simulation, all particles were moved
        # into the image plane.

        # Define detector x axis as being perpendicular to both the detector
        # vector and the grid z-axis
        nx = np.cross(np.array([0,0,1]), self.det_n)
        ny = np.cross(nx, self.det_n)

        # Determine locations of points in the detector plane using unit
        # vectors
        xloc = np.dot(self.r, nx)
        yloc = np.dot(self.r, ny)

        if size is None:
            # If a detector size is not given, choose lengths based on the
            # dimensions of the grid
            w = self.mag*np.max([ np.max(np.abs(self.xaxis.si.value)),
                        np.max(np.abs(self.yaxis.si.value)),
                        np.max(np.abs(self.zaxis.si.value)) ])

            size = np.array([[-w, w], [-w,w]])

        # If #bins is not set, make a guess
        if bins is None:
            bins = [100, 100]

        # Generate the histogram
        intensity, h, v = np.histogram2d(xloc.si.value, yloc.si.value,
                                         range=size, bins=bins)

        return h*u.m, v*u.m, intensity

