"""
Routines for the analysis of proton radiographs. These routines can be broadly
classified as either creating synthetic radiographs from prescribed fields or
methods of 'inverting' experimentally created radiographs to reconstruct the
original fields (under some set of assumptions).
"""

__all__ = [
    "AbstractField",
    "NoFields",
    "ElectrostaticGaussianSphere",
    "AxiallyMagnetizedCylinder",
    "ElectrostaticPlanarShock",
    "example_fields",
    "SyntheticProtonRadiograph",
]

import astropy.constants as const
import astropy.units as u
import numpy as np
import scipy.interpolate as interp
import warnings

from abc import ABC, abstractmethod
from scipy.special import erf as erf


class AbstractField(ABC):
    """
    Base class for example fields for testing SyntheticProtonRadiograph
    functions.
    """

    def __init__(self, grid, emax=0 * u.V / u.m, bmax=0 * u.T, regular_grid=None):
        r"""
        Initialize a field object

        Parameters
        ----------
        grid : `~astropy.units.Quantity` ndarray, (nx,ny,nz,3)
            Positions of grid points in space

        emax : `~astropy.units.Quantity`
            Maximum electric field to normalize example field to. The default
            is 0 * u.V / u.m.

        bmax : `~astropy.units.Quantity`
            Maximum magnetic field to normalize the example field to. The
            default is 0 * u.T.

        regular_grid : bool
            If None, the grid will be tested to determine whether or not it
            is regularly spaced. If True or False, the grid will be assumed
            to be regularly or irregularly spaced respectively without testing.

        Returns
        -------
        None.

        """
        self.grid = grid
        self._emax = emax
        self._bmax = bmax

        self.E = np.zeros(self.grid.shape) * u.V / u.m
        self.B = np.zeros(self.grid.shape) * u.T

        # If the regular grid flag is not set, auto detect whether the grid
        # is uniform
        if regular_grid is None:
            self.regular_grid = _detect_regular_grid(grid)
        else:
            self.regular_grid = regular_grid

        self.L = np.abs(np.max(grid) - np.min(grid))
        self.xarr = grid[..., 0]
        self.yarr = grid[..., 1]
        self.zarr = grid[..., 2]

        if self.regular_grid:
            self.xaxis = self.xarr[:, 0, 0]
            self.yaxis = self.yarr[0, :, 0]
            self.zaxis = self.zarr[0, 0, :]

        # Calculate radius arrays in spherical and cylindrical coordinates
        # for use in generating different field structures
        self.radius = np.sqrt(self.xarr ** 2 + self.yarr ** 2 + self.zarr ** 2)
        self.pradius = np.sqrt(self.xarr ** 2 + self.yarr ** 2)

        # Check that the model selected can be used (eg. no gradients on a
        # non-uniform grid)
        self._validate()

        # Generate the fields.
        self._gen_fields()

        # Normalize the fields to the emax and bmax values set
        self._norm()

    @abstractmethod
    def _validate(self):
        """
        Raise errors when the specified parameteters are not compatible with
        the field subclass that is being instantiated.
        """
        raise NotImplementedError

    @abstractmethod
    def _gen_fields(self):
        """
        Fill the empty E and B arrays with the example fields.
        """
        raise NotImplementedError

    def _norm(self):
        """
        Normalize the fields to the emax and bmax values provided.
        """
        max_E = np.max(self.E)
        max_B = np.max(self.B)

        if max_E != 0:
            self.E = self._emax * (self.E / max_E).to(u.dimensionless_unscaled)

        if max_B != 0:
            self.B = self._bmax * (self.B / max_B).to(u.dimensionless_unscaled)

    @property
    def bmax(self):
        """
        Get bmax
        """
        return self._bmax

    @bmax.setter
    def bmax(self, val):
        """
        Set bmax and re-normalize fields to the new value.
        """
        self._bmax = val
        # Normalize the fields to the new value
        self.norm()

    @property
    def emax(self):
        """
        Get emax
        """
        return self._emax

    @emax.setter
    def emax(self, val):
        """
        Set emax and re-normalize fields to the new value.
        """
        self._emax = val
        # Normalize the fields to the new values
        self.norm()


class NoFields(AbstractField):
    r"""
    Empty E and B arrays.
    """

    def _gen_fields(self):
        pass

    def _validate(self):
        pass



class ElectrostaticGaussianSphere(AbstractField):
    r"""
    A radial, spherically symmetric electric field produced by a spherical
    blob of potential with a Gaussian radial distribution:

    .. math::
        \phi = e^{-(r/a)^2}
    .. math::
        E = -\nabla \phi  (r < L/2)
    .. math::
        E = 0 (r < L/2)
    .. math::
        B = 0

    Where :math:`r` is the radius, :math:`L` is the field grid length scale,
    and :math:`a=L/3`.

    Since a gradient is used in the calculation,
    this option is not compatible with a non-uniform grid.
    """

    def _gen_fields(self):
        a = self.L / 3
        potential = np.exp(-(self.radius ** 2) / a ** 2) * u.V

        Ex, Ey, Ez = np.gradient(potential, self.xaxis, self.yaxis, self.zaxis)

        self.E[:, :, :, 0] = -1 * np.where(self.radius < 0.5 * self.L, Ex, 0)
        self.E[:, :, :, 1] = -1 * np.where(self.radius < 0.5 * self.L, Ey, 0)
        self.E[:, :, :, 2] = -1 * np.where(self.radius < 0.5 * self.L, Ez, 0)

    def _validate(self):
        if not self.regular_grid:
            raise ValueError(
                "ElectrostaticGaussianSphere can only be created "
                "on a regularly spaced grid."
            )


class AxiallyMagnetizedCylinder(AbstractField):
    r"""
    A cylinder of constant magnetic field aligned with the z-axis:

    .. math::
        E = 0
    .. math::
        B = 1 (\rho < a)
    .. math::
        B = 0 (\rho > a)

    Where :math:`\rho` is the cylinder radius, :math:`L` is the field grid
    length scale, and :math:`a=L/4`.
    """

    def _gen_fields(self):
        a = self.L / 4
        self.B[:, :, :, 2] = np.where(self.pradius < a, 400 * u.T, 0 * u.T)

    def _validate(self):
        pass


class ElectrostaticPlanarShock(AbstractField):
    r"""
    A model of an electrostatic planar shock. The discontinuity is located
    at z=0 and has a Gaussian distribution in the xy plane:

    .. math::
        \phi = (1 - \Gamma(z/\delta))e^(-(\rho/a)^2)
    .. math::
        E = - \nabla \phi
    .. math::
        B = 0

    Where :math:`\rho` is the cylindrical radius, :math:`a=Max(\rho)`,
    :math:`\delta=a/120 is the shock width`,
    and :math:`\Gamma()` is the Gauss error function.

    Since a gradient is used in the calculation,
    this option is not compatible with a non-uniform grid.
    """

    def _gen_fields(self):
        a = np.max(self.pradius) / 2
        delta = a / 120

        potential = (
            (1 - erf(self.zarr / delta)) * np.exp(-((self.pradius / a) ** 2)) * u.V
        )

        Ex, Ey, Ez = np.gradient(potential, self.xaxis, self.yaxis, self.zaxis)
        self.E[..., 0] = -Ex
        self.E[..., 1] = -Ey
        self.E[..., 2] = -Ez

    def _validate(self):
        if not self.regular_grid:
            raise ValueError(
                "ElectrostaticPlanarShock can only be created "
                "on a regularly spaced grid."
            )


def example_fields(
    grid=None,
    model="electrostatic gaussian sphere",
    regular_grid=True,
    num=(100, 100, 100),
    length=1 * u.mm,
    emax=1e9 * u.V / u.m,
    bmax=100 * u.T,
):
    r"""
    This function generates example fields based on analytical models for
    testing or demonstrating the proton radiography module.

    Parameters
    ----------
    grid : `~astropy.units.Quantity` array shape (nx,ny,nz,3), optional
        A grid of positions on which to calculate fields. If no grid is
        specified, one will be created.

    model : str, optional
        The field model to load. The default represents the electrostatic
        E-field created by a spherical Gaussian potential. The full list of
        options is:

        * "no fields": A grid with empty field arrays (E=B=0)

        * "electrostatic gaussian sphere": A radial, spherically symmetric
            electric field sphere produced by a sphere of potential with a
            Gaussian distribution. Since a gradient is used in the calculation,
            this option is not compatible with a non-uniform grid.

        * "axial magnetic field": A cylinder of magnetic field oriented along
            the z-axis, with E=0.

        * "electrostatic planar shock": A model of an electrostatic planar
            shock. The discontinuity is located at z=0 and has a Gaussian
            distribution in the xy plane.

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
            * If a single value is given, then the length of each dimension
                will be set to (-length, length)

            * If an array of shape (3) is given, each dimension is set to be
                symmetric using each value, eg. xdim = [-length[0], length[0]]...

            * If an array of shape (2,3) is given, then L[:,i] is the min
                and max of the ith dimension.

    emax: `~astropy.units.Quantity`, optional
        Scale E-field to this maximum value. Default is 1e9 V/m

    bmax: `~astropy.units.Quantity`, optional
        Scale B-field to this maximum value. Default is 100 T

    Returns
    -------
    grid : `~astropy.units.Quantity` array, shape (nx,ny,nz,3)
        The grid of positions corresponding to the field grids.

    E : `~astropy.units.Quantity` array, shape (nx,ny,nz,3)
        Electric field array in units of V/m

    B : `~astropy.units.Quantity` array, shape (nx,ny,nz,3)
        Magnetic field array in units of Tesla

    """

    # If no grid is specified, create a grid
    if grid is None:
        grid = _make_grid(num, length, regular_grid=regular_grid)

    # Load the model class for the test example chosen
    models = {
        "no fields": NoFields,
        "electrostatic gaussian sphere": ElectrostaticGaussianSphere,
        "electrostatic planar shock": ElectrostaticPlanarShock,
        "axial magnetic field": AxiallyMagnetizedCylinder,
    }

    # Generate the fields by instantiating the test field object
    fields = models[model](grid, emax=emax, bmax=bmax, regular_grid=regular_grid)

    return fields.grid, fields.E, fields.B


def _make_grid(num, length, regular_grid=True):
    r"""
    Creates a grid given a number of options

    Parameters
    ----------

    num : int or tuple of 3 ints, optional
        The number of points along each axis

    length : ndarray of quantities
        The length of each dimension, which can be given in several ways:

            * If a single value is given, then the length of each dimension
            will be set to (-length, length)

            * If an array of shape (3) is given, each dimension is set to be
                symmetric using each value, eg. xdim = [-length[0], length[0]]...

            * If an array of shape (2,3) is given, then L[:,i] is the min
                and max of the ith dimension.

    regular_grid : bool
        If True, generate a regularly spaced grid. If False, generate a grid
        with irregular spacing. Default is True.

    Returns
    -------
    grid : ndarray (nx,ny,nz,3)
    """

    # Process different options for inputs
    if isinstance(num, int):
        num = (num, num, num)

    L = np.zeros([2, 3]) * length.unit
    if length.size == 1:
        L[:, 0] = np.array([-length.value, length.value]) * length.unit
        L[:, 1] = np.array([-length.value, length.value]) * length.unit
        L[:, 2] = np.array([-length.value, length.value]) * length.unit
    elif length.size == 3:
        L[:, 0] = np.array([-length[0].value, length[0].value]) * length.unit
        L[:, 1] = np.array([-length[1].value, length[1].value]) * length.unit
        L[:, 2] = np.array([-length[2].value, length[2].value]) * length.unit
    else:
        raise ValueError("Dimensions of length in _create_grid are not valid.")

    # Create a blank grid
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

    return grid


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


def _nearest_neighbor(array):
    r"""
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


def _vector_intersects_grid(p1, p2, grid):
    r"""
    Returns True if the vector from p1 to p2 intersects the grid. Otherwise,
    returns false. Assumes the grid is at least approximately axis-aligned.
    This is a standard ray-box intersection algorithm.
    """
    p1, p2, grid = p1.si.value, p2.si.value, grid.si.value
    # Caclulate the minimum and maximum of each
    Ax, Bx = np.min(grid[..., 0]), np.max(grid[..., 0])
    Ay, By = np.min(grid[..., 1]), np.max(grid[..., 1])
    Az, Bz = np.min(grid[..., 2]), np.max(grid[..., 2])
    A = np.array([Ax, Ay, Az])
    B = np.array([Bx, By, Bz])

    # Calculate the equation of the line from p1 to p2 such that
    # r = p1 + t*D
    D = p2 - p1

    # Calculate the intersection points. These operations are just vectorized
    # for convenience. Ignore div-by-zero: outputting infty's here is fine.
    with np.errstate(divide="ignore"):
        Tmin = (A - p1) / D
        Tmax = (B - p1) / D

    Tmin = np.max(Tmin)
    Tmax = np.min(Tmax)

    return Tmin < Tmax


def _detect_regular_grid(grid, tol=1e-6):
    """
    Determine whether a grid is regular (uniformly spaced) by computing the
    variance of the grid gradients.
    """
    variance = np.zeros([3])
    dx = np.gradient(grid[..., 0], axis=0)
    variance[0] = np.std(dx) / np.mean(dx)
    dy = np.gradient(grid[..., 1], axis=1)
    variance[1] = np.std(dy) / np.mean(dy)
    dz = np.gradient(grid[..., 2], axis=2)
    variance[2] = np.std(dz) / np.mean(dz)

    return np.allclose(variance, 0.0, atol=tol)


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

        # Validate input arrays
        arr = {"grid": grid, "E": E, "B": B}
        for x in arr.keys():
            if not np.isfinite(arr[x]).all():
                raise ValueError(
                    f"Input arrays must be finite: {x} contains "
                    "either NaN or infinite values."
                )

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

        # Check that source-detector vector actually passes through the grid
        if not _vector_intersects_grid(self.source, self.detector, self.grid):
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

        self.regular_grid = _detect_regular_grid(self.grid)
        if self.regular_grid:
            self._log("Auto-detected a regularly spaced grid.")
        else:
            self._log("Auto-detected an irregularly spaced grid.")

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
        dist = np.min(np.linalg.norm(self.grid - self.source, axis=3))

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

        If the vol_weighted_fields flag is set, the fields returned are an
        average of the fields at the gridpoints surrounding the
        interpolated particle grid point, weighted by the distance of each grid
        point from the particle position.

        If no flags are set, the E and B fields are returned at the
        interpolated grid point.

        At the end, particles that are off the grid (as determined by the
        on_grid array) have their fields set to zero.
        """

        if self.vol_weighted_fields:
            # Initialize arrays
            l = len(self.xi)
            dist = np.zeros([27,l])
            E_est = np.zeros([27, l,3])*self.E.unit
            B_est = np.zeros([27, l,3])*self.B.unit
            E = np.zeros([l,3])*self.E.unit
            B = np.zeros([l,3])*self.B.unit
            nx,ny,nz,nax = self.grid.shape

            # Loop through all of the points surrounding (and including) the
            # interpolated gridpoint
            ind = 0
            for dx in np.array([-1, 0, 1]).astype(np.int32):
                for dy in np.array([-1, 0, 1]).astype(np.int32):
                    for dz in np.array([-1, 0, 1]).astype(np.int32):

                        x = self.xi + dx
                        y = self.yi + dy
                        z = self.zi + dz

                        # Figure out if any indices are out-of-bounds
                        valid = ((x>=0) & (x<nx)
                                 & (y>=0) & (y<ny)
                                 & (z>=0) & (z<nz))
                        out = np.where(valid == False)

                        # Set out-of-bound indices to something inbounds to
                        # prevent errors.
                        x[out], y[out], z[out] = 0,0,0

                        dist [ind,:] = np.linalg.norm(self.grid[x,y,z,:] -
                                                      self.r[self.gi, :],
                                                      axis=1)

                        # Set out-of-bound distances to infinity, since
                        # no fields are to be found there.
                        dist[ind, out] = np.infty

                        # Record the field at each point
                        E_est[ind,:,:] = self.E[x, y, z, :]
                        B_est[ind,:,:] = self.B[x, y, z, :]

                        # Increment the index counter
                        ind += 1

            # Weight by distance
            weights = 1/dist

            # Average the fields according to the weights.
            # TODO: Find a more elegant way to write this? Problem is that
            # dimensions of weights is not the same as E,B
            E[:,0] = np.average(E_est[:,:,0], weights=weights, axis=0)
            E[:,1] = np.average(E_est[:,:,1], weights=weights, axis=0)
            E[:,2] = np.average(E_est[:,:,2], weights=weights, axis=0)

            B[:,0] = np.average(B_est[:,:,0], weights=weights, axis=0)
            B[:,1] = np.average(B_est[:,:,1], weights=weights, axis=0)
            B[:,2] = np.average(B_est[:,:,2], weights=weights, axis=0)

        else:
            E = self.E[self.xi, self.yi, self.zi, :]
            B = self.B[self.xi, self.yi, self.zi, :]


        # Set fields for off-grid particles
        E = E * np.outer(self.on_grid, np.ones(3))
        B = B * np.outer(self.on_grid, np.ones(3))

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
            self.r[self.gi, :] - self.grid[self.xi, self.yi, self.zi, :], axis=1
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
        vol_weighted_fields = False,
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

        Returns
        -------
        None.

        """
        # Load inputs
        self.nparticles = int(nparticles)
        self.dt = dt
        self.dt_range = dt_range
        self.vol_weighted_fields = vol_weighted_fields

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
                    np.max(np.abs(self.grid[..., 0].value)),
                    np.max(np.abs(self.grid[..., 1].value)),
                    np.max(np.abs(self.grid[..., 2].value)),
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
