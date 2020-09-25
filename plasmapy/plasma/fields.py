"""
Defines the PosGrid and AbstractField classes and provides subclasses
representing different field configurations
"""

__all__ = [
    "PosGrid",
    "AbstractField",
    "NoFields",
    "ElectrostaticGaussianSphere",
    "AxiallyMagnetizedCylinder",
    "ElectrostaticPlanarShock",
    "example_fields",
]

import astropy.constants as const
import astropy.units as u
import numpy as np

from abc import ABC, abstractmethod
from scipy.special import erf as erf


class PosGrid:
    """
    Represents a 3D position grid. If a grid is passed in through the grid
    keyword, the object is directily initialized from this grid and other
    keywords are ignored. Otherwise, a grid is constructed based on the
    keywords described below.

    Parameters
    ----------
    grid : `~astropy.units.Quantity`, shape (nx,ny,nz,3)
        An array giving the positions of each grid point. Units must be
        convertable to meters. If this keyword is set, this grid is used
        to initialize the object and all other keywords are ignored.

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

    """

    def __init__(self, **kwargs):

        # If grid has been provided, initialize object with the provided
        # grid
        if "grid" in kwargs.keys():
            self.grid = kwargs["grid"]

            if not isinstance(self.grid, u.Quantity):
                raise ValueError(
                    "Grid must be an astropy.units.Quantity, but "
                    f"type given was {type(self.grid)}"
                )

            if not len(self.grid.shape) == 4 or self.grid.shape[3] != 3:
                raise ValueError(
                    "Grid must have shape [nx,ny,nz,3], but "
                    f"grid given has shape {self.grid.shape}"
                )

        # Otherwise create a grid
        else:
            self._make_grid(**kwargs)

        # If regular grid is set, use that value
        # otherwise set to None and grid regularity will be auto-detected
        # if needed.
        if "regular_grid" in kwargs.keys():
            self._regular_grid = kwargs["regular_grid"]
        else:
            self._regular_grid = None

        # Properties to be calculated as needed
        self._nearest_neighbor = None

    @property
    def regular_grid(self):
        """
        Get value of regular_grid
        If None, calculate
        """
        if self._regular_grid is None:
            self._detect_regular_grid()
        return self._regular_grid

    @property
    def nearest_neighbor(self):
        """
        Get value of nearest_neighbor
        If None, calculate
        """
        if self._nearest_neighbor is None:
            self._calculate_nearest_neighbor()
        return self._nearest_neighbor

    @property
    def xarr(self):
        """
        Returns the 3D array of x values
        """
        return self.grid[..., 0]

    @property
    def yarr(self):
        """
        Returns the 3D array of y values
        """
        return self.grid[..., 1]

    @property
    def zarr(self):
        """
        Returns the 3D array of z values
        """
        return self.grid[..., 2]

    @property
    def xaxis(self):
        """
        Get x-axis (only valid for a uniform grid)
        """
        if self._regular_grid is True:
            return self.grid[:, 0, 0, 0]
        else:
            raise ValueError(
                "xaxis property is only well-defined for " "uniformly-spaced grids."
            )

    @property
    def yaxis(self):
        """
        Get y-axis (only valid for a uniform grid)
        """
        if self._regular_grid is True:
            return self.grid[0, :, 0, 1]
        else:
            raise ValueError(
                "yaxis property is only well-defined for " "uniformly-spaced grids."
            )

    @property
    def zaxis(self):
        """
        Get z-axis (only valid for a uniform grid)
        """
        if self._regular_grid is True:
            return self.grid[0, 0, :, 2]
        else:
            raise ValueError(
                "zaxis property is only well-defined for " "uniformly-spaced grids."
            )

    @property
    def dx(self):
        """
        Calculate dx (only valid for a uniform grid)
        """
        if self._regular_grid is True:
            return np.mean(np.gradient(self.xaxis))
        else:
            raise ValueError(
                "dx property is only well-defined for " "uniformly-spaced grids."
            )

    @property
    def dy(self):
        """
        Calculate dy (only valid for a uniform grid)
        """
        if self._regular_grid is True:
            return np.mean(np.gradient(self.yaxis))
        else:
            raise ValueError(
                "dy property is only well-defined for " "uniformly-spaced grids."
            )

    @property
    def dz(self):
        """
        Calculate dz (only valid for a uniform grid)
        """
        if self._regular_grid is True:
            return np.mean(np.gradient(self.zaxis))
        else:
            raise ValueError(
                "dz property is only well-defined for " "uniformly-spaced grids."
            )

    @property
    def shape(self):
        """
        Shortcut method for retriving the shape of the grid
        """
        return self.grid.shape

    @property
    def unit(self):
        """
        Shortcut method for retriving the unit of the grid
        """
        return self.grid.unit

    def _make_grid(self, **kwargs):
        """
        Create a grid based on keywords set
        """

        if "num" in kwargs.keys():
            num = kwargs["num"]
        else:
            num = 100

        if "length" in kwargs.keys():
            length = kwargs["length"]
        else:
            length = 1 * u.cm

        if "regular_grid" in kwargs.keys():
            regular_grid = kwargs["regular_grid"]
        else:
            regular_grid = True

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

        self._regular_grid = regular_grid
        self.grid = grid

    def _detect_regular_grid(self, tol=1e-6):
        """
        Determine whether a grid is regular (uniformly spaced) by computing the
        variance of the grid gradients.
        """
        variance = np.zeros([3])
        dx = np.gradient(self.grid[..., 0], axis=0)
        variance[0] = np.std(dx) / np.mean(dx)
        dy = np.gradient(self.grid[..., 1], axis=1)
        variance[1] = np.std(dy) / np.mean(dy)
        dz = np.gradient(self.grid[..., 2], axis=2)
        variance[2] = np.std(dz) / np.mean(dz)

        self._regular_grid = np.allclose(variance, 0.0, atol=tol)

    def _calculate_nearest_neighbor(self):
        r"""
        Creates a Quantity array of the distance to the nearest neighbor from each
        point of shape [nx,ny,nz]. This output can then be used to define the local
        scalar grid resolution on an irregular grid.
        """
        nx, ny, nz, x = self.grid.shape

        dist = np.zeros([nx, ny, nz, 26]) * self.grid.unit

        ind = 0
        for x in [-1, 0, 1]:
            for y in [-1, 0, 1]:
                for z in [-1, 0, 1]:

                    # Skip the zero shift case.
                    if x == 0 and y == 0 and z == 0:
                        continue

                    # Shift the array
                    shifted_grid = np.roll(self.grid, shift=(x, y, z), axis=(0, 1, 2))
                    # Compute the distance between the shifted array points and the
                    # previous array points
                    dist[..., ind] = np.linalg.norm(self.grid - shifted_grid, axis=3)
                    # Increment counter
                    ind += 1

        # Return the minimum distance between all 27 points
        self._nearest_neighbor = np.min(dist, axis=3)

    def vector_intersects(self, p1, p2):
        r"""
        Returns True if the vector from p1 to p2 intersects the grid. Otherwise,
        returns false. Assumes the grid is at least approximately axis-aligned.
        This is a standard ray-box intersection algorithm.
        """
        p1, p2, grid = p1.si.value, p2.si.value, self.grid.si.value
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


class AbstractField(ABC):
    """
    Base class for field grids.
    """

    def __init__(self, grid, emax=0 * u.V / u.m, bmax=0 * u.T):
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

        Returns
        -------
        None.

        """

        # If grid is not already a PosGrid object, try making it one
        if not isinstance(grid, PosGrid):
            self.grid = PosGrid(grid=grid)
        else:
            self.grid = grid

        self._emax = emax
        self._bmax = bmax

        self.E = np.zeros(self.grid.shape) * u.V / u.m
        self.B = np.zeros(self.grid.shape) * u.T

        # Calculate a characteristic length scale from the grid to scale
        # field features with
        self.L = np.abs(np.max(self.grid.grid) - np.min(self.grid.grid))

        # Calculate radius arrays in spherical and cylindrical coordinates
        # for use in generating different field structures
        self.radius = np.sqrt(
            self.grid.xarr ** 2 + self.grid.yarr ** 2 + self.grid.zarr ** 2
        )
        self.pradius = np.sqrt(self.grid.xarr ** 2 + self.grid.yarr ** 2)

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
        self._norm()

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
        self._norm()


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
        E = -\nabla \phi \text{, for }  (r < L/2)
    .. math::
        E = 0 \text{, for } (r < L/2)
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

        Ex, Ey, Ez = np.gradient(
            potential, self.grid.xaxis, self.grid.yaxis, self.grid.zaxis
        )

        self.E[:, :, :, 0] = -1 * np.where(self.radius < 0.5 * self.L, Ex, 0)
        self.E[:, :, :, 1] = -1 * np.where(self.radius < 0.5 * self.L, Ey, 0)
        self.E[:, :, :, 2] = -1 * np.where(self.radius < 0.5 * self.L, Ez, 0)

    def _validate(self):
        if not self.grid.regular_grid:
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
        B = 1  \text{, for } (\rho < a)
    .. math::
        B = 0 \text{, for } (\rho > a)

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
            (1 - erf(self.grid.zarr / delta)) * np.exp(-((self.pradius / a) ** 2)) * u.V
        )

        Ex, Ey, Ez = np.gradient(
            potential, self.grid.xaxis, self.grid.yaxis, self.grid.zaxis
        )
        self.E[..., 0] = -Ex
        self.E[..., 1] = -Ey
        self.E[..., 2] = -Ez

    def _validate(self):
        if not self.grid.regular_grid:
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
    This function generates example fields based on analytical models as
    input for other functions.

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
        grid = PosGrid(num=num, length=length, regular_grid=regular_grid)
    elif not isinstance(grid, PosGrid):
        grid = PosGrid(grid)

    # Load the model class for the test example chosen
    models = {
        "no fields": NoFields,
        "electrostatic gaussian sphere": ElectrostaticGaussianSphere,
        "electrostatic planar shock": ElectrostaticPlanarShock,
        "axial magnetic field": AxiallyMagnetizedCylinder,
    }

    # Generate the fields by instantiating the test field object
    fields = models[model](grid, emax=emax, bmax=bmax)

    return fields.grid, fields.E, fields.B
