"""
Defines the AbstractGrid class and child classes
"""

__all__ = [
    "AbstractGrid",
    "CartesianGrid",
    "NonUniformCartesianGrid",
    "CylindricalGrid",
    "SphericalGrid",
]

import astropy.units as u
import numpy as np
import scipy.interpolate as interp

from abc import ABC
from typing import Union


class AbstractGrid(ABC):
    """
    Abstract grid represents a 3D grid of positions. The grid is stored as an
    np.ndarray, while the units associated with each dimension are stored
    separately.
    """

    def __init__(self, *seeds, num=100, units=None, **kwargs):

        if len(seeds) == 1:
            # If one input is given, assume it is a grid
            self._load_grid(seeds[0], units=units)

        elif len(seeds) == 2:
            self._make_grid(seeds[0], seeds[1], num=num, units=units, **kwargs)

        else:
            raise TypeError(
                f"{self.__class__.__name__} takes 1 or 2 "
                f"positional arguments but {len(seeds)} were given"
            )

        # Initialize some variables
        self._is_uniform_grid = None
        self._interpolator = None

    @property
    def grid(self):
        """Grid of positional values."""
        return self._grid

    @property
    def is_uniform_grid(self):
        """
        Value of is_uniform_grid
        If None, calculate
        """
        if self._is_uniform_grid is None:
            self._detect_is_uniform_grid()
        return self._is_uniform_grid

    @property
    def interpolator(self):
        """
        A nearest-neighbor interpolator that returns the nearest grid index
        to a position
        """
        if self._interpolator is None:
            self._make_interpolator()
        return self._interpolator

    def _make_interpolator(self):
        """
        Initializes a nearest-neighbor interpolator that returns the nearest
        grid indices for a given position.
        """
        # Create a grid of indices for use in interpolation
        n0, n1, n2, ndim = self.shape
        indgrid = np.indices([n0, n1, n2])
        indgrid = np.moveaxis(indgrid, 0, -1)

        # Create an input array of grid positions
        pts = (
            self.ax0.si.value,
            self.ax1.si.value,
            self.ax2.si.value,
        )

        # Create the interpolator
        # Note that this regular grid interpolator requires points be a regular
        # grid but does not require uniform spacing

        self._interpolator = interp.RegularGridInterpolator(
            pts, indgrid, method="nearest", bounds_error=False, fill_value=np.nan
        )

    def interpolate_indices(self, pos: Union[np.ndarray, u.Quantity], *args):
        r"""
        Interpolate the nearest grid indices to a position using a
        nearest-neighbor interpolator

        Parameters
        ----------
        pos : np.ndarray or u.Quantity array, shape (n,3)
            An array of positions in space, where the second dimension
            corresponds to the three dimensions of the grid. If an np.ndarray
            is provided, units will be assumed to match those of the grid.

        """
        # Condition pos
        # If a single point was given, add empty dimension
        if pos.ndim == 1:
            pos = np.reshape(pos, [1, 3])
        pos2 = np.zeros(pos.shape)
        # Convert position to SI and then strip units
        if hasattr(pos, "unit"):
            pos2 = pos.si.value
        else:
            for i in range(3):
                pos2[:, i] = (pos[:, i] * self.units[i]).si.value

        # Interpolate indices
        i = self.interpolator(pos2)
        # Convert any non-NaN values to ints
        i = i.astype(np.int32)

        return i

    def nearest_neighbor_interpolator(self, pos: Union[np.ndarray, u.Quantity], *args):
        r"""
        Interpolate values on the grid using a nearest-neighbor scheme with
        no higher-order weighting.

        Parameters
        ----------
        pos : np.ndarray or u.Quantity array, shape (n,3)
            An array of positions in space, where the second dimension
            corresponds to the three dimensions of the grid. If an np.ndarray
            is provided, units will be assumed to match those of the grid.

        args : np.ndarray or u.Quantity objects of shape (n0,n1,2) or (n0,n1,n2,3)
            Arrays of values at each grid vertex to be interpolated onto the
            positions given.

        """
        # pos is validated in interpolate_indices

        # Validate args
        # must be np.ndarray or u.Quantity arrays of same shape as grid
        for arg in args:
            if (
                isinstance(arg, (u.Quantity, np.ndarray))
                and arg.shape[0:-1] == self.grid_shape
            ):
                pass
            else:
                raise ValueError(
                    "Argument must be a u.Quantity or np.ndarray "
                    "with the same shape as the grid, but argument "
                    f"given had shape {arg.shape} and type "
                    f"{type(arg)}"
                )

        i = self.interpolate_indices(pos)

        output = []
        for arg in args:
            values = arg[i[:, 0], i[:, 1], i[:, 2], ...]
            output.append(values)

        return tuple(output)

    def volume_averaged_interpolator(self, pos: Union[np.ndarray, u.Quantity], *args):
        r"""
        Interpolate values on the grid using a volume-averaged scheme with
        no higher-order weighting.

        Parameters
        ----------
        pos : np.ndarray or u.Quantity array, shape (n,3)
            An array of positions in space, where the second dimension
            corresponds to the three dimensions of the grid. If an np.ndarray
            is provided, units will be assumed to match those of the grid.

        args : np.ndarray or u.Quantity objects of shape (n0,n1,2) or (n0,n1,n2,3)
            Arrays of values at each grid vertex to be interpolated onto the
            positions given.

        """

        raise NotImplemented(
            "Volume-averaged interpolator is not yet " "implemented for this grid type."
        )

    @property
    def shape(self):
        return self._grid.shape

    @property
    def grid_shape(self):
        """
        The shape of the grid without the number of dimensions,
        eg. (n0, n1, n2)
        """
        return self._grid.shape[0:-1]

    @property
    def unit(self):
        """
        The unit for the entire grid. Only valid if all dimensions of the
        grid have the same units: otherwise, an exception is raised.
        """
        if self._units[0] == self._units[1] and self._units[0] == self._units[2]:
            return self._units[0]
        else:
            raise ValueError(
                "Array dimensions do not all have the same " f"units: {self._units}"
            )

    @property
    def units(self):
        """The list of units corresponding to the dimensions"""
        return self._units

    @property
    def unit0(self):
        """Unit of dimension 1"""
        return self._units[0]

    @property
    def unit1(self):
        """Unit of dimension 2"""
        return self._units[1]

    @property
    def unit2(self):
        """Unit of dimension 3"""
        return self._units[2]

    @property
    def ax0(self):
        """Axis 1"""
        return self._grid[:, 0, 0, 0] * self.unit0

    @property
    def ax1(self):
        """Axis 2"""
        return self._grid[0, :, 0, 1] * self.unit1

    @property
    def ax2(self):
        """Axis 3"""
        return self._grid[0, 0, :, 2] * self.unit2

    @property
    def dax0(self):
        """
        Grid step size along axis 1
        Only valid if grid is uniform: otherwise an exception is raised
        """
        if self.is_uniform_grid:
            return np.mean(np.gradient(self.ax0))
        else:
            raise ValueError(
                "The grid step size properties are only valid on "
                "uniformly spaced grids."
            )

    @property
    def dax1(self):
        """
        Grid step size along axis 2
        Only valid if grid is uniform: otherwise an exception is raised
        """
        if self.is_uniform_grid:
            return np.mean(np.gradient(self.ax1))
        else:
            raise ValueError(
                "The grid step size properties are only valid on "
                "uniformly spaced grids."
            )

    @property
    def dax2(self):
        """
        Grid step size along axis 3
        Only valid if grid is uniform: otherwise an exception is raised
        """
        if self.is_uniform_grid:
            return np.mean(np.gradient(self.ax2))
        else:
            raise ValueError(
                "The grid step size properties are only valid on "
                "uniformly spaced grids."
            )

    @property
    def arr0(self):
        """Array of positions in dimension 1"""
        return self._grid[..., 0] * self.unit0

    @property
    def arr1(self):
        """Array of positions in dimension 2"""
        return self._grid[..., 1] * self.unit1

    @property
    def arr2(self):
        """Array of positions in dimension 3"""
        return self._grid[..., 2] * self.unit2

    def _validate(self):
        """
        Checks to make sure that the grid parameters are
        consistent with the coordinate system and units selected

        """
        return True

    def _load_grid(self, grid: Union[np.ndarray, u.Quantity], units=None):
        """
        Initialize the grid object from a user-supplied grid

        Parameters
        ----------
        grid : np.ndarray or u.Quantity array, shape (n0, n1, n2, 3)
            A grid of positions in 3D space. If an np.ndarray is provided,
            the units keyword must be set to provide units for each dimension.
            If a u.Quantity array is provided, its units will be used.

        units : A u.core.Unit object or a list of three of the same
            Units to be applied to each dimension. If only one unit is given,
            the same unit will be applied to all dimensions.

        Returns
        -------
        None.

        """

        # Determine units from grid or units keyword (in that order)
        if hasattr(grid, "unit"):
            self._units = [grid.unit] * 3
        elif units is not None:
            self._units = units
        else:
            raise ValueError(
                "Input grid must be either a u.Quantity array "
                "or an np.ndarray with units provided in the "
                " units keyword. Type of grid provided was "
                f"{type(grid)}."
            )

        # If quantity array was given, strip units
        if hasattr(grid, "value"):
            self._grid = grid.value
        else:
            self._grid = grid

        # Check to make sure that the object created satisfies any
        # requirements: eg. units correspond to the coordinate system
        self._validate()

    def _make_grid(
        self,
        start: Union[int, float, u.Quantity],
        stop: Union[int, float, u.Quantity],
        num=100,
        units=None,
        **kwargs,
    ):
        """
        Creates a grid based on start, stop, and num values in a manner
        that mirrors the interface of the np.linspace function.

        Parameters
        ----------
        start : number (int,float,u.Quantity) or a list of three of the same
            Starting values for each dimension. If one value is given,
            the same value will be used for all three dimensions.

        stop : number (int,float,u.Quantity) or a list of three of the same
            End values for each dimension. If one value is given,
            the same value will be used for all three dimensions.

        num : int or list of three ints, optional
            The number of points in each dimension. If a single integer is
            given, the same number of points will be used in each dimension.
            The default is 100.

        units: u.core.Unit object or array of three of the same
            Units for each array dimension. Overridden by units on stop
            and start, if present.

        **kwargs: Additional arguments
            Any additional arguments will be passed directly to np.linspace()

        Returns
        -------
        None.

        """
        if units is None:
            units = [u.dimensionless_unscaled] * 3

        # If single values are given, expand to a list of appropriate length
        if isinstance(stop, (int, float, u.Quantity)):
            stop = [stop] * 3
        if isinstance(start, (int, float, u.Quantity)):
            start = [start] * 3
        if isinstance(num, (int, float, u.Quantity)):
            num = [num] * 3
        if isinstance(units, u.core.Unit):
            units = [units] * 3

        # Check to make sure all lists now contain three values
        # (throws exception if user supplies a list of two, say)
        var = {"stop": stop, "start": start, "num": num, "units": units}
        for k in var.keys():
            if len(var[k]) != 3:
                raise ValueError(
                    f"{k} must be either a single value or a "
                    "list of three values, but "
                    f"({len(var[k])} values were given)."
                )

        # Extract units from input arrays (if they are there), then
        # remove the units from those arrays
        for i in range(3):
            # Determine unit for dimension
            if hasattr(stop[i], "unit"):
                units[i] = stop[i].unit
            elif hasattr(start[i], "unit"):
                units[i] = start[i].unit
            else:
                # Use units provided by units keyword
                pass

            # Attempt to convert start and stop to unit chosen
            try:
                start[i] = start[i].to(units[i])
            except u.UnitConversionError:
                raise ValueError(
                    f"Units of {start[i]} and " f" {units[i]} are not compatible"
                )
            except AttributeError:
                # Assume exception was raised because value is not a u.Quantity
                start[i] *= units[i]

            try:
                stop[i] = stop[i].to(units[i])
            # Note that UnitConversionError is impossible for stop
            except AttributeError:
                # Assume exception was raised because value is not a u.Quantity
                stop[i] *= units[i]

            # strip units
            stop[i] = stop[i].value
            start[i] = start[i].value

        # Create coordinate mesh
        arr0, arr1, arr2 = self._make_mesh(start, stop, num, **kwargs)
        n0, n1, n2 = arr0.shape

        grid = np.zeros([n0, n1, n2, 3])
        grid[..., 0] = arr0
        grid[..., 1] = arr1
        grid[..., 2] = arr2

        self._units = units
        self._grid = grid

        # Check to make sure that the object created satisfies any
        # requirements: eg. units correspond to the coordinate system
        self._validate()

    def _make_mesh(self, start, stop, num, **kwargs):
        """
        Creates mesh as part of _make_grid(). Separated into its own function
        so it can be re-implemented to make non-uniformly spaced meshes
        """
        # Construct the axis arrays
        ax0 = np.linspace(start[0], stop[0], num=num[0], **kwargs)
        ax1 = np.linspace(start[1], stop[1], num=num[1], **kwargs)
        ax2 = np.linspace(start[2], stop[2], num=num[2], **kwargs)

        # Construct the coordinate arrays
        arr0, arr1, arr2 = np.meshgrid(ax0, ax1, ax2, indexing="ij")

        return arr0, arr1, arr2

    def _detect_is_uniform_grid(self, tol=1e-6):
        """
        Determine whether a grid is uniform (uniformly spaced) by computing the
        variance of the grid gradients.
        """
        variance = np.zeros([3])
        dx = np.gradient(self.arr0, axis=0)
        variance[0] = np.std(dx) / np.mean(dx)
        dy = np.gradient(self.arr1, axis=1)
        variance[1] = np.std(dy) / np.mean(dy)
        dz = np.gradient(self.arr2, axis=2)
        variance[2] = np.std(dz) / np.mean(dz)

        self._is_uniform_grid = np.allclose(variance, 0.0, atol=tol)


class CartesianGrid(AbstractGrid):
    """
    A uniformly spaced Cartesian grid.
    """

    def _validate(self):
        # Check that all units are lengths
        for i in range(3):
            try:
                self.units[i].to(u.m)
            except u.UnitConversionError:
                raise ValueError(
                    "Units of grid are not valid for a Cartesian "
                    f"grid: {self.units}."
                )

    @property
    def grid(self):
        """Grid of positional values. CartesianGrid's have the unique
        property that the grid returned has units, since all dimensions are
        of the same unit."""
        return self._grid * self.unit

    @property
    def x_arr(self):
        """
        The 3D array of x values
        """
        return self.arr0

    @property
    def y_arr(self):
        """
        The 3D array of y values
        """
        return self.arr1

    @property
    def z_arr(self):
        """
        The 3D array of z values
        """
        return self.arr2

    @property
    def x_axis(self):
        """
        The x-axis
        """
        return self.ax0

    @property
    def y_axis(self):
        """
        The y-axis (only valid for a uniform grid)
        """
        return self.ax1

    @property
    def z_axis(self):
        """
        The z-axis (only valid for a uniform grid)
        """
        return self.ax2

    @property
    def d_x(self):
        """
        Calculated dx (only valid for a uniform grid)
        """
        return self.dax0

    @property
    def d_y(self):
        """
        Calculated dy (only valid for a uniform grid)
        """
        return self.dax1

    @property
    def d_z(self):
        """
        Calculated dz (only valid for a uniform grid)
        """
        return self.dax2

    @property
    def distance_from_origin(self):
        """
        The 3D array of the radial position of each point from the
        origin
        """
        return np.sqrt(self.x_axis ** 2 + self.y_axis ** 2 + self.z_axis ** 2)

    def volume_averaged_interpolator(self, pos: Union[np.ndarray, u.Quantity], *args):
        r"""
        Interpolate values on the grid using a nearest-neighbor scheme with
        no higher-order weighting.

        Parameters
        ----------
        pos :  u.Quantity array, shape (n,3)
            An array of positions in space, where the second dimension
            corresponds to the three dimensions of the grid. If an np.ndarray
            is provided, units will be assumed to match those of the grid.

        args : np.ndarray or u.Quantity objects of shape (n0,n1,2) or (n0,n1,n2,3)
            Arrays of values at each grid vertex to be interpolated onto the
            positions given.

        """

        # Condition pos
        # If a single point was given, add empty dimension
        if pos.ndim == 1:
            pos = np.reshape(pos, [1, 3])

        # Convert position to u.Quantiy
        if not hasattr(pos, "unit"):
            pos *= self.unit

        # Interpolate the indices
        i = self.interpolate_indices(pos)
        nparticles = i.shape[0]

        # Calculate the grid positions for each particle as interpolated
        # by the nearest neighbor interpolator
        gpos = self.grid[i[:, 0], i[:, 1], i[:, 2], :]

        # Determine the points bounding the grid cell containing the
        # particle
        x0 = np.where(pos[:, 0] > gpos[:, 0], i[:, 0], i[:, 0] - 1)
        x1 = x0 + 1
        y0 = np.where(pos[:, 1] > gpos[:, 1], i[:, 1], i[:, 1] - 1)
        y1 = y0 + 1
        z0 = np.where(pos[:, 2] > gpos[:, 2], i[:, 2], i[:, 2] - 1)
        z1 = z0 + 1

        # Calculate the cell volume
        cell_vol = self.d_x * self.d_y * self.d_z
        n0, n1, n2 = self.grid_shape

        # Create a list of empty arrays to hold final results
        output = []
        for i, arg in enumerate(args):
            ndim = arg.shape[-1]
            output.append(np.zeros([nparticles, ndim]) * arg.unit)

        # Go through all of the vertices around the position and volume-
        # weight the values
        for x in [x0, x1]:
            for y in [y0, y1]:
                for z in [z0, z1]:

                    # Determine if gridpoint is within bounds
                    valid = (
                        (x >= 0) & (x < n0) & (y >= 0) & (y < n1) & (z >= 0) & (z < n2)
                    )
                    out = np.where(valid == False)

                    # Distance from grid vertex to particle position
                    d = np.abs(self.grid[x, y, z, :] - pos)

                    # Fraction of cell volume that is closest to the
                    # current point
                    weight = (d[:, 0] * d[:, 1] * d[:, 2]) / cell_vol
                    weight = weight.to(u.dimensionless_unscaled)
                    weight[out] = 0
                    weight = np.outer(weight, np.ones(3))

                    # For each argument, include the contributed by this
                    # grid vertex
                    for i, arg in enumerate(args):
                        output[i] += weight * arg[x, y, z, ...]

        return tuple(output)


class NonUniformCartesianGrid(CartesianGrid):
    """
    A Cartesian grid in which the _make_mesh method produces a non-uniformly
    spaced grid.
    """

    def _make_mesh(self, start, stop, num, **kwargs):
        """
        Creates mesh as part of _make_grid(). Separated into its own function
        so it can be re-implemented to make non-uniform grids.
        """
        # Construct the axis arrays
        ax0 = np.sort(np.random.uniform(low=start[0], high=stop[0], size=num[0]))

        ax1 = np.sort(np.random.uniform(low=start[1], high=stop[1], size=num[1]))

        ax2 = np.sort(np.random.uniform(low=start[2], high=stop[2], size=num[2]))

        # Construct the coordinate arrays
        arr0, arr1, arr2 = np.meshgrid(ax0, ax1, ax2, indexing="ij")

        return arr0, arr1, arr2


class CylindricalGrid(AbstractGrid):
    """
    A grid with dimensions (rho, theta, z) and units (length, angle, length)
    """

    def _validate(self):
        # Check that units are consistent with cylindrical coordinates
        units = [u.m, u.rad, u.m]
        for i, unit in enumerate(units):
            try:
                self.units[i].to(unit)
            except u.UnitConversionError:
                raise ValueError(
                    "Units of grid are not valid for a Cylindrical grid: "
                    "dimensions must be (length, angle, length) but given "
                    f"dimensions are {self.units}."
                )

        # Check that angular dimensions have no values the maximum angle for
        # that dimension
        if np.max(self.theta_arr.to(u.rad).value) > np.pi * 2:
            raise ValueError(
                "Theta dimension of cylindrical grid can contain no values "
                "greater than 2pi radians."
            )

    @property
    def rho_arr(self):
        """
        The 3D array of rho values
        """
        return self.arr0

    @property
    def theta_arr(self):
        """
        The 3D array of theta values
        """
        return self.arr1

    @property
    def z_arr(self):
        """
        The 3D array of z values
        """
        return self.arr2

    @property
    def rho_axis(self):
        """
        The rho-axis  (only valid for a uniform grid)
        """
        return self.ax0

    @property
    def theta_axis(self):
        """
        The theta-axis  (only valid for a uniform grid)
        """
        return self.ax1

    @property
    def z_axis(self):
        """
        The z-axis  (only valid for a uniform grid)
        """
        return self.ax2

    @property
    def d_rho(self):
        """
        Calculated drho  (only valid for a uniform grid)
        """
        return self.dax0

    @property
    def d_theta(self):
        """
        Calculated dtheta (only valid for a uniform grid)
        """
        return self.dax1

    @property
    def d_z(self):
        """
        Calculated dz (only valid for a uniform grid)
        """
        return self.dax2


class SphericalGrid(AbstractGrid):
    """
    A grid with dimensions (r, theta, phi) and units (length, angle, angle)
    """

    def _validate(self):
        # Check that units are consistent with cylindrical coordinates
        units = [u.m, u.rad, u.rad]
        for i, unit in enumerate(units):
            try:
                self.units[i].to(unit)
            except u.UnitConversionError:
                raise ValueError(
                    "Units of grid are not valid for a Cylindrical grid: "
                    "dimensions must be (length, angle, length) but given "
                    f"dimensions are {self.units}."
                )

        # Check that angular dimensions have no values the maximum angle for
        # that dimension
        if np.max(self.theta_arr.to(u.rad).value) > np.pi:
            raise ValueError(
                "Theta dimension of spherical grid can contain no values "
                "greater than pi radians."
            )

        if np.max(self.phi_arr.to(u.rad).value) > np.pi * 2:
            raise ValueError(
                "Phi dimension of spherical grid can contain no values "
                "greater than 2pi radians."
            )

    @property
    def r_arr(self):
        """
        The 3D array of r values
        """
        return self.arr0

    @property
    def theta_arr(self):
        """
        The 3D array of theta values
        """
        return self.arr1

    @property
    def phi_arr(self):
        """
        The 3D array of phi values
        """
        return self.arr2

    @property
    def r_axis(self):
        """
        The r-axis  (only valid for a uniform grid)
        """
        return self.ax0

    @property
    def theta_axis(self):
        """
        The theta-axis  (only valid for a uniform grid)
        """
        return self.ax1

    @property
    def phi_axis(self):
        """
        The phi-axis  (only valid for a uniform grid)
        """
        return self.ax2

    @property
    def d_r(self):
        """
        Calculated dr  (only valid for a uniform grid)
        """
        return self.dax0

    @property
    def d_theta(self):
        """
        Calculated dtheta (only valid for a uniform grid)
        """
        return self.dax1

    @property
    def d_phi(self):
        """
        Calculated dphi (only valid for a uniform grid)
        """
        return self.dax2
