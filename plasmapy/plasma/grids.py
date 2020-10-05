"""
Defines the AbstractGrid class and child classes
"""

__all__ = [
    "AbstractGrid",
    "CartesianGrid",
    "IrregularCartesianGrid",
    "CylindricalGrid",
    "SphericalGrid",
]

import astropy.units as u
import numpy as np

from abc import ABC
from typing import Union


class AbstractGrid(ABC):
    """
    Abstract grid represents a 3D grid of positions. The grid is stored as an
    np.ndarray, while the units associated with each dimension are stored
    separately.
    """

    # TODO: add appropriate typing here on start, stop, num
    # start/stop can be -> number(int,float), ndarray, u.quantity, list of numbers
    # or u.quantities units can be astropy.units.core.Unit or list of same
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

        # Setup method contains other initialization stuff
        self._setup()

    def _setup(self):
        # Properties
        self._coordinate_system = None
        self._regular_grid = None

    @property
    def grid(self):
        """Grid of positional values."""
        return self._grid

    @property
    def regular_grid(self):
        """
        Value of regular_grid
        If None, calculate
        """
        if self._regular_grid is None:
            self._detect_regular_grid()
        return self._regular_grid

    @property
    def shape(self):
        return self._grid.shape

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
        Only valid if grid is regular: otherwise an exception is raised
        """
        if self.regular_grid:
            return np.mean(np.gradient(self.ax0)) * self.unit0
        else:
            raise ValueError(
                "The grid step size properties are only valid on "
                "regularly spaced grids."
            )

    @property
    def dax1(self):
        """
        Grid step size along axis 2
        Only valid if grid is regular: otherwise an exception is raised
        """
        if self.regular_grid:
            return np.mean(np.gradient(self.ax1)) * self.unit1
        else:
            raise ValueError(
                "The grid step size properties are only valid on "
                "regularly spaced grids."
            )

    @property
    def dax2(self):
        """
        Grid step size along axis 3
        Only valid if grid is regular: otherwise an exception is raised
        """
        if self.regular_grid:
            return np.mean(np.gradient(self.ax2)) * self.unit2
        else:
            raise ValueError(
                "The grid step size properties are only valid on "
                "regularly spaced grids."
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

        # If single values are given, expand to a list of appropriate length
        if isinstance(units, u.core.Unit):
            units = [units] * 3

        # Determine units from grid or units keyword (in that order)
        if hasattr(grid, "unit"):
            self._units = [grid.unit] * 3
        elif units is not None:
            self._units = units
        else:
            raise ValueError(
                "Input grid must be either a u.Quantity array "
                "or an np.ndarray with units provided in the "
                " units keyword."
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
        units=[u.dimensionless_unscaled] * 3,
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

        # Construct the axis arrays
        ax0 = np.linspace(start[0], stop[0], num=num[0], **kwargs)
        ax1 = np.linspace(start[1], stop[1], num=num[1], **kwargs)
        ax2 = np.linspace(start[2], stop[2], num=num[2], **kwargs)

        # Construct the coordinate arrays and grid
        arr0, arr1, arr2 = np.meshgrid(ax0, ax1, ax2, indexing="ij")
        grid = np.zeros([ax0.size, ax1.size, ax2.size, 3])
        grid[..., 0] = arr0
        grid[..., 1] = arr1
        grid[..., 2] = arr2

        self._units = units
        self._grid = grid

        # Check to make sure that the object created satisfies any
        # requirements: eg. units correspond to the coordinate system
        self._validate()

    def _detect_regular_grid(self, tol=1e-6):
        """
        Determine whether a grid is regular (uniformly spaced) by computing the
        variance of the grid gradients.
        """
        variance = np.zeros([3])
        dx = np.gradient(self.arr0, axis=0)
        variance[0] = np.std(dx) / np.mean(dx)
        dy = np.gradient(self.arr1, axis=1)
        variance[1] = np.std(dy) / np.mean(dy)
        dz = np.gradient(self.arr2, axis=2)
        variance[2] = np.std(dz) / np.mean(dz)

        self._regular_grid = np.allclose(variance, 0.0, atol=tol)


class CartesianGrid(AbstractGrid):
    """
    A regularly spaced Cartesian grid.
    """

    def _validate(self):
        # Check that all units are lengths
        for i in range(3):
            try:
                (self.units[i].to(u.m))
            except u.UnitConversionError:
                raise ValueError(
                    "Units of grid are not valid for a Cartesian "
                    f"grid: {self.units}."
                )

    @property
    def grid(self):
        """Grid of positional values. CartesianGrid has the unique property
        that the grid returned is dimensional, since all axes have the same
        units."""
        return self._grid * self.unit

    @property
    def xarr(self):
        """
        The 3D array of x values
        """
        return self.arr0

    @property
    def yarr(self):
        """
        The 3D array of y values
        """
        return self.arr1

    @property
    def zarr(self):
        """
        The 3D array of z values
        """
        return self.arr2

    @property
    def xaxis(self):
        """
        The x-axis
        """
        return self.ax0

    @property
    def yaxis(self):
        """
        The y-axis (only valid for a uniform grid)
        """
        return self.ax1

    @property
    def zaxis(self):
        """
        The z-axis (only valid for a uniform grid)
        """
        return self.ax2

    @property
    def dx(self):
        """
        Calculated dx (only valid for a uniform grid)
        """
        return self.dax0

    @property
    def dy(self):
        """
        Calculated dy (only valid for a uniform grid)
        """
        return self.dax1

    @property
    def dz(self):
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
        return np.sqrt(self.xaxis ** 2 + self.yaxis ** 2 + self.zaxis ** 2)


class IrregularCartesianGrid(CartesianGrid):
    """
    A Cartesian grid in which the _make_grid method produces an irregularly
    spaced grid (rather than a uniform one)
    """

    pass


class CylindricalGrid(AbstractGrid):
    """
    A grid with dimensions (r, theta, z) and units (length, angle, length)
    """

    pass


class SphericalGrid(AbstractGrid):
    """
    A grid with dimensions (r, theta, phi) and units (length, angle, angle)
    """

    pass
