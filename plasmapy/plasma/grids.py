"""
Defines the AbstractGrid class and child classes
"""

__all__ = [
    "AbstractGrid",
    "CartesianGrid",
    "NonUniformCartesianGrid",
]

import astropy.units as u
import numpy as np
import pandas as pd
import scipy.interpolate as interp

import warnings
import xarray as xr

from abc import ABC
from typing import Union

# These standard keys are used throughout PlasmaPy to refer to certain
# plasma quantities. This dictionary also provides the expected unit.
recognized_keys = {
    "x": ("x spatial position", u.m),
    "y": ("y spatial position", u.m),
    "z": ("z spatial position", u.m),
    "rho": ("Mass density", u.kg / u.m ** 3),
    "E_x": ("Electric field (x component)", u.V / u.m),
    "E_y": ("Electric field (y component)", u.V / u.m),
    "E_z": ("Electric field (z component)", u.V / u.m),
    "B_x": ("Magnetic field (x component)", u.T),
    "B_y": ("Magnetic field (y component)", u.T),
    "B_z": ("Magnetic field (z component)", u.T),
}


def _detect_is_uniform_grid(pts0, pts1, pts2, tol=1e-6):
    r"""
    Determine whether a grid is uniform (uniformly spaced) by computing the
    variance of the grid gradients.
    """
    variance = np.zeros([3])
    dx = np.gradient(pts0, axis=0)
    variance[0] = np.std(dx) / np.mean(dx)
    dy = np.gradient(pts1, axis=1)
    variance[1] = np.std(dy) / np.mean(dy)
    dz = np.gradient(pts2, axis=2)
    variance[2] = np.std(dz) / np.mean(dz)

    return np.allclose(variance, 0.0, atol=tol)


class AbstractGrid(ABC):
    r"""
    Abstract grid represents a 3D grid of positions. The grid is stored as an
    np.ndarray, while the units associated with each dimension are stored
    separately.
    """

    def __init__(self, *seeds, num=100, **kwargs):


        # Initialize some variables
        self._is_uniform_grid = None
        self._interpolator = None
        self._grids = None  # [nx,ny,nz] x 3
        self._grid = None

        # If three inputs are given, assume it's a user-provided grid
        if len(seeds) == 3:
            self._load_grid(seeds[0], seeds[1], seeds[2])

        # If two inputs are given, assume they are start and stop arrays
        # to create a new grid
        elif len(seeds) == 2:
            self._make_grid(seeds[0], seeds[1], num=num, **kwargs)


        else:
            raise TypeError(
                f"{self.__class__.__name__} takes 2 or 3 "
                f"positional arguments but {len(seeds)} were given"
            )

    def _validate(self):
        r"""
        Checks to make sure that the grid parameters are
        consistent with the coordinate system and units selected

        """
        return True

    # *************************************************************************
    # Fundamental properties of the grid
    # *************************************************************************

    def __repr__(self):

        line_sep = "-----------------------------\n"

        s = f"*** Grid summary ***\n{type(self)}\n"
        s += f"Shape: {self.shape}\nUnits: {self.units}\n"

        if self.is_uniform_grid:
            s += (
                "Uniformly Spaced, dx,dy,dz = "
                f"({self.dax0:.3f},{self.dax1:.3f},{self.dax2:.3f})\n"
            )
        else:
            s += "Non-Uniform Spacing\n"

        keys = list(self.ds.data_vars)
        rkeys = [k for k in keys if k in recognized_keys]
        nrkeys = [k for k in keys if k not in recognized_keys]

        s += line_sep + "Recognized Quantities:\n"
        if len(rkeys) == 0:
            s += "-None-\n"
        else:
            for key in rkeys:
                unit = self.ds[key].attrs["unit"]
                s += f"-> {key} ({unit})\n"

        s += line_sep + "Unrecognized Quantities:\n"
        if len(nrkeys) == 0:
            s += "-None-\n"
        else:
            for key in nrkeys:
                unit = self.ds[key].attrs["unit"]
                s += f"-> {key} ({unit})\n"

        return s

    @property
    def shape(self):
        r""" Shape of the grid"""
        if self.is_uniform_grid:
            return (self.ax0.size, self.ax1.size, self.ax2.size)
        else:
            return self.ds.coords["ax0"].shape

    @property
    def grids(self):
        r"""Grids of vertex positions"""
        if self._grids is None:

            if self.is_uniform_grid:
                pts0, pts1, pts2 = np.meshgrid(
                    self.ax0, self.ax1, self.ax2, indexing="ij"
                )
                self._grids = (pts0, pts1, pts2)
            else:
                self._grids = (
                    self.ds["ax0"].data,
                    self.ds["ax1"].data,
                    self.ds["ax2"].data,
                )

        return self._grids

    @property
    def grid(self):
        r"""A single grid of vertex positions"""

        if self._grid is None:
            pts0, pts1, pts2 = self.grids
            if self.is_uniform_grid:
                n0, n1, n2 = pts0.shape
                grid = np.zeros([n0, n1, n2, 3])
            else:
                n = pts0.size
                grid = np.zeros([n, 3])

            grid[..., 0] = pts0
            grid[..., 1] = pts1
            grid[..., 2] = pts2
            self._grid = grid

        return self._grid

    @property
    def pts0(self):
        r"""Array of positions in dimension 1"""
        return self.grids[0] * self.unit0

    @property
    def pts1(self):
        r"""Array of positions in dimension 2"""
        return self.grids[1] * self.unit1

    @property
    def pts2(self):
        r"""Array of positions in dimension 3"""
        return self.grids[2] * self.unit2

    @property
    def units(self):
        r"""Returns a list of the units of each dimension"""
        return self.ds.attrs["axis_units"]

    @property
    def unit0(self):
        r"""Unit of dimension 1"""
        return self.units[0]

    @property
    def unit1(self):
        r"""Unit of dimension 2"""
        return self.units[1]

    @property
    def unit2(self):
        r"""Unit of dimension 3"""
        return self.units[2]

    @property
    def unit(self):
        r"""
        The unit for the entire grid. Only valid if all dimensions of the
        grid have the same units: otherwise, an exception is raised.
        """
        if self.units[0] == self.units[1] and self.units[0] == self.units[2]:
            return self.units[0]
        else:
            raise ValueError(
                "Array dimensions do not all have the same " f"units: {self.units}"
            )

    # *************************************************************************
    # 1D axes and step sizes (valid only for uniform grids)
    # *************************************************************************

    @property
    def ax0(self):
        r"""
        Axis 1: Only valid if grid is uniform: otherwise an exception is raised
        """

        if self.is_uniform_grid:
            return self.ds.coords["ax0"].values * self.unit0
        else:
            raise ValueError(
                "The axis properties are only valid on " "uniformly spaced grids."
            )

    @property
    def ax1(self):
        r"""Axis 2: Only valid if grid is uniform: otherwise an exception is raised"""
        if self.is_uniform_grid:
            return self.ds.coords["ax1"].values * self.unit1
        else:
            raise ValueError(
                "The axis properties are only valid on " "uniformly spaced grids."
            )

    @property
    def ax2(self):
        r"""
        Axis 3
        Only valid if grid is uniform: otherwise an exception is raised
        """
        if self.is_uniform_grid:
            return self.ds.coords["ax2"].values * self.unit2
        else:
            raise ValueError(
                "The axis properties are only valid on " "uniformly spaced grids."
            )

    @property
    def dax0(self):
        r"""
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
        r"""
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
        r"""
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

    # *************************************************************************
    # Loading and creating grids
    # *************************************************************************

    def _load_grid(
        self, pts0: u.Quantity, pts1: u.Quantity, pts2: u.Quantity, **kwargs
    ):
        r"""
        Initialize the grid object from a user-supplied grid

        Parameters
        ----------
        grid{0,1,2} : u.Quantity array, shape (n0, n1, n2)
            Grids of coordinate positions.

        **kwargs: u.Quantity array, shape (n0, n1, n2)
            Quantities defined on the grid

        Returns
        -------
        None.

        """

        # Validate input
        if not (pts0.shape == pts1.shape and pts0.shape == pts2.shape):
            raise ValueError(
                "Provided arrays of grid points are of unequal "
                f"shape: pts0 = {pts0.shape}, "
                f"pts1 = {pts1.shape}, "
                f"pts2 = {pts2.shape}."
            )

        self.is_uniform_grid = _detect_is_uniform_grid(pts0, pts1, pts2)

        # Create dataset
        self.ds = xr.Dataset()

        self.ds.attrs["axis_units"] = [pts0.unit, pts1.unit, pts2.unit]
        if self.is_uniform_grid:
            self.ds.coords["ax0"] = pts0[:, 0, 0]
            self.ds.coords["ax1"] = pts1[0, :, 0]
            self.ds.coords["ax2"] = pts2[0, 0, :]

        else:
            mdx = pd.MultiIndex.from_arrays(
                [pts0.flatten(), pts1.flatten(), pts2.flatten()],
                names=["ax0", "ax1", "ax2"],
            )
            self.ds.coords["ax"] = mdx

        # Add quantities
        for qk in kwargs.keys():
            q = kwargs[qk]

            self.add_quantity(qk, q)

        # Check to make sure that the object created satisfies any
        # requirements: eg. units correspond to the coordinate system
        self._validate()

    def add_quantity(self, key: str, quantity: u.Quantity):
        r"""
        Adds a quantity to the dataset as a new DataArray

        Parameters
        ----------
        key : str
            Key string that labels the quantity to be added.
        quantity : u.Quantity
            A u.Quantity object whose shape must match that of the grid

        Returns
        -------
        None.

        """

        # Check key against a list of "known" keys with pre-defined
        # meanings (eg. E_x, n_e) and raise a warning if a "non-standard"
        # key is being used so the user is aware.
        if key in recognized_keys.keys():
            try:
                quantity.to(recognized_keys[key][1])
            except u.UnitConversionError:
                raise ValueError(
                    f"Units provided for {key} ({quantity.unit}) "
                    "are not compatible with the correct units "
                    f"for that recognized key ({recognized_keys[key]})."
                )

        else:
            warnings.warn(
                f"Warning: {key} is not recognized quantity key", stacklevel=2
            )

        if self.is_uniform_grid:
            axes = ["ax0", "ax1", "ax2"]
        # If grid is non-uniform, flatten quantity
        else:
            quantity = quantity.flatten()
            axes = ["ax"]

        if quantity.shape != self.shape:
            raise ValueError(
                f"Shape of quantity '{key}' {quantity.shape} "
                f"does not match the grid shape {self.shape}."
            )

        data = xr.DataArray(quantity, dims=axes, attrs={"unit": quantity.unit})
        self.ds[key] = data


    def add_quantities(self, keys: list, quantities: list):
        r"""
        Adds a list of keys and quantities to the grid. See "add_quantity"

        Parameters
        ----------
        keys : list
            A list of str keys for quantities to add
        quantities : list
            A list of quantities to add. Each must be an array that matches
            the size of the grid

        Returns
        -------
        None.

        """
        if len(keys) != len(quantities):
            raise ValueError("Number of keys must equal number of quantites")

        for key, quantity in zip(keys, quantities):
            self.add_quantity(key, quantity)

    def _make_grid(
        self,
        start: Union[int, float, u.Quantity],
        stop: Union[int, float, u.Quantity],
        num=100,
        units=None,
        **kwargs,
    ):
        r"""
        Creates a grid based on start, stop, and num values in a manner
        that mirrors the interface of the np.linspace function.

        Parameters
        ----------
        start : number (u.Quantity) or a list of three of the same
            Starting values for each dimension. If one value is given,
            the same value will be used for all three dimensions.

        stop : number (u.Quantity) or a list of three of the same
            End values for each dimension. If one value is given,
            the same value will be used for all three dimensions.

        num : int or list of three ints, optional
            The number of points in each dimension. If a single integer is
            given, the same number of points will be used in each dimension.
            The default is 100.

        **kwargs: Additional arguments
            Any additional arguments will be passed directly to np.linspace()

        Returns
        -------
        None.

        """

        # If array of quantities are given instead of a list, convert
        if isinstance(stop, u.Quantity) and stop.size == 3:
            stop = list(stop)
        elif isinstance(stop, u.Quantity) and stop.size == 1:
            stop = [stop] * 3

        if isinstance(start, u.Quantity) and start.size > 1:
            start = list(start)
        elif isinstance(start, u.Quantity) and start.size == 1:
            start = [start] * 3

        if isinstance(num, (int, float)):
            num = [int(num)] * 3

        # Check to make sure all lists now contain three values
        # (throws exception if user supplies a list of two, say)
        var = {"stop": stop, "start": start, "num": num}
        for k in var.keys():
            if len(var[k]) != 3:
                raise ValueError(
                    f"{k} must be either a single value or a "
                    "list of three values, but "
                    f"({len(var[k])} values were given)."
                )

        # Extract units from input arrays (if they are there), then
        # remove the units from those arrays
        units = []
        for i in range(3):
            # Determine unit for dimension
            unit = start[i].unit
            units.append(unit)

            # Attempt to convert stop unit to start unit
            try:
                stop[i] = stop[i].to(unit)

            except u.UnitConversionError:
                raise ValueError(
                    f"Units of {stop[i]} and " f" {unit} are not compatible"
                )
            except AttributeError:
                raise AttributeError(
                    "Start and stop values must be u.Quantity instances."
                )

            # strip units
            stop[i] = stop[i].value
            start[i] = start[i].value

        # Create coordinate mesh
        pts0, pts1, pts2 = self._make_mesh(start, stop, num, **kwargs)

        # Load into the dataset using the _load_grid function
        self._load_grid(
            pts0 * units[0], pts1 * units[1], pts2 * units[2],
        )

    def _make_mesh(self, start, stop, num, **kwargs):
        r"""
        Creates mesh as part of _make_grid(). Separated into its own function
        so it can be re-implemented to make non-uniformly spaced meshes
        """
        # Construct the axis arrays
        ax0 = np.linspace(start[0], stop[0], num=num[0], **kwargs)
        ax1 = np.linspace(start[1], stop[1], num=num[1], **kwargs)
        ax2 = np.linspace(start[2], stop[2], num=num[2], **kwargs)

        # Construct the coordinate arrays
        pts0, pts1, pts2 = np.meshgrid(ax0, ax1, ax2, indexing="ij")

        return pts0, pts1, pts2

    # *************************************************************************
    # Interpolators
    # *************************************************************************

    @property
    def interpolator(self):
        r"""
        A nearest-neighbor interpolator that returns the nearest grid index
        to a position
        """
        if self._interpolator is None:
            if self.is_uniform_grid:
                self._make_uniform_grid_interpolator()
            else:
                self._make_nonuniform_grid_interpolator()

        return self._interpolator

    def _make_uniform_grid_interpolator(self):
        r"""
        Initializes a nearest-neighbor interpolator that returns the nearest
        grid indices for a given position (given in SI units).

        This function works on a uniformly spaced grid
        """
        # Create a grid of indices for use in interpolation
        n0, n1, n2 = self.shape
        indgrid = np.indices([n0, n1, n2])
        indgrid = np.moveaxis(indgrid, 0, -1)

        # Create an input array of grid positions
        pts = (
            self.ax0.si.value,
            self.ax1.si.value,
            self.ax2.si.value,
        )

        self._interpolator = interp.RegularGridInterpolator(
            pts, indgrid, method="nearest", bounds_error=False, fill_value=np.nan
        )

    def _make_nonuniform_grid_interpolator(self):
        r"""
        Initializes a nearest-neighbor interpolator that returns the nearest
        grid indices for a given position (given in SI units).

        This function works on unstructured (non-uniform) data
        """

        # Make an array of point positions
        pts0, pts1, pts2 = self.pts0.si.value, self.pts1.si.value, self.pts2.si.value
        pts = np.array([pts0, pts1, pts2])
        pts = np.moveaxis(pts, 0, 1)

        # Create a flat array of indices corresponding to those positions
        indgrid = np.arange(self.shape[0])

        self._interpolator = interp.NearestNDInterpolator(pts, indgrid)

    def interpolate_indices(self, pos: Union[np.ndarray, u.Quantity]):
        r"""
        Interpolate the nearest grid indices to a position using a
        nearest-neighbor interpolator

        Parameters
        ----------
        pos : np.ndarray or u.Quantity array, shape (n,3)
            An array of positions in space, where the second dimension
            corresponds to the three dimensions of the grid. If an np.ndarray
            is provided, units will be assumed to match those of the grid.

        Returns
        -------

        i : np.ndarray, shape (n,3)
            An array of indices corresponding to the positions such that
            i[n,:] = ix,iy,iz such that grid[ix,iy,iz,:] ~ pos[n,:]

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

        # TODO: Check interpolated positions and reject any (set to NaN)
        # that are above a certain tolerance distance?
        # currently the nonuniform interpolator can't tell when a value
        # is out of bounds...

        # Note: i contains nan values which must be replaced with 0's with
        # appropriate units in the second layer interpolator functions.

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

        *args : str
            Strings that correspond to DataArrays in the dataset

        """
        # pos is validated in interpolate_indices

        # Validate args
        # must be np.ndarray or u.Quantity arrays of same shape as grid
        key_list = list(self.ds.data_vars)
        for arg in args:

            if not arg in key_list:
                raise KeyError(
                    "Quantity arguments must correspond to "
                    "DataArrays in the DataSet. "
                    f"{arg} was not found. "
                    f"Existing keys are: {key_list}"
                )

        # Interpolate the nearest-neighbor indices
        i = self.interpolate_indices(pos)

        # Get the indices that are equal to nan (fill values), then set
        # their values to 0. They will be over-written after the interpolation

        # Nan array is shape [n] and is 1 if none of the indices for a
        # position are NaN, and 0 otherwise.

        # i has different shape for non-uniform grids
        if self.is_uniform_grid:
            nan_mask = np.where(np.isnan(np.sum(i, axis=1)), 0, 1)
        else:
            nan_mask = np.where(np.isnan(i), 0, 1)

        # Replace all NaNs temporarily with 0
        i = np.where(np.isnan(i), 0, i)
        i = i.astype(np.int32)  # Cast as integers

        # Fetch the values at those indices from each quantity
        output = []
        for arg in args:
            # Read the values from the dataset
            if self.is_uniform_grid:
                values = self.ds[arg].values[i[:, 0], i[:, 1], i[:, 2]]
            else:
                values = self.ds[arg].values[i]
            # Apply the NaN mask (set any values that were out of bounds
            # to zero)
            values *= nan_mask
            values = np.squeeze(values)
            values *= self.ds[arg].attrs["unit"]
            output.append(values)

        if len(output) == 1:
            return output[0]
        else:
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

        *args : str
            Strings that correspond to DataArrays in the dataset

        """

        raise NotImplementedError(
            "Volume-averaged interpolator is not yet " "implemented for this grid type."
        )


class CartesianGrid(AbstractGrid):
    r"""
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

        *args : str
            Strings that correspond to DataArrays in the dataset

        """

        # Condition pos
        # If a single point was given, add empty dimension
        if pos.ndim == 1:
            pos = np.reshape(pos, [1, 3])

        # Convert position to u.Quantiy
        if not hasattr(pos, "unit"):
            pos *= self.unit

        # Validate args
        # must be np.ndarray or u.Quantity arrays of same shape as grid
        key_list = list(self.ds.data_vars)
        for arg in args:

            if not arg in key_list:
                raise KeyError(
                    "Quantity arguments must correspond to "
                    "DataArrays in the DataSet. "
                    f"{arg} was not found. "
                    f"Existing keys are: {key_list}"
                )

        # Interpolate the indices
        i = self.interpolate_indices(pos)
        nparticles = i.shape[0]

        # Get the indices that are equal to nan (fill values), then set
        # their values to 0. They will be over-written after the interpolation

        # Nan array is shape [n] and is 1 if none of the indices for a
        # position are NaN, and 0 otherwise.
        nan_mask = np.where(np.isnan(np.sum(i, axis=1)), 0, 1)
        # Replace all NaNs temporarily with 0
        i = np.where(np.isnan(i), 0, i)
        i = i.astype(np.int32)  # Cast as integers

        # Calculate the grid positions for each particle as interpolated
        # by the nearest neighbor interpolator
        xpos = self.ax0[i[:, 0]]
        ypos = self.ax1[i[:, 1]]
        zpos = self.ax2[i[:, 2]]

        # Determine the points bounding the grid cell containing the
        # particle
        x0 = np.where(pos[:, 0] > xpos, i[:, 0], i[:, 0] - 1)
        x1 = x0 + 1
        y0 = np.where(pos[:, 1] > ypos, i[:, 1], i[:, 1] - 1)
        y1 = y0 + 1
        z0 = np.where(pos[:, 2] > zpos, i[:, 2], i[:, 2] - 1)
        z1 = z0 + 1

        # Calculate the cell volume
        cell_vol = self.dax0 * self.dax1 * self.dax2
        n0, n1, n2 = self.shape

        # Create a list of empty arrays to hold final results
        output = []
        for i, arg in enumerate(args):
            output.append(np.zeros([nparticles]) * self.ds[arg].attrs["unit"])

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
                    grid_pos = (
                        np.array([self.ax0[x], self.ax1[y], self.ax2[z]]) * self.unit
                    )
                    grid_pos = np.moveaxis(grid_pos, 0, -1)

                    d = np.abs(grid_pos - pos)

                    # Fraction of cell volume that is closest to the
                    # current point
                    weight = (d[:, 0] * d[:, 1] * d[:, 2]) / cell_vol
                    weight = weight.to(u.dimensionless_unscaled)
                    weight[out] = 0

                    # For each argument, include the contributed by this
                    # grid vertex
                    for i, arg in enumerate(args):
                        values = self.ds[arg].values[x, y, z]
                        # Apply nan_mask to set out-of-bounds values to 0
                        values *= nan_mask

                        values *= self.ds[arg].attrs["unit"]
                        output[i] += weight * values

        if len(output) == 1:
            return output[0]
        else:
            return tuple(output)


class NonUniformCartesianGrid(CartesianGrid):
    r"""
    A Cartesian grid in which the _make_mesh method produces a non-uniformly
    spaced grid.
    """

    def _make_mesh(self, start, stop, num, **kwargs):
        r"""
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


def example_grid(name, L=1 * u.cm, num=100):

    grid = CartesianGrid(-L, L, num=num)

    if name == "electrostatic_gaussian_sphere":
        a = L / 3
        b = L / 2
        radius = np.linalg.norm(grid.grid * grid.unit, axis=3)
        arg = (radius / a).to(u.dimensionless_unscaled)
        potential = np.exp(-(arg ** 2)) * u.V

        Ex, Ey, Ez = np.gradient(potential, grid.dax0, grid.dax1, grid.dax2)

        Ex = -np.where(radius < b, Ex, 0)
        Ey = -np.where(radius < b, Ey, 0)
        Ez = -np.where(radius < b, Ez, 0)

        grid.add_quantities(["E_x", "E_y", "E_z"], [Ex, Ey, Ez])

    return grid
