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
from cached_property import cached_property
from collections import namedtuple
from scipy.spatial import distance
from typing import List, Union

from plasmapy.utils.decorators.helpers import modify_docstring


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
    Abstract grid represents a 3D grid of positions. The grid is stored
    as an `~numpy.ndarray`, while the units associated with each
    dimension are stored separately.


    There are two preferred methods to creating a grid object:

    1. Initializing the grid by providing three 3D
       `~astropy.units.Quantity` arrays of positions along each axis
       (``xpoints``, ``ypoints``, ``zpoints``)

       .. code-block:: python

          AbstractGrid(xpoints, ypoints, zpoints)

    2. A new grid can also be created using a syntax similar to
       `numpy.linspace` by providing two three-element
       `~astropy.units.Quantity` arrays of start and stop values and
       setting the ``num`` keyword to the number of points along each axis.

       .. code-block:: python

          AbstractGrid(
              start=[x0, y0, z0],
              stop=[x1, y1, z1],
              num=[Nx, Ny, Nz],
              **kwargs
          )

       In this case, any additional keyword arguments ``**kwargs`` provided
       will be passed directly to `~numpy.linspace`.

    """

    def __init__(self, *seeds, num=100, **kwargs):

        # Initialize some variables
        self._interpolator = None
        self._is_uniform = None

        # If three inputs are given, assume it's a user-provided grid
        if len(seeds) == 3:
            self._load_grid(seeds[0], seeds[1], seeds[2])

        # If two inputs are given, assume they are start and stop arrays
        # to create a new grid
        # kwargs are passed to np.linspace in _make_grid()
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
        consistent with the coordinate system and units selected.
        """
        return True

    # A named tuple describing a key recognized by PlasmaPy to correspond to
    # a particular physical quantity
    RecognizedQuantity = namedtuple(
        "RecognizedQuantities", ["key", "description", "unit"]
    )

    # These standard keys are used to refer to certain
    # physical quantities. This dictionary also provides the expected unit.
    _recognized_quantities_list = [
        RecognizedQuantity("x", "x spatial position", u.m),
        RecognizedQuantity("y", "y spatial position", u.m),
        RecognizedQuantity("z", "z spatial position", u.m),
        RecognizedQuantity("rho", "Mass density", u.kg / u.m ** 3),
        RecognizedQuantity("E_x", "Electric field (x component)", u.V / u.m),
        RecognizedQuantity("E_y", "Electric field (y component)", u.V / u.m),
        RecognizedQuantity("E_z", "Electric field (z component)", u.V / u.m),
        RecognizedQuantity("B_x", "Magnetic field (x component)", u.T),
        RecognizedQuantity("B_y", "Magnetic field (y component)", u.T),
        RecognizedQuantity("B_z", "Magnetic field (z component)", u.T),
        RecognizedQuantity("phi", "Electric Scalar Potential", u.V),
    ]

    # Create a dict of recognized quantities for fast access by key
    _recognized_quantities = {}
    for _rq in _recognized_quantities_list:
        _recognized_quantities[_rq.key] = _rq

    @property
    def recognized_quantities(self):
        r"""
        A dictionary of standard key names representing particular physical
        quantities. Using these keys allows these
        quantities to be recognized automatically by other PlasmaPy functions.
        Each entry contains a tuple containing a description and the unit
        associated with the quantity.
        """
        return self._recognized_quantities

    def require_quantities(self, req_quantities, replace_with_zeros=False):
        r"""
        Check to make sure that a list of required quantities are present.
        Optionally, can create missing quantities and fill them with
        an array of zeros.

        Parameters
        ----------
        req_quantities : `list` of `str`
            A list of quantity keys that are required.

        replace_with_zeros : `bool`, optional
            If true, missing quantities will be replaced with an array
            of zeros. If false, an exception will be raised instead.
            The default is False.

        Raises
        ------
        KeyError
            If ``replace_with_zeros`` is `False` and a required quantity
            is missing.

        KeyError
            If ``replace_with_zeros`` is `True` but the
            `~astropy.units.Quantity` is not in the list of recognized
            quantities. In this case the units for the quantity are unknown,
            so an array of zeros cannot be constructed.
        """
        for rq in req_quantities:

            # Error check that grid contains E and B variables required
            if rq not in self.quantities:

                # If missing, warn user and then replace with an array of zeros
                if replace_with_zeros:
                    warnings.warn(
                        f"{rq} is not specified for the provided grid."
                        "This quantity will be assumed to be zero.",
                        RuntimeWarning,
                    )

                    if rq in self.recognized_quantities.keys():
                        unit = self.recognized_quantities[rq].unit
                    else:
                        raise KeyError(
                            f"{rq} is not a recognized key, and "
                            "so cannot be automatically assumed "
                            "to be zero."
                        )

                    arg = {rq: np.zeros(self.shape) * unit}
                    self.add_quantities(**arg)

                else:
                    raise KeyError(
                        f"{rq} is not specified for the provided "
                        "grid but is required."
                    )

    # *************************************************************************
    # Fundamental properties of the grid
    # *************************************************************************

    def __repr__(self):

        line_sep = "-----------------------------\n"
        shape = list(self.shape)
        coords = list(self.ds.coords.keys())
        ax_units = self.units
        ax_dtypes = [self.ds[i].dtype for i in coords]

        coord_lbls = [str(i) + ": " + str(j) for i, j in zip(coords, shape)]

        s = f"*** Grid Summary ***\n{type(self)}\n"

        s += f"Dimensions: ({', '.join(coord_lbls)})\n"

        if self.is_uniform:
            s += (
                "Uniformly Spaced: (dax0, dax1, dax2) = "
                f"({self.dax0:.3f}, {self.dax1:.3f}, {self.dax2:.3f})\n"
            )
        else:
            s += "Non-Uniform Spacing\n"

        s += line_sep + "Coordinates:\n"
        for i in range(len(self.shape)):
            s += f"\t-> {coords[i]} ({ax_units[i]}) {ax_dtypes[i]} ({shape[i]},)\n"

        keys = self.quantities
        rkeys = [k for k in keys if k in list(self.recognized_quantities.keys())]
        nrkeys = [k for k in keys if k not in list(self.recognized_quantities.keys())]

        s += line_sep + "Recognized Quantities:\n"
        if len(rkeys) == 0:
            s += "-None-\n"
        else:
            for key in rkeys:
                unit = self.ds[key].attrs["unit"]
                dtype = self.ds[key].dtype
                shape = self.ds[key].shape
                s += f"\t-> {key} ({unit}) {dtype} {shape} \n"

        s += line_sep + "Unrecognized Quantities:\n"
        if len(nrkeys) == 0:
            s += "-None-\n"
        else:
            for key in nrkeys:
                unit = self.ds[key].attrs["unit"]
                dtype = self.ds[key].dtype
                shape = self.ds[key].shape
                s += f"\t-> {key} ({unit}) {dtype} {shape} \n"

        return s

    def __getitem__(self, key):
        """
        Given a key, return the corresponding array as a `~astropy.units.Quantity`

        Returning with ``copy=False`` means that the array returned is a direct
        reference to the underlying DataArray, so changes made will be reflected
        in the underlying DataArray.
        """
        return u.Quantity(self.ds[key].data, self.ds[key].attrs["unit"], copy=False)

    @property
    def is_uniform(self) -> bool:
        """
        A boolean value reflecting whether or not the grid points are
        uniformly spaced.
        """

        if self._is_uniform is None:  # coverage: ignore
            raise ValueError(
                "The `is_uniform` attribute is not accessible "
                "before a grid has been loaded."
            )

        return self._is_uniform

    @property
    def shape(self):
        r"""Shape of the grid"""
        if self.is_uniform:
            return (self.ax0.size, self.ax1.size, self.ax2.size)
        else:
            return self.ds.coords["ax0"].shape

    @property
    def grids(self):
        r"""
        Three grids of vertex positions (in each coordinate), each having
        shape (N0, N1, N2).
        """
        if self.is_uniform:
            pts0, pts1, pts2 = np.meshgrid(self.ax0, self.ax1, self.ax2, indexing="ij")
            _grids = (pts0, pts1, pts2)
        else:
            _grids = (
                self.ds["ax0"].data * self.unit0,
                self.ds["ax1"].data * self.unit1,
                self.ds["ax2"].data * self.unit2,
            )

        return _grids

    @property
    def grid(self):
        r"""
        A single grid of vertex positions of shape (N0, N1, N2, 3).

        Only defined for grids for which the ``unit`` property is defined.
        """
        pts0, pts1, pts2 = self.grids
        if self.is_uniform:
            n0, n1, n2 = pts0.shape
            grid = np.zeros([n0, n1, n2, 3]) * self.unit
        else:
            n = pts0.size
            grid = np.zeros([n, 3]) * self.unit

        grid[..., 0] = pts0
        grid[..., 1] = pts1
        grid[..., 2] = pts2
        _grid = grid

        return _grid

    @property
    def pts0(self):
        r"""Array of positions in dimension 1."""
        return self.grids[0]

    @property
    def pts1(self):
        r"""Array of positions in dimension 2."""
        return self.grids[1]

    @property
    def pts2(self):
        r"""Array of positions in dimension 3."""
        return self.grids[2]

    @property
    def units(self) -> list:
        r"""A `list` of the units of each dimension."""
        return self.ds.attrs["axis_units"]

    @property
    def unit0(self):
        r"""Unit of dimension 1."""
        return self.units[0]

    @property
    def unit1(self):
        r"""Unit of dimension 2."""
        return self.units[1]

    @property
    def unit2(self):
        r"""Unit of dimension 3."""
        return self.units[2]

    @property
    def unit(self):
        r"""
        The unit for the entire grid. Only valid if all dimensions of the
        grid have the same units.

        Raises
        ------
        `ValueError`
            If all grid dimensions do not have identical units.
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
    def si_scale_factors(self) -> List[float]:
        """
        3-element list containing unitless scale factors for converting
        the corresponding axis from its stored units to SI.
        """
        return self._si_factors

    def _get_ax(self, *, axis: int, si: bool = False):
        """
        Helper function for retrieving axis values.

        Parameters
        ----------
        axis: int
            Axis index for axis of interest (e.g. ``0`` for :attr:`ax0`).
        si: bool
            (Default: `False`) If `True` then convert the axis values to
            their SI equivalent.

        Returns
        -------
        ~numpy.ndarray or `~astropy.units.Quantity`
            If ``si==True`` then return a unitless `~numpy.ndarray`.
            If ``si==False`` then return a `~astropy.units.Quantity`
            array with the original units.

        Raises
        ------
        ValueError
            If the grid is not uniform.
        """
        ax_name = f"ax{axis}"

        if not self.is_uniform:
            raise ValueError(
                "The axis properties are only valid on uniformly spaced grids."
            )

        vals = self.ds.coords[ax_name].values
        if si:
            return vals * self.si_scale_factors[axis]
        else:
            return vals * self.units[axis]

    def _get_dax(self, *, axis: int, si: bool = False):
        """
        Helper function for calculating grid spacing.

        See Also
        --------
        plasmapy.plasma.grids.AbstractGrid._get_ax
        """
        ax = self._get_ax(axis=axis, si=si)
        return np.mean(np.gradient(ax))

    @property
    def _ax0_si(self):
        """
        The :attr:`ax0` axis without units, but scaled such that its values
        correspond to SI values.

        Only valid for uniform grids.
        """
        return self._get_ax(axis=0, si=True)

    @property
    def ax0(self):
        r"""
        First axis of the grid.

        Only valid for uniform grids.
        """
        return self._get_ax(axis=0)

    @property
    def _ax1_si(self):
        """
        The :attr:`ax1` axis without units, but scaled such that its values
        correspond to SI values.

        Only valid for uniform grids.
        """
        return self._get_ax(axis=1, si=True)

    @property
    def ax1(self):
        r"""
        Second axis of the grid.

        Only valid for uniform grids.
        """
        return self._get_ax(axis=1)

    @property
    def _ax2_si(self):
        """
        The :attr:`ax2` axis without units, but scaled such that its values
        correspond to SI values.

        Only valid for uniform grids.
        """
        return self._get_ax(axis=2, si=True)

    @property
    def ax2(self):
        r"""
        Third axis of the grid.

        Only valid for uniform grids.
        """
        return self._get_ax(axis=2)

    @property
    def _dax0_si(self):
        """
        Grid step size along axis :attr:`ax0` without units and scaled such
        that its values correspond to SI values.

        Only valid for uniform grids.
        """
        return self._get_dax(axis=0, si=True)

    @property
    def dax0(self):
        r"""
        Grid step size along axis :attr:`ax0`.

        Only valid for uniform grids.
        """
        return self._get_dax(axis=0)

    @property
    def _dax1_si(self):
        """
        Grid step size along axis :attr:`ax1` without units and scaled such
        that its values correspond to SI values.

        Only valid for uniform grids.
        """
        return self._get_dax(axis=1, si=True)

    @property
    def dax1(self):
        r"""
        Grid step size along axis :attr:`ax1`.

        Only valid for uniform grids.
        """
        return self._get_dax(axis=1)

    @property
    def _dax2_si(self):
        """
        Grid step size along axis :attr:`ax2` without units and scaled such
        that its values correspond to SI values.

        Only valid for uniform grids.
        """
        return self._get_dax(axis=2, si=True)

    @property
    def dax2(self):
        r"""
        Grid step size along axis :attr:`ax2`.

        Only valid for uniform grids.
        """
        return self._get_dax(axis=2)

    @property
    def grid_resolution(self):
        r"""
        A scalar estimate of the grid resolution.

        For uniform grids, this is the minima of [dax0, dax1, dax2].

        For non-uniform grids, it is the closest spacing between any two points.
        """

        if self.is_uniform:
            return min(self.dax0, self.dax1, self.dax2)
        else:
            distances = distance.cdist(self.grid, self.grid)
            np.fill_diagonal(distances, np.inf)
            return np.min(distances)

    # *************************************************************************
    # Loading and creating grids
    # *************************************************************************

    def _load_grid(
        self,
        pts0: u.Quantity,
        pts1: u.Quantity,
        pts2: u.Quantity,
    ):
        r"""
        Initialize the grid object from a user-supplied grid.

        Parameters
        ----------
        grid{0,1,2} : `~astropy.units.Quantity` array, shape (n0, n1, n2)
            Grids of coordinate positions.

        **kwargs: `~astropy.units.Quantity` array, shape (n0, n1, n2)
            Quantities defined on the grid.
        """

        # Validate input
        if not (pts0.shape == pts1.shape and pts0.shape == pts2.shape):
            raise ValueError(
                "Provided arrays of grid points are of unequal "
                f"shape: pts0 = {pts0.shape}, "
                f"pts1 = {pts1.shape}, "
                f"pts2 = {pts2.shape}."
            )

        self._is_uniform = _detect_is_uniform_grid(pts0, pts1, pts2)

        # Create dataset
        self.ds = xr.Dataset()

        self.ds.attrs["axis_units"] = [pts0.unit, pts1.unit, pts2.unit]

        # Store the conversion factors for each axis to SI
        self._si_factors = [
            pts0.unit.si.scale,
            pts1.unit.si.scale,
            pts2.unit.si.scale,
        ]

        if self.is_uniform:
            self.ds.coords["ax0"] = pts0[:, 0, 0]
            self.ds.coords["ax1"] = pts1[0, :, 0]
            self.ds.coords["ax2"] = pts2[0, 0, :]

        else:
            mdx = pd.MultiIndex.from_arrays(
                [pts0.flatten(), pts1.flatten(), pts2.flatten()],
                names=["ax0", "ax1", "ax2"],
            )
            self.ds.coords["ax"] = mdx

        # Check to make sure that the object created satisfies any
        # requirements: eg. units correspond to the coordinate system
        self._validate()

    def add_quantities(self, **kwargs):
        r"""
        Adds a quantity to the dataset as a new DataArray.

        Parameters
        ----------
        key, array pairs as keyword arguments
            The key will be used as the dataset key, while the array holds the
            quantity.
        """

        for key in kwargs.keys():
            quantity = kwargs[key]

            # Check key against a list of "known" keys with pre-defined
            # meanings (eg. E_x, n_e) and raise a warning if a "non-standard"
            # key is being used so the user is aware.
            if key in self.recognized_quantities.keys():
                try:
                    quantity.to(self.recognized_quantities[key].unit)
                except u.UnitConversionError:
                    raise ValueError(
                        f"Units provided for {key} ({quantity.unit}) "
                        "are not compatible with the correct units "
                        f"for that recognized key ({self.recognized_quantities[key]})."
                    )

            else:
                warnings.warn(
                    f"Warning: {key} is not recognized quantity key", stacklevel=2
                )

            if self.is_uniform:
                dims = ["ax0", "ax1", "ax2"]
                coords = {
                    "ax0": self.ds.coords["ax0"],
                    "ax1": self.ds.coords["ax1"],
                    "ax2": self.ds.coords["ax2"],
                }
            # If grid is non-uniform, flatten quantity
            else:
                quantity = quantity.flatten()
                dims = ["ax"]
                coords = {"ax": self.ds.coords["ax"]}

            if quantity.shape != self.shape:
                raise ValueError(
                    f"Shape of quantity '{key}' {quantity.shape} "
                    f"does not match the grid shape {self.shape}."
                )

            data = xr.DataArray(
                quantity, dims=dims, coords=coords, attrs={"unit": quantity.unit}
            )
            self.ds[key] = data

    @property
    def quantities(self):
        r"""
        A list of the keys corresponding to the quantities currently
        defined on the grid.
        """
        return list(self.ds.data_vars)

    def _make_grid(
        self,
        start: Union[int, float, u.Quantity],
        stop: Union[int, float, u.Quantity],
        num=100,
        units=None,
        **kwargs,
    ):
        r"""
        Creates a grid based on ``start``, ``stop``, and ``num`` values
        in a manner that mirrors the interface of the `numpy.linspace`
        function.

        Parameters
        ----------
        start : number (`~astropy.units.Quantity`) or array of three of the same
            Starting values for each dimension. If one value is given,
            the same value will be used for all three dimensions.

        stop : number (`~astropy.units..Quantity`) or array of three of the same
            End values for each dimension. If one value is given,
            the same value will be used for all three dimensions.

        num : `int` or `list` of three `int` objects, optional
            The number of points in each dimension. If a single integer is
            given, the same number of points will be used in each dimension.
            The default is 100.

        **kwargs: Additional arguments
            Any additional arguments will be passed directly to
            `numpy.linspace`.
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
            pts0 * units[0],
            pts1 * units[1],
            pts2 * units[2],
        )

    def _make_mesh(self, start, stop, num, **kwargs):
        r"""
        Creates mesh as part of _make_grid(). Separated into its own function
        so it can be re-implemented to make non-uniformly spaced meshes.
        """
        # Construct the axis arrays
        ax0 = np.linspace(start[0], stop[0], num=num[0], **kwargs)
        ax1 = np.linspace(start[1], stop[1], num=num[1], **kwargs)
        ax2 = np.linspace(start[2], stop[2], num=num[2], **kwargs)

        # Construct the coordinate arrays
        pts0, pts1, pts2 = np.meshgrid(ax0, ax1, ax2, indexing="ij")

        return pts0, pts1, pts2

    # *************************************************************************
    # Methods
    # *************************************************************************

    def on_grid(self, pos):
        r"""
        Given a list of positions, determines which are in the region
        bounded by the grid points.

        For non-uniform grids, "on grid" is defined as being bounded by
        grid points in all axes.

        Parameters
        ----------
        pos : `~numpy.ndarray` or `~astropy.units.Quantity` array, shape (n,3)
            An array of positions in space, where the second dimension
            corresponds to the three dimensions of the grid.

        """

        if hasattr(pos, "unit"):
            pos = pos.si.value

        # Find the bounds
        if self.is_uniform:
            ax0_min, ax0_max = np.min(self.ax0.si.value), np.max(self.ax0.si.value)
            ax1_min, ax1_max = np.min(self.ax1.si.value), np.max(self.ax1.si.value)
            ax2_min, ax2_max = np.min(self.ax2.si.value), np.max(self.ax2.si.value)

        else:
            pts0, pts1, pts2 = self.grids
            ax0_min, ax0_max = np.min(self.pts0).si.value, np.max(self.pts0).si.value
            ax1_min, ax1_max = np.min(self.pts1).si.value, np.max(self.pts1).si.value
            ax2_min, ax2_max = np.min(self.pts2).si.value, np.max(self.pts2).si.value

        # Check each point elementwise against the bounds
        on_grid = (
            np.greater(ax0_min, pos[:, 0]).astype(np.int8)
            + np.less(ax0_max, pos[:, 0]).astype(np.int8)
            + np.greater(ax1_min, pos[:, 1]).astype(np.int8)
            + np.less(ax1_max, pos[:, 1]).astype(np.int8)
            + np.greater(ax2_min, pos[:, 2]).astype(np.int8)
            + np.less(ax2_max, pos[:, 2]).astype(np.int8)
        )

        return np.where(on_grid == 0, True, False)

    def vector_intersects(self, p1, p2):
        r"""
        `True` if the vector from ``p1`` to ``p2`` intersects the grid,
        and `False` otherwise.

        This is a standard ray-box intersection algorithm.
        """
        p1, p2 = p1.si.value, p2.si.value
        # Caclulate the minimum and maximum of each
        Ax, Bx = np.min(self.pts0.si.value), np.max(self.pts0.si.value)
        Ay, By = np.min(self.pts1.si.value), np.max(self.pts1.si.value)
        Az, Bz = np.min(self.pts2.si.value), np.max(self.pts2.si.value)
        A = np.array([Ax, Ay, Az])
        B = np.array([Bx, By, Bz])

        # Calculate the equation of the line from p1 to p2 such that
        # r = p1 + t*D
        D = np.abs(p2 - p1)

        # Calculate the intersection points. These operations are just vectorized
        # for convenience. Ignore div-by-zero: outputting infty's here is fine.
        with np.errstate(divide="ignore"):
            Tmin = (A - p1) / D
            Tmax = (B - p1) / D

        Tmin = np.max(Tmin)
        Tmax = np.min(Tmax)

        return Tmin < Tmax

    # *************************************************************************
    # Interpolators
    # *************************************************************************

    # This property holds the list of quantity keys currently being interpolated
    # It's used in the following cached properties
    _interp_args = []

    @cached_property
    def _interp_quantities(self):
        r"""Create a dimensionless array of quantities to be interpolated."""
        nargs = len(self._interp_args)
        # Load the arrays to be interpolated from and their units
        if self.is_uniform:
            nx, ny, nz = self.shape
            _interp_quantities = np.zeros([nx, ny, nz, nargs])
        else:
            npoints = self.shape[0]
            _interp_quantities = np.zeros([npoints, nargs])

        for j, arg in enumerate(self._interp_args):
            _interp_quantities[..., j] = self.ds[arg].values

        return _interp_quantities

    @cached_property
    def _interp_units(self):
        r"""
        Create a list of the units corresponding to the last dimension
        in the `_interp_quantities` array.
        """
        _interp_units = []
        for j, arg in enumerate(self._interp_args):
            _interp_units.append(self.ds[arg].attrs["unit"])

        return _interp_units

    @property
    def interpolator(self):
        r"""
        A nearest-neighbor interpolator that returns the nearest grid index
        to a position.
        """
        if self._interpolator is None:
            if self.is_uniform:
                self._make_uniform_grid_interpolator()
            else:
                self._make_nonuniform_grid_interpolator()

        return self._interpolator

    def _make_uniform_grid_interpolator(self):
        r"""
        Initializes a nearest-neighbor interpolator that returns the nearest
        grid indices for a given position (given in SI units).

        This function works on a uniformly spaced grid.
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

        This function works on unstructured (non-uniform) data.
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
        nearest-neighbor interpolator.

        Parameters
        ----------
        pos : `~numpy.ndarray` or `~astropy.units.Quantity` array, shape (n,3)
            An array of positions in space, where the second dimension
            corresponds to the three dimensions of the grid. If a
            `~numpy.ndarray` is provided, units will be assumed to match
            those of the grid.

        Returns
        -------
        i : `~numpy.ndarray`, shape (n,3)
            An array of indices corresponding to the positions such that
            ``i[n,:] = ix,iy,iz`` such that ``grid[ix,iy,iz,:]`` ∼ ``pos[n,:]``

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

    def nearest_neighbor_interpolator(
        self, pos: Union[np.ndarray, u.Quantity], *args, persistent=False
    ):
        r"""
        Interpolate values on the grid using a nearest-neighbor scheme with
        no higher-order weighting.

        Parameters
        ----------
        pos : `~numpy.ndarray` or `~astropy.units.Quantity` array, shape (n,3)
            An array of positions in space, where the second dimension
            corresponds to the three dimensions of the grid. If an
            `~numpy.ndarray` is provided, units will be assumed to match
            those of the grid.

        *args : `str`
            Strings that correspond to DataArrays in the dataset

        persistent : `bool`
            If `True`, the interpolator will assume the grid and its
            contents have not changed since the last interpolation. This
            substantially speeds up the interpolation when many
            interpolations are performed on the same grid in a loop.
            ``persistent`` overrides to `False` if the arguments list
            has changed since the last call.
        """
        # pos is validated in interpolate_indices

        # Validate args
        # must be np.ndarray or u.Quantity arrays of same shape as grid
        for arg in args:

            if arg not in self.quantities:
                raise KeyError(
                    "Quantity arguments must correspond to "
                    "DataArrays in the DataSet. "
                    f"{arg} was not found. "
                    f"Existing keys are: {self.quantities}"
                )

        # If persistent, double check the arguments list hasn't changed
        # If they have, run as non-persistent this time
        if persistent and args != self._interp_args:
            persistent = False

        # Update _interp_args variable
        self._interp_args = args

        # Interpolate the nearest-neighbor indices
        i = self.interpolate_indices(pos)
        nargs = len(args)

        # Get the indices that are equal to nan (fill values), then set
        # their values to 0. They will be over-written after the interpolation

        # Nan array is shape [n] and is 1 if none of the indices for a
        # position are NaN, and 0 otherwise.

        # i has different shape for non-uniform grids
        if self.is_uniform:
            nan_mask = np.where(np.isnan(np.sum(i, axis=1)), 0, 1)
        else:
            nan_mask = np.where(np.isnan(i), 0, 1)

        # Replace all NaNs temporarily with 0
        i = np.where(np.isnan(i), 0, i)
        i = i.astype(np.int32)  # Cast as integers

        # If not persistent, clear the cached properties so they are re-created
        # when called below
        if not persistent:
            try:
                del self._interp_quantities
            except AttributeError:
                pass
            try:
                del self._interp_units
            except AttributeError:
                pass

        # Fetch the values at those indices from each quantity
        if self.is_uniform:
            values = self._interp_quantities[i[:, 0], i[:, 1], i[:, 2], :]
        else:
            values = self._interp_quantities[i]

        # Apply the NaN mask (set any values that were out of bounds
        # to zero)
        values *= np.outer(nan_mask, np.ones(nargs))

        # Split output array into arrays with units
        # Apply units to output arrays
        output = []
        for i in range(nargs):
            output.append(values[:, i] * self._interp_units[i])

        if len(output) == 1:
            return output[0]
        else:
            return tuple(output)

    def volume_averaged_interpolator(
        self, pos: Union[np.ndarray, u.Quantity], *args, persistent=False
    ):
        r"""
        Interpolate values on the grid using a volume-averaged scheme with
        no higher-order weighting.

        Parameters
        ----------
        pos : `~numpy.ndarray` or `~astropy.units.Quantity` array, shape (n,3)
            An array of positions in space, where the second dimension
            corresponds to the three dimensions of the grid. If a
            `~numpy.ndarray` is provided, units will be assumed to match
            those of the grid.

        *args : `str`
            Strings that correspond to DataArrays in the dataset

        persistent : `bool`
            If `True`, the interpolator will assume the grid and its
            contents have not changed since the last interpolation. This
            substantially speeds up the interpolation when many
            interpolations are performed on the same grid in a loop.
            ``persistent`` overrides to `False` if the arguments list
            has changed since the last call.
        """

        raise NotImplementedError(
            "Volume-averaged interpolator is not yet " "implemented for this grid type."
        )


class CartesianGrid(AbstractGrid):
    r"""A uniformly spaced Cartesian grid."""

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

    @modify_docstring(prepend=AbstractGrid.volume_averaged_interpolator.__doc__)
    def volume_averaged_interpolator(
        self, pos: Union[np.ndarray, u.Quantity], *args, persistent=False
    ):
        """
        Notes
        -----

        This interpolator approximates the value of a quantity at a given
        interpolation point using a weighted sum of the values at the eight grid
        vertices that surround the point. The weighting factors are calculated by
        defining a volume :math:`dx \\times dy \\times dz`
        (where :math:`dx`, :math:`dy`, and :math:`dz` are the grid
        spacings in each direction) around each grid vertex and around the
        interpolation point. The contribution of each grid vertex is then
        weighted by the fraction of the volume surrounding the interpolation
        point that overlaps the volume surrounding that vertex. This effectively
        introduces a linear interpolation between grid vertices.

        This implementation of this algorithm assumes that the grid is uniformly
        spaced and Cartesian.
        """

        if isinstance(pos, u.Quantity):
            pos = pos.to(u.m).value
        elif self.unit != u.m:
            pos *= self.unit.si.scale

        # If a single point was given, add empty dimension
        if pos.ndim == 1:
            pos = np.reshape(pos, [1, 3])

        # -- Validate args --
        # must be np.ndarray or u.Quantity arrays of same shape as grid
        for arg in args:
            if arg not in self.quantities:
                raise KeyError(
                    "Quantity arguments must correspond to "
                    "DataArrays in the DataSet. "
                    f"{arg} was not found. "
                    f"Existing keys are: {self.quantities}"
                )

        # If persistent, double check the arguments list hasn't changed
        # If they have, run as non-persistent this time
        if persistent and args != self._interp_args:
            persistent = False

        # Update _interp_args variable
        self._interp_args = args

        # If not persistent, clear the cached properties so they are re-created
        # when called below
        if not persistent:
            try:
                del self._interp_quantities
            except AttributeError:
                pass
            try:
                del self._interp_units
            except AttributeError:
                pass

        # -- begin averaging --

        nparticles = pos.shape[0]
        nargs = len(args)

        # Load grid attributes (so this isn't repeated)
        ax0, ax1, ax2 = self._ax0_si, self._ax1_si, self._ax2_si
        dx, dy, dz = self._dax0_si, self._dax1_si, self._dax2_si
        n0, n1, n2 = self.shape

        # find cell nearest to each position
        nearest_neighbor_index = np.zeros((nparticles, 3), dtype=np.int32)
        nearest_neighbor_index[..., 0] = np.abs(pos[:, 0, None] - ax0).argmin(axis=1)
        nearest_neighbor_index[..., 1] = np.abs(pos[:, 1, None] - ax1).argmin(axis=1)
        nearest_neighbor_index[..., 2] = np.abs(pos[:, 2, None] - ax2).argmin(axis=1)

        # Create a mask for positions that are off the grid. The values at
        # these points will be set to zero later.
        mask_particle_off = (
            (pos[:, 0] < ax0.min())
            | (pos[:, 0] > ax0.max())
            | (pos[:, 1] < ax1.min())
            | (pos[:, 1] > ax1.max())
            | (pos[:, 2] < ax2.min())
            | (pos[:, 2] > ax2.max())
        )

        # Get the physical positions for the nearest neighbor cell
        # for the each particle
        xpos = ax0[nearest_neighbor_index[:, 0]]
        ypos = ax1[nearest_neighbor_index[:, 1]]
        zpos = ax2[nearest_neighbor_index[:, 2]]
        nearest_neighbor_pos = np.array([xpos, ypos, zpos]).swapaxes(0, 1)

        # Determine the indices for the grid cells bounding the particle
        # - The imaginary cell centered on the particle will overlap with
        #   1 to 8 surrounding cells
        # - typically this is 8 but will be 4, 2, or 1 when the particle is
        #   near the boundary of the grid
        bounding_cell_indices = np.empty((nparticles, 8, 3), dtype=np.int32)
        lower_indices = np.where(
            pos >= nearest_neighbor_pos,
            nearest_neighbor_index,
            nearest_neighbor_index - 1,
        )

        # populate x indices
        bounding_cell_indices[:, 0:4, 0] = np.tile(
            lower_indices[:, 0], (4, 1)
        ).swapaxes(0, 1)
        bounding_cell_indices[:, 4:, 0] = bounding_cell_indices[:, 0:4, 0] + 1

        # populate y indices
        bounding_cell_indices[:, [0, 1, 4, 5], 1] = np.tile(
            lower_indices[:, 1], (4, 1)
        ).swapaxes(0, 1)
        bounding_cell_indices[:, [2, 3, 6, 7], 1] = (
            bounding_cell_indices[:, [0, 1, 4, 5], 1] + 1
        )

        # populate z indices
        bounding_cell_indices[:, 0::2, 2] = np.tile(
            lower_indices[:, 2], (4, 1)
        ).swapaxes(0, 1)
        bounding_cell_indices[:, 1::2, 2] = bounding_cell_indices[:, 0::2, 2] + 1

        # Create a mask for cells whose locations would be off the grid,
        # which occurs when the interpolation point is near the edge of the
        # grid. These values will be weighted as zero later.
        mask_cell_off = (
            (bounding_cell_indices < 0).any(axis=2)
            | (bounding_cell_indices[:, :, 0] >= n0)
            | (bounding_cell_indices[:, :, 1] >= n1)
            | (bounding_cell_indices[:, :, 2] >= n2)
        )

        # Zero any out of bounds indices so IndexError is not raised
        # during indexing.  This means an incorrect value will be retrieved
        # but will not be used because of the zero weighting and the
        # off the grid mask
        bounding_cell_indices[mask_cell_off, :] = 0

        # Calculate the volume of the overlap between the point volume
        # and the volume of each of the surrounding vertices
        lx = dx - np.abs(pos[:, None, 0] - ax0[bounding_cell_indices[..., 0]])
        ly = dy - np.abs(pos[:, None, 1] - ax1[bounding_cell_indices[..., 1]])
        lz = dz - np.abs(pos[:, None, 2] - ax2[bounding_cell_indices[..., 2]])
        bounding_cell_weights = lx * ly * lz

        # Set the weight for any off-grid vertices (cell or particle) to zero
        bounding_cell_weights[mask_cell_off] = 0.0
        bounding_cell_weights[mask_particle_off, ...] = 0.0
        norms = np.sum(bounding_cell_weights, axis=1)
        mask_norm_zero = norms == 0.0
        bounding_cell_weights[~mask_norm_zero] = (
            bounding_cell_weights[~mask_norm_zero, ...] / norms[~mask_norm_zero, None]
        )

        # Get the values of each of the interpolated quantities at each
        # of the bounding vertices
        vals = self._interp_quantities[
            bounding_cell_indices[..., 0],
            bounding_cell_indices[..., 1],
            bounding_cell_indices[..., 2],
            :,
        ]
        # Construct a weighted average of the interpolated quantities
        weighted_ave = np.sum(bounding_cell_weights[..., None] * vals, axis=1)
        weighted_ave[mask_particle_off, :] = np.nan

        # Split output array into arrays with units
        # Apply units to output arrays
        output = []
        for arg in range(nargs):
            output.append(weighted_ave[..., arg] * self._interp_units[arg])

        if len(output) == 1:
            return output[0]
        else:
            return tuple(output)


class NonUniformCartesianGrid(AbstractGrid):
    r"""
    A Cartesian grid in which the ``_make_mesh`` method produces a
    non-uniformly spaced grid.
    """

    def _validate(self):
        """Check that all units are lengths."""
        for i in range(3):
            try:
                self.units[i].to(u.m)
            except u.UnitConversionError:
                raise ValueError(
                    "Units of grid are not valid for a Cartesian "
                    f"grid: {self.units}."
                )

    def _make_mesh(self, start, stop, num, **kwargs):
        r"""
        Creates mesh as part of ``_make_grid()``. Separated into its own
        function so it can be re-implemented to make non-uniform grids.
        """
        # Construct the axis arrays
        ax0 = np.sort(np.random.uniform(low=start[0], high=stop[0], size=num[0]))

        ax1 = np.sort(np.random.uniform(low=start[1], high=stop[1], size=num[1]))

        ax2 = np.sort(np.random.uniform(low=start[2], high=stop[2], size=num[2]))

        # Construct the coordinate arrays
        arr0, arr1, arr2 = np.meshgrid(ax0, ax1, ax2, indexing="ij")

        return arr0, arr1, arr2
