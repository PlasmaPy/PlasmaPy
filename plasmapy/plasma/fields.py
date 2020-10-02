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

import astropy.units as u
import numpy as np

from abc import ABC, abstractmethod
from scipy.special import erf as erf



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
