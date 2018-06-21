import h5py
import numpy as np
import astropy.units as u

from plasmapy.classes import GenericPlasma
from plasmapy.utils import OpenPMDError

import os


# This is the order what OpenPMD uses to store unit
# dimensions for a record.
_UNITS = (u.meter,
          u.kilogram,
          u.second,
          u.ampere,
          u.Kelvin,
          u.mol,
          u.candela)


def _fetch_units(openPMD_dims):
    """
    Converts a collection of OpenPMD dimensions to astropy.units.
    """

    units = u.dimensionless_unscaled
    for factor, unit in zip(openPMD_dims, _UNITS):
        units *= (unit ** factor)
    units, *_ = units.compose()
    return units


class HDF5Reader(GenericPlasma):
    def __init__(self, hdf5, **kwargs):
        """
        Core class for accessing various attributes on HDF5 files that
        are based on OpenPMD standards.

        Attributes
        ----------
        electric_field : `astropy.units.Quantity`
            An (x, y, z) array containing electric field data.
        charge_density : `astropy.units.Quantity`
            An array containing charge density data.

        Parameters
        ----------
        hdf5 : `str`
            Path to HDF5 file.
        """

        if not os.path.isfile(hdf5):
            raise FileNotFoundError(f"Could not find file: '{hdf5}'")

        h5 = h5py.File(hdf5)
        try:
            openPMD = h5.attrs["openPMDextension"]
        except KeyError:
            openPMD = False

        if not openPMD:
            raise OpenPMDError("Input HDF5 file does not go on with "
                               "standards defined by OpenPMD")

        self.subname = tuple(h5['data'])[0]
        self.h5 = h5

    @property
    def electric_field(self):
        path = 'data/' + self.subname + '/fields/E'
        if path in self.h5:
            units = _fetch_units(self.h5[path].attrs["unitDimension"])
            return np.array((self.h5[path]['x'],
                             self.h5[path]['y'],
                             self.h5[path]['z'])) * units
        else:
            raise AttributeError("No electric field data available "
                                 "in HDF5 file")

    @property
    def charge_density(self):
        path = 'data/' + self.subname + '/fields/rho'
        if path in self.h5:
            units = _fetch_units(self.h5[path].attrs["unitDimension"])
            return np.array(self.h5[path]) * units
        else:
            raise AttributeError("No charge density data available "
                                 "in HDF5 file")

    @classmethod
    def is_datasource_for(cls, **kwargs):
        if not "hdf5" in kwargs:
            return False

        hdf5 = kwargs.get("hdf5")
        openPMD = kwargs.get("openPMD")

        isfile = os.path.isfile(hdf5)
        if not isfile:
            raise FileNotFoundError(f"Could not find file: '{hdf5}'")

        if not "openPMD" in kwargs and isfile:
            h5 = h5py.File(hdf5)
            try:
                openPMD = h5.attrs["openPMDextension"]
            except KeyError:
                openPMD = False

        return openPMD
