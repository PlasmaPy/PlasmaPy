import h5py
import numpy as np
import astropy.units as u

from plasmapy.classes import GenericPlasma
from plasmapy.utils import OpenPMDError

import os


_UNITS = (u.meter,
          u.kilogram,
          u.second,
          u.ampere,
          u.Kelvin,
          u.mol,
          u.candela)


def _fetch_units(openPMD_dims):
    units = u.dimensionless_unscaled
    for factor, unit in zip(openPMD_dims, _UNITS):
        units *= (unit ** factor)
    units, *_ = units.compose()
    return units


class _ElectricField:
    def __init__(self, E):
        self.units = _fetch_units(E.attrs["unitDimension"])
        self.E = E

    @property
    def x(self):
        return np.array(self.E['x']) * self.units

    @property
    def y(self):
        return np.array(self.E['y']) * self.units

    @property
    def z(self):
        return np.array(self.E['z']) * self.units


class HDF5Reader(GenericPlasma):
    def __init__(self, hdf5, **kwargs):
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
            return _ElectricField(self.h5[path])
        else:
            raise AttributeError

    @property
    def charge_density(self):
        path = 'data/' + self.subname + '/fields/rho'
        if path in self.h5:
            units = _fetch_units(self.h5[path].attrs["unitDimension"])
            return np.array(self.h5[path]) * units
        else:
            raise AttributeError

    @classmethod
    def is_datasource_for(cls, **kwargs):
        if not "hdf5" in kwargs:
            return False

        hdf5 = kwargs.get("hdf5")
        openPMD = kwargs.get("openPMD")

        if not "openPMD" in kwargs and os.path.isfile(hdf5):
            h5 = h5py.File(hdf5)
            try:
                openPMD = h5.attrs["openPMDextension"]
            except KeyError:
                openPMD = False

        return openPMD
