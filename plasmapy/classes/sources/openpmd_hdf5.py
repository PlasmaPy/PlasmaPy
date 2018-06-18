import h5py
import numpy as np
import astropy.units as u

from plasmapy.classes import GenericPlasma
from plasmapy.utils import OpenPMDError

import os


class _ElectricField:
    def __init__(self, E):
        self.E = E

    @property
    def x(self):
        return np.array(self.E['x']) * u.V / u.m

    @property
    def y(self):
        return np.array(self.E['y']) * u.V / u.m

    @property
    def z(self):
        return np.array(self.E['z']) * u.V / u.m


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
            return np.array(self.h5[path]) * u.C / u.m**3
        else:
            raise AttributeError

    @property
    def x(self):
        path = 'data/' + self.subname + '/fields/E/x'
        if path in self.h5:
            return self.h5[path]
        else:
            raise AttributeError

    @property
    def y(self):
        raise AttributeError

    @property
    def z(self):
        path = 'data/' + self.subname + '/fields/E/z'
        if path in self.h5:
            return self.h5[path]
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
