import h5py

from plasmapy.classes import GenericPlasma
from plasmapy.utils import OpenPMDError

import os


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
        self.h5 = h5

    @property
    def electric_field(self):
        raise NotImplementedError

    @property
    def charge_density(self):
        raise NotImplementedError

    @property
    def x(self):
        raise NotImplementedError

    @property
    def y(self):
        raise AttributeError

    @property
    def z(self):
        raise NotImplementedError

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
