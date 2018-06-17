# import h5py
# import openPMD

from plasmapy.classes import GenericPlasma


class HDF5Reader(GenericPlasma):
    def __init__(self, hdf5):
        # Load HDF5 using h5py or openPMD library
        # and do other stuff here.

        self.hdf5 = hdf5

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
        if len(kwargs) == 1:
            match = kwargs.get('hdf5')
        else:
            match = False
        return match
