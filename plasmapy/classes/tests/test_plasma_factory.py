import astropy.units as u
import numpy as np
import os

import plasmapy.classes


class TestPlasma(object):
    def test_patters(self):
        # Input data whose specific subclass cannot be known
        generic = plasmapy.classes.Plasma(blablobleh='spam')
        assert isinstance(generic, plasmapy.classes.GenericPlasma)

    def test_Plasma3D(self):
        # Input data for Plasma3D
        three_dimensional = plasmapy.classes.Plasma(domain_x=np.linspace(0, 1, 3) * u.m,
                                                    domain_y=np.linspace(0, 1, 3) * u.m,
                                                    domain_z=np.linspace(0, 1, 3) * u.m)
        assert isinstance(three_dimensional, plasmapy.classes.sources.Plasma3D)
        assert isinstance(three_dimensional, plasmapy.classes.BasePlasma)

    def test_PlasmaBlob(self):
        # Input data for PlasmaBlob
        T_e = 25 * 15e3 * u.K
        n_e = 1e26 * u.cm ** -3
        Z = 2.0
        particle = 'p'
        blob = plasmapy.classes.Plasma(T_e=T_e,
                                       n_e=n_e,
                                       Z=Z,
                                       particle=particle)
        assert isinstance(blob, plasmapy.classes.sources.PlasmaBlob)
        assert isinstance(blob, plasmapy.classes.BasePlasma)

    def test_HDF5Reader(self):
        h5 = plasmapy.classes.Plasma(hdf5=os.path.join(
                                         os.path.dirname(__file__) or '.',
                                         "../sources/tests/data",
                                         "data00000255.h5"))
        assert isinstance(h5, plasmapy.classes.sources.HDF5Reader)
        assert isinstance(h5, plasmapy.classes.BasePlasma)
