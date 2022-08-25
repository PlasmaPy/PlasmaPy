import astropy.units as u
import numpy as np
import pytest

import plasmapy.plasma

from plasmapy.particles.data.test import data_dir


@pytest.fixture(scope="module")
def h5(request):
    h5 = plasmapy.plasma.Plasma(hdf5=data_dir / "data00000255.h5")
    yield h5
    h5.close()


class TestPlasma:
    def test_patters(self):
        # Input data whose specific subclass cannot be known
        generic = plasmapy.plasma.Plasma(blablobleh="spam")
        assert isinstance(generic, plasmapy.plasma.GenericPlasma)

    def test_Plasma3D(self):
        # Input data for Plasma3D
        three_dimensional = plasmapy.plasma.Plasma(
            domain_x=np.linspace(0, 1, 3) * u.m,
            domain_y=np.linspace(0, 1, 3) * u.m,
            domain_z=np.linspace(0, 1, 3) * u.m,
        )
        assert isinstance(three_dimensional, plasmapy.plasma.sources.Plasma3D)
        assert isinstance(three_dimensional, plasmapy.plasma.BasePlasma)

    def test_PlasmaBlob(self):
        # Input data for PlasmaBlob
        T_e = 25 * 15e3 * u.K
        n_e = 1e26 * u.cm**-3
        Z = 2.0
        particle = "p"
        blob = plasmapy.plasma.Plasma(T_e=T_e, n_e=n_e, Z=Z, particle=particle)
        assert isinstance(blob, plasmapy.plasma.sources.PlasmaBlob)
        assert isinstance(blob, plasmapy.plasma.BasePlasma)
