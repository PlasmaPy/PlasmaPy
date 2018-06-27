from plasmapy.classes.sources import openpmd_hdf5
from plasmapy.utils import OpenPMDError
from plasmapy.data.test import rootdir

from astropy import units as u
from typing import Union, Tuple, List
import os
import pytest


class TestOpenPMD2D:
    """Test 2D HDF5 dataset based on OpenPMD."""
    # Downloaded from
    # https://github.com/openPMD/openPMD-example-datasets/blob/draft/example-2d.tar.gz
    # per the Creative Commons Zero v1.0 Universal license
    h5 = openpmd_hdf5.HDF5Reader(hdf5=os.path.join(rootdir, "data00000255.h5"))

    def test_has_electric_field_with_units(self):
        assert self.h5.electric_field.to(u.V / u.m)

    def test_correct_shape_electric_field(self):
        assert self.h5.electric_field.shape == (3, 51, 201)

    def test_has_charge_density_with_units(self):
        assert self.h5.charge_density.to(u.C / u.m**3)

    def test_correct_shape_charge_density(self):
        assert self.h5.charge_density.shape == (51, 201)


class TestOpenPMD3D:
    """Test 3D HDF5 dataset based on OpenPMD."""
    # Downloaded from
    # https://github.com/openPMD/openPMD-example-datasets/blob/draft/example-3d.tar.gz
    # per the Creative Commons Zero v1.0 Universal license
    h5 = openpmd_hdf5.HDF5Reader(hdf5=os.path.join(rootdir, "data00000100.h5"))

    def test_has_electric_field_with_units(self):
        assert self.h5.electric_field.to(u.V / u.m)

    def test_correct_shape_electric_field(self):
        assert self.h5.electric_field.shape == (3, 26, 26, 201)

    def test_has_charge_density_with_units(self):
        assert self.h5.charge_density.to(u.C / u.m**3)

    def test_correct_shape_charge_density(self):
        assert self.h5.charge_density.shape == (26, 26, 201)


units_test_table = [
    ((1., 1., 0., -1., 0., 0., 2.),
     u.m * u.kg / u.amp * u.cd ** 2),
    ((1, 0, 1, 2, 0, 0, 0),
     u.m * u.s * u.amp ** 2),
    ([-3.,  0.,  1.,  1.,  0.,  0.,  0.],
     u.coulomb / u.m**3),
    ([2, 1, -3, -2, 0, 0, 0],
     u.ohm)
]


@pytest.mark.parametrize("openPMD_dims, expected", units_test_table)
def test_fetch_units(openPMD_dims, expected: Union[Tuple, List]):
    units = openpmd_hdf5._fetch_units(openPMD_dims)
    assert units == expected


def test_unavailable_hdf5():
    with pytest.raises(FileNotFoundError):
        openpmd_hdf5.HDF5Reader(hdf5="this_file_does_not_exist.h5")


def test_non_openpmd_hdf5():
    with pytest.raises(OpenPMDError):
        openpmd_hdf5.HDF5Reader(hdf5=os.path.join(rootdir, "blank.h5"))
