import pytest

from astropy import units as u
from typing import List, Tuple, Union

import plasmapy.plasma

from plasmapy.particles.data.test import data_dir
from plasmapy.plasma.exceptions import DataStandardError
from plasmapy.plasma.sources import openpmd_hdf5


@pytest.fixture(scope="module")
def h5_2d(request):
    h5 = openpmd_hdf5.HDF5Reader(hdf5=data_dir / "data00000255.h5")
    yield h5
    h5.close()


@pytest.fixture(scope="module")
def h5_3d(request):
    h5 = openpmd_hdf5.HDF5Reader(hdf5=data_dir / "data00000100.h5")
    yield h5
    h5.close()


@pytest.fixture(scope="module")
def h5_theta(request):
    h5 = openpmd_hdf5.HDF5Reader(hdf5=data_dir / "data00000200.h5")
    yield h5
    h5.close()


@pytest.mark.slow
class TestOpenPMD2D:
    """Test 2D HDF5 dataset based on OpenPMD."""

    # Downloaded from
    # https://github.com/openPMD/openPMD-example-datasets/blob/draft/example-2d.tar.gz
    # per the Creative Commons Zero v1.0 Universal license

    def test_has_electric_field_with_units(self, h5_2d):
        assert isinstance(h5_2d.electric_field.to(u.V / u.m), u.Quantity)

    def test_correct_shape_electric_field(self, h5_2d):
        assert h5_2d.electric_field.shape == (3, 51, 201)

    def test_has_charge_density_with_units(self, h5_2d):
        # this should simply pass without exception
        h5_2d.charge_density.to(u.C / u.m**3)

    def test_correct_shape_charge_density(self, h5_2d):
        assert h5_2d.charge_density.shape == (51, 201)

    def test_has_magnetic_field(self, h5_2d):
        with pytest.raises(AttributeError):
            h5_2d.magnetic_field

    def test_has_electric_current(self, h5_2d):
        with pytest.raises(AttributeError):
            h5_2d.electric_current


class TestOpenPMD3D:
    """Test 3D HDF5 dataset based on OpenPMD."""

    # Downloaded from
    # https://github.com/openPMD/openPMD-example-datasets/blob/draft/example-3d.tar.gz
    # per the Creative Commons Zero v1.0 Universal license

    def test_has_electric_field_with_units(self, h5_3d):
        assert isinstance(h5_3d.electric_field.to(u.V / u.m), u.Quantity)

    def test_correct_shape_electric_field(self, h5_3d):
        assert h5_3d.electric_field.shape == (3, 26, 26, 201)

    def test_has_charge_density_with_units(self, h5_3d):
        assert isinstance(h5_3d.charge_density.to(u.C / u.m**3), u.Quantity)

    def test_correct_shape_charge_density(self, h5_3d):
        assert h5_3d.charge_density.shape == (26, 26, 201)

    def test_has_magnetic_field(self, h5_3d):
        with pytest.raises(AttributeError):
            h5_3d.magnetic_field

    def test_has_electric_current(self, h5_3d):
        with pytest.raises(AttributeError):
            h5_3d.electric_current


@pytest.mark.slow
class TestOpenPMDThetaMode:
    """Test thetaMode HDF5 dataset based on OpenPMD."""

    # Downloaded from
    # https://github.com/openPMD/openPMD-example-datasets/blob/draft/example-thetaMode.tar.gz
    # per the Creative Commons Zero v1.0 Universal license

    def test_has_electric_field_with_units(self, h5_theta):
        assert isinstance(h5_theta.electric_field.to(u.V / u.m), u.Quantity)

    def test_correct_shape_electric_field(self, h5_theta):
        assert h5_theta.electric_field.shape == (3, 3, 51, 201)

    def test_has_charge_density_with_units(self, h5_theta):
        assert isinstance(h5_theta.charge_density.to(u.C / u.m**3), u.Quantity)

    def test_correct_shape_charge_density(self, h5_theta):
        assert h5_theta.charge_density.shape == (3, 51, 201)

    def test_has_magnetic_field_with_units(self, h5_theta):
        assert isinstance(h5_theta.magnetic_field.to(u.T), u.Quantity)

    def test_correct_shape_magnetic_field(self, h5_theta):
        assert h5_theta.magnetic_field.shape == (3, 3, 51, 201)

    def test_has_electric_current_with_units(self, h5_theta):
        assert isinstance(
            h5_theta.electric_current.to(u.A * u.kg / u.m**3), u.Quantity
        )

    def test_correct_shape_electric_current(self, h5_theta):
        assert h5_theta.electric_current.shape == (3, 3, 51, 201)


units_test_table = [
    ((1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0), u.m * u.kg / u.amp * u.cd**2),
    ((1, 0, 1, 2, 0, 0, 0), u.m * u.s * u.amp**2),
    ([-3.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0], u.coulomb / u.m**3),
    ([2, 1, -3, -2, 0, 0, 0], u.ohm),
]


@pytest.mark.parametrize("openPMD_dims, expected", units_test_table)
@pytest.mark.slow
def test_fetch_units(openPMD_dims, expected: Union[Tuple, List]):
    units = openpmd_hdf5._fetch_units(openPMD_dims)
    assert units == expected


def test_unavailable_hdf5():
    with pytest.raises(FileNotFoundError):
        openpmd_hdf5.HDF5Reader(hdf5="this_file_does_not_exist.h5")


def test_non_openpmd_hdf5():
    with pytest.raises(DataStandardError):
        openpmd_hdf5.HDF5Reader(hdf5=data_dir.joinpath("blank.h5"))


def test_HDF5Reader(h5_2d):
    assert isinstance(h5_2d, plasmapy.plasma.sources.HDF5Reader)
    assert isinstance(h5_2d, plasmapy.plasma.BasePlasma)
