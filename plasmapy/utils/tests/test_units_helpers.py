import astropy.units as u
import astropy.units.physical as physical
import pytest

from plasmapy.utils.units_helpers import get_physical_type_dict


@pytest.fixture
def pre_physical_type_dict():
    units = (u.m, u.m ** -3, u.m * u.s)
    return {unit.physical_type: 5 * unit for unit in units}


def test_get_physical_type(pre_physical_type_dict):
    units = pre_physical_type_dict.values()
    new_physical_type_dict = get_physical_type_dict(units)
    assert pre_physical_type_dict == new_physical_type_dict
