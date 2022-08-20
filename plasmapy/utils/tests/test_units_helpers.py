"""Tests of `plasmapy.utils._units_helpers`."""

import astropy.units as u
import pytest

from astropy.constants import c, m_e
from collections import namedtuple

from plasmapy.utils._units_helpers import _get_physical_type_dict


def test_get_physical_type_dict_specific_example():
    units = [u.m, u.m**-3, u.m * u.s]
    quantities = [5 * unit for unit in units]
    expected = {quantity.unit.physical_type: quantity for quantity in quantities}
    new_physical_type_dict = _get_physical_type_dict(quantities)
    assert new_physical_type_dict == expected


velocity = u.get_physical_type("velocity")
mass = u.get_physical_type("mass")
dimensionless = u.get_physical_type(1)

test_case = namedtuple("case", ["collection", "kwargs", "expected"])


test_cases = [
    test_case(
        collection=(c, m_e),
        kwargs={},
        expected={velocity: c, mass: m_e},
    ),
    test_case(
        collection=(c, 5),
        kwargs={"only_quantities": True, "numbers_become_quantities": True},
        expected={velocity: c, dimensionless: u.Quantity(5)},
    ),
    test_case(
        collection=(c, 5),
        kwargs={"only_quantities": False, "numbers_become_quantities": True},
        expected={velocity: c, dimensionless: u.Quantity(5)},
    ),
    test_case(
        collection=(c, 5),
        kwargs={"only_quantities": True, "numbers_become_quantities": False},
        expected={velocity: c},
    ),
    test_case(
        collection=(c, 5),
        kwargs={"only_quantities": False, "numbers_become_quantities": False},
        expected={velocity: c, dimensionless: 5},
    ),
    test_case(
        collection=("...", "....", set()),
        kwargs={"only_quantities": False, "numbers_become_quantities": False},
        expected={},
    ),
]


@pytest.mark.parametrize("collection, kwargs, expected", test_cases)
def test_get_physical_type_dict(collection, kwargs, expected):
    physical_type_dict = _get_physical_type_dict(collection, **kwargs)
    assert physical_type_dict == expected


test_cases_exceptions = [
    test_case(
        collection=(u.m, 5 * u.m),
        kwargs={},
        expected=ValueError,
    ),
]


@pytest.mark.parametrize("collection, kwargs, expected", test_cases_exceptions)
def test_get_physical_type_dict_exceptions(collection, kwargs, expected):
    with pytest.raises(expected):
        _get_physical_type_dict(collection, **kwargs)
