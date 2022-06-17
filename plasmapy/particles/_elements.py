"""
Module for loading atomic data for elements from
:file:`plasmapy/particles/data/elements.json`.

The periodic tabla data is from: https://periodic.lanl.gov/index.shtml
"""
__all__ = [
    "element_obj_hook",
    "data_about_elements",
    "atomic_numbers_to_symbols",
    "element_names_to_symbols",
]

import astropy.units as u
import json
import pkgutil

from dataclasses import dataclass


@dataclass
class PeriodicTable:
    """A data container for the periodic table information for an element."""

    group: int
    period: int
    block: str
    category: str


def element_obj_hook(obj):
    """Provide an ``object_hook`` designed for `json.load` and `json.loads`."""
    return obj["value"] * u.Unit(obj["unit"]) if "unit" in obj else obj


# this code was used to create the JSON file as per vn-ki on Matrix:
# https://matrix.to/#/!hkWCiyhQyxiYJlUtKF:matrix.org/
#    $1554667515670438wIKlP:matrix.org?via=matrix.org&via=cadair.com
#
# def plasma_default(obj):
#     if isinstance(obj, u.Quantity):
#         return {
#             "unit": obj.unit.name,
#             "value": obj.value,
#         }
#
# with open("elements.json", "w") as f:
#    json.dump(_Elements, f, default=plasma_default, indent=2)


data_about_elements = json.loads(
    pkgutil.get_data("plasmapy", "particles/data/elements.json"),
    object_hook=element_obj_hook,
)


atomic_numbers_to_symbols = {
    elemdict["atomic number"]: symb for (symb, elemdict) in data_about_elements.items()
}

element_names_to_symbols = {
    elemdict["element name"]: symb for (symb, elemdict) in data_about_elements.items()
}
