"""
Module for loading atomic data for elements from
:file:`plasmapy/particles/data/elements.json`.

The periodic tabla data is from: http://periodic.lanl.gov/index.shtml

.. attention::
    This module is not part of PlasmaPy's public API.
"""
__all__ = []

import astropy.units as u
import collections
import json
import pkgutil

_PeriodicTable = collections.namedtuple(
    "periodic_table", ["group", "category", "block", "period"]
)


def _element_obj_hook(obj):
    if "unit" in obj:
        return obj["value"] * u.Unit(obj["unit"])
    return obj


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


_elements = json.loads(
    pkgutil.get_data("plasmapy", "particles/data/elements.json"),
    object_hook=_element_obj_hook,
)


_atomic_numbers_to_symbols = {
    elemdict["atomic number"]: symb for (symb, elemdict) in _elements.items()
}

_element_names_to_symbols = {
    elemdict["element name"]: symb for (symb, elemdict) in _elements.items()
}
