"""
Create a dictionary containing basic information for isotopes and
neutrons.
"""

import json
import pkgutil

import astropy.units as u

# this code was used to create the JSON file as per vn-ki on Riot:
# https://matrix.to/#/!hkWCiyhQyxiYJlUtKF:matrix.org/
#    $1554667515670438wIKlP:matrix.org?via=matrix.org&via=cadair.com
#
# def _isotope_default(obj):
#     if isinstance(obj, u.Quantity):
#         return {
#             "unit": obj.unit.name,
#             "value": obj.value,
#         }
# with open("isotopes.json", "w") as f:
#     json.dump(_Isotopes, f, default=plasma_default, indent=2)


def _isotope_obj_hook(obj):
    if "unit" in obj:
        return obj["value"] * u.Unit(obj["unit"])
    return obj


_Isotopes = json.loads(
    pkgutil.get_data("plasmapy", "particles/data/isotopes.json"),
    object_hook=_isotope_obj_hook,
)
