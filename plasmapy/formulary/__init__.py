"""
The `~plasmapy.formulary` subpackage contains commonly used formulae
from plasma science.
"""
# __all__ & __aliases__ will be auto populated below
__all__ = []
__aliases__ = []

from plasmapy.formulary.braginskii import *
from plasmapy.formulary.collisions import *
from plasmapy.formulary.dielectric import *
from plasmapy.formulary.dimensionless import *
from plasmapy.formulary.distribution import *
from plasmapy.formulary.drifts import *
from plasmapy.formulary.ionization import *
from plasmapy.formulary.magnetostatics import *
from plasmapy.formulary.mathematics import *
from plasmapy.formulary.parameters import *
from plasmapy.formulary.quantum import *
from plasmapy.formulary.relativity import *

# auto populate __all__
for obj_name in list(globals()):
    if not (obj_name.startswith("__") or obj_name.endswith("__")):
        __all__.append(obj_name)
__all__.sort()

# auto populate __aliases__
for modname in (
    "braginskii",
    "collisions",
    "dielectric",
    "dimensionless",
    "distribution",
    "drifts",
    "ionization",
    "magnetostatics",
    "mathematics",
    "parameters",
    "quantum",
    "relativity",
):
    try:
        obj = globals()[modname]
        __aliases__.extend(obj.__aliases__)
    except (KeyError, AttributeError):
        pass
__aliases__.sort()

# cleanup namespace
del modname, obj, obj_name
