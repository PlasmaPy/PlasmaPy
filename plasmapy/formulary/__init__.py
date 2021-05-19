"""
The `~plasmapy.formulary` subpackage contains commonly used formulae
from plasma science.
"""
# __all__, __aliases__, & __lite_funcs__ will be auto populated below
__all__ = []
__aliases__ = []
__lite_funcs__ = []

from .braginskii import *
from .collisions import *
from .dielectric import *
from .dimensionless import *
from .distribution import *
from .drifts import *
from .ionization import *
from .magnetostatics import *
from .mathematics import *
from .parameters import *
from .quantum import *
from .relativity import *

# auto populate __all__
for obj_name in list(globals()):
    if not (obj_name.startswith("__") or obj_name.endswith("__")):
        __all__.append(obj_name)
__all__.sort()

# auto populate __aliases__ & __lite_funcs__
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
    except KeyError:
        continue

    try:
        __aliases__.extend(obj.__aliases__)
    except AttributeError:
        pass

    try:
        __lite_funcs__.extend(obj.__lite_funcs__)
    except AttributeError:
        pass

# cleanup namespace
del modname, obj, obj_name
