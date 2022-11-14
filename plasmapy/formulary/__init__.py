"""
The `~plasmapy.formulary` subpackage contains commonly used formulae
from plasma science.
"""
__all__ = []
__aliases__ = []
__lite_funcs__ = []

from plasmapy.formulary.braginskii import *
from plasmapy.formulary.collisions import *
from plasmapy.formulary.densities import *
from plasmapy.formulary.dielectric import *
from plasmapy.formulary.dimensionless import *
from plasmapy.formulary.distribution import *
from plasmapy.formulary.drifts import *
from plasmapy.formulary.frequencies import *
from plasmapy.formulary.ionization import *
from plasmapy.formulary.lengths import *
from plasmapy.formulary.magnetostatics import *
from plasmapy.formulary.mathematics import *
from plasmapy.formulary.misc import *
from plasmapy.formulary.quantum import *
from plasmapy.formulary.radiation import *
from plasmapy.formulary.relativity import *
from plasmapy.formulary.speeds import *

# auto populate __all__
import inspect  # isort: skip

for obj_name in list(globals()):
    if not (
        obj_name.startswith("__") or obj_name.endswith("__")
    ) and not inspect.ismodule(globals()[obj_name]):
        __all__.append(obj_name)
__all__.sort()

# Put non-formulary imports here so that they don't get included in __all__

import contextlib  # isort: skip

# auto populate __aliases__ & __lite_funcs__
for modname in (
    "braginskii",
    "collisions",
    "densities",
    "dielectric",
    "dimensionless",
    "distribution",
    "drifts",
    "frequencies",
    "ionization",
    "lengths",
    "magnetostatics",
    "mathematics",
    "misc",
    "parameters",
    "quantum",
    "radiation",
    "relativity",
    "speeds",
):
    try:
        obj = globals()[modname]
    except KeyError:  # coverage: ignore
        continue

    with contextlib.suppress(AttributeError):
        __aliases__.extend(obj.__aliases__)

    with contextlib.suppress(AttributeError):
        __lite_funcs__.extend(obj.__lite_funcs__)

__aliases__ = list(sorted(set(__aliases__)))
__lite_funcs__ = list(sorted(set(__lite_funcs__)))

# cleanup namespace
del contextlib, inspect, modname, obj, obj_name
