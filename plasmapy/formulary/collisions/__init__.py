"""
The `~plasmapy.formulary.collisions` subpackage contains commonly
used collisional formulae from plasma science.
"""
__all__ = []

import inspect

from plasmapy.formulary.collisions.coulomb import *
from plasmapy.formulary.collisions.dimensionless import *
from plasmapy.formulary.collisions.frequencies import *
from plasmapy.formulary.collisions.lengths import *
from plasmapy.formulary.collisions.misc import *

# auto populate __all__
for obj_name in list(globals()):
    if not (
        obj_name.startswith("__") or obj_name.endswith("__")
    ) and not inspect.ismodule(globals()[obj_name]):
        __all__.append(obj_name)
__all__.sort()

del inspect
