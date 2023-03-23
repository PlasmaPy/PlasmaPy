"""
The `~plasmapy.formulary.collisions.helio` subpackage contains
commonly used heliosphere plasma science, primarily the solar wind.
"""
__all__ = []

import inspect

from plasmapy.formulary.collisions.helio.collisional_analysis import *

# auto populate __all__
for obj_name in list(globals()):
    if not (
        obj_name.startswith("__") or obj_name.endswith("__")
    ) and not inspect.ismodule(globals()[obj_name]):
        __all__.append(obj_name)
__all__.sort()

del inspect