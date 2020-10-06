"""
The `~plasmapy.dispersion` subpackage contains functions associated with the 
dispersion solver and other relevant quantities
"""
# __all__ will be auto populated below
__all__ = []

from .dispersionfunction import plasma_dispersion_func, plasma_dispersion_func_deriv

# auto populate __all__
for obj_name in list(globals()):
    if not (obj_name.startswith("__") or obj_name.endswith("__")):
        __all__.append(obj_name)
__all__.sort()
