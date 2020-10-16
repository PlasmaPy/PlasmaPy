"""
<<<<<<< HEAD
The `~plasmapy.dispersion` subpackage contains functionality associated with
plasma dispersion relations, solvers and analytical solutions.
=======
The `~plasmapy.dispersion` subpackage contains functions associated with the of
dispersion solver and other relevant quantities
>>>>>>> dsolver
"""
__all__ = [
    "plasma_dispersion_func",
    "plasma_dispersion_func_deriv",
]

from plasmapy.dispersion.dispersionfunction import (
    plasma_dispersion_func,
    plasma_dispersion_func_deriv,
)
<<<<<<< HEAD
=======

from plasmapy.dispersion.two_fluid_dispersion_solver import (
    two_fluid_dispersion_solution,
)

# auto populate __all__
for obj_name in list(globals()):
    if not (obj_name.startswith("__") or obj_name.endswith("__")):
        __all__.append(obj_name)
__all__.sort()
>>>>>>> dsolver
