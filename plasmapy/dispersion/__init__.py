"""
The `~plasmapy.dispersion` subpackage contains functionality associated with
plasma dispersion relations, solvers and analytical solutions.
"""
__all__ = [
    "plasma_dispersion_func",
    "plasma_dispersion_func_deriv",
]

from plasmapy.dispersion.dispersionfunction import (
    plasma_dispersion_func,
    plasma_dispersion_func_deriv,
)
