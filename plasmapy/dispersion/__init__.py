"""
The `~plasmapy.dispersion` subpackage contains functionality associated
with plasma dispersion relations, including numerical solvers and
`~plasmapy.dispersion.analytical` solutions.

.. attention::

   |expect-api-changes|
"""
__all__ = [
    "analytical",
    "dispersion_functions",
    "numerical",
    "plasma_dispersion_func",
    "plasma_dispersion_func_deriv",
]

from plasmapy.dispersion import analytical, dispersion_functions, numerical
from plasmapy.dispersion.dispersion_functions import (
    plasma_dispersion_func,
    plasma_dispersion_func_deriv,
)
