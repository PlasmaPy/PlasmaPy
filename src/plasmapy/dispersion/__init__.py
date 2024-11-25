"""
The `~plasmapy.dispersion` subpackage contains functionality associated
with plasma dispersion relations, including numerical solvers and
`~plasmapy.dispersion.analytical` solutions. Plasa dispersion dispersion 
refers to the dependence of a wave’s phase velocity and group velocity on 
its frequency. This is a crucial aspect when studying wave propagation in 
plasmas, as it influences phenomena like wave-particle interactions, energy
transport, and the stability of plasma configurations.

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
