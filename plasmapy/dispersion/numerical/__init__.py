"""
The `~plasmapy.dispersion.numerical` subpackage contains functionality
associated with numerical dispersion solvers.
"""
__all__ = ["hollweg", "kinetic_alfven_"]

from plasmapy.dispersion.numerical.hollweg_ import hollweg
from plasmapy.dispersion.numerical.kinetic_alfven_ import kinetic_alfven
