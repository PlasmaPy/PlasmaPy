"""
The `~plasmapy.dispersion.numerical` subpackage contains functionality
associated with numerical dispersion solvers.
"""
__all__ = ["hollweg","stix"]

from plasmapy.dispersion.numerical.hollweg_ import hollweg
from plasmapy.dispersion.numerical.stix_ import stix
