"""
The `~plasmapy.dispersion.numerical` subpackage contains functionality
associated with numerical dispersion solvers.

.. attention::

   |expect-api-changes|
"""
__all__ = [
    "hirose",
    "hollweg",
]

from plasmapy.dispersion.numerical.hirose_ import hirose
from plasmapy.dispersion.numerical.hollweg_ import hollweg
