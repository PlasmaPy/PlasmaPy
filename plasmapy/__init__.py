
"""
PlasmaPy is a community-developed and community-driven core Python
package for plasma physics.
"""

import sys
from warnings import warn

if sys.version_info[:2] < (3, 6):
    warn("PlasmaPy does not support Python 3.5 and below")

from . import analytic_functions
