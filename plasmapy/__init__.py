from .utils import import_helpers

__name__ = 'plasmapy'

__doc__ = ('A community-developed and community-driven open source core '
           'Python package for plasma physics.')

try:
    from .version import version as __version__
except ImportError:
    __version__ = '0.1.0dev0'

__minimum_python_version__ = '3.6'

# The following information is duplicated in requirements/base.txt

__minimum_versions__ = {
    'numpy': '1.13',
    'astropy': '2.0',
    'scipy': '0.19',
    }

import_helpers.check_python(__minimum_python_version__)
import_helpers.check_versions(__minimum_versions__)

try:
    from .classes import Plasma
    from . import classes
    from . import constants
    from . import atomic
    from . import mathematics
    from . import physics
    from . import utils
except ImportError as exc:
    raise ImportError('Unable to load PlasmaPy subpackages.') from exc

# A more extensive and thoughtful method for cleaning up our top-level
# namespace is in Astropy's __init__.py (see also pull request #210).

del (import_helpers, __minimum_python_version__, __minimum_versions__)
