from .utils import import_helpers

__name__ = 'plasmapy'

__doc__ = ('A community-developed and community-driven open source core '
           'Python package for plasma physics.')

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''

__minimum_python_version__ = '3.6'

__minimum_versions__ = {
    'numpy': '1.13',
    'astropy': '2.0',
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

# Clean up the top-level namespace by deleting everything that (1)
# isn't listed in keep, (2) is a dunder, and (3) isn't a submodule of
# this package.  This was adapted from astropy's __init__.py.

__keep__ = ['Plasma']

from types import ModuleType as __module_type__

for varname in dir():

    in_keep = varname in __keep__

    is_dunder = varname.startswith('__') and varname.endswith('__')

    # When using relative imports like ``from .. import config``, the
    # ``config`` variable is automatically created in the namespace of
    # whatever module ``..`` resolves to (in this case astropy).  This
    # happens a few times just in the module setup above.  This allows
    # the cleanup to keep any public submodules of the package.

    is_subpackage = (varname[0] != '_' and
                     isinstance(locals()[varname], __module_type__) and
                     locals()[varname].__name__.startswith(__name__ + '.'))

    if not (in_keep or is_dunder or is_subpackage):
        del locals()[varname]

del (import_helpers, __minimum_python_version__, __minimum_versions__,
     varname, __module_type__, __keep__, in_keep, is_dunder, is_subpackage)
