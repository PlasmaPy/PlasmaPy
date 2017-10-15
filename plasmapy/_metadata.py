##
# Package metadata
##

import ah_bootstrap

from astropy_helpers.git_helpers import get_git_devstr

# Name
name = 'plasmapy'

# PlasmaPy uses Semantic Versioning of the form: MAJOR.MINOR.PATCH
#
#  - The MAJOR version changes when there are backwards incompatible changes
#  - The MINOR version changes when backwards compatible functionality is added
#  - The PATCH version changes when the public API remains the same
#
# During initial development releases (with MAJOR = 0), backwards compatibility
# does not need to be maintained when MINOR is incremented.
#
# While a new version is being developed, '.dev' followed by the commit number
# will be appended to the version string.

version = '0.1.0.dev'

release = 'dev' in version

if release:
    version += get_git_devstr(False)

# Long description / docstring
description = """
PlasmaPy is a community-developed and community-driven core Python
package for plasma physics.
"""

# Author(s)
author = 'The PlasmaPy Community'
