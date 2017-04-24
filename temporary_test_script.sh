"""A temporary script to run PlasmaPy tests.  This will become
obsolete and should be removed when we set up continuous integration
testing.  From the command line, run:

    python test_script.sh

Here, 'python' may need to be replaced by 'python3'.  Remember:
PlasmaPy will guarantee support for Python 3.6 and above.

"""

import plasmapy as pl

# analytic subpackage tests

pl.analytic.tests.test_plasma_dispersion_func()
pl.analytic.tests.test_plasma_dispersion_func_deriv()
