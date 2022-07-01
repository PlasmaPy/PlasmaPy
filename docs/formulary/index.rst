.. _formulary:

********************************
Formulary (`plasmapy.formulary`)
********************************

`plasmapy.formulary` provides theoretical formulas for calculation of
physical quantities helpful for plasma physics.

.. table::
   :widths: 5 16

   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Classical Transport <braginskii>          | `plasmapy.formulary.braginskii`         |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Collisions <collisions>                   | `plasmapy.formulary.collisions`         |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Dielectrics <dielectric>                  | `plasmapy.formulary.dielectric`         |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Dimensionless <dimensionless>             | `plasmapy.formulary.dimensionless`      |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Distribution Functions <distribution>     | `plasmapy.formulary.distribution`       |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Drifts <drifts>                           | `plasmapy.formulary.drifts`             |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Frequencies <frequencies>                 | `plasmapy.formulary.frequencies`        |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Ionization <ionization>                   | `plasmapy.formulary.ionization`         |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Lengths <lengths>                         | `plasmapy.formulary.lengths`            |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Magnetostatics <magnetostatics>           | `plasmapy.formulary.magnetostatics`     |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Mathematics <mathematics>                 | `plasmapy.formulary.mathematics`        |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Miscellaneous Parameters <misc>           | `plasmapy.formulary.misc`               |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Plasma Parameters <parameters>            | `plasmapy.formulary.parameters`         |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Quantum Relations <quantum>               | `plasmapy.formulary.quantum`            |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Radiation <radiation>                     | `plasmapy.formulary.radiation`          |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Relativistic Relations <relativity>       | `plasmapy.formulary.relativity`         |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+
   | .. toctree:: Speeds <speeds>                           | `plasmapy.formulary.speeds`             |
   |    :maxdepth: 1                                        |                                         |
   +--------------------------------------------------------+-----------------------------------------+

The subpackage makes heavy use of |Quantity| for handling conversions
between different unit systems. This is especially important for
electron-volts, commonly used in plasma physics to denote temperature,
although it is technically a unit of energy.

Most functions expect |Quantity| objects as inputs, however some will
use the `~plasmapy.utils.decorators.validators.validate_quantities`
decorator to automatically cast arguments to |Quantity| objects with the
appropriate units. If that happens, you will be notified via a
`astropy.units.UnitsWarning`.

Please note that well-maintained physical constant data with units and
uncertainties can be found in `astropy.constants`.

Examples
========

For a general overview of how unit-based input works, take a look at the
following examples:

.. nbgallery::
    :caption: Examples
    :glob:

    ../notebooks/formulary/*

API
===

.. automodapi:: plasmapy.formulary
