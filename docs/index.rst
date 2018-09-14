:tocdepth: 3

.. _plasmapy-documentation:

######################
PlasmaPy Documentation
######################

.. image:: _static/graphic-circular.png
   :alt: PlasmaPy logo
   :align: right
   :scale: 80%

`PlasmaPy <http://www.plasmapy.org/>`_ is an open source
community-developed core `Python <https://www.python.org/>`_ 3.6+
package for plasma physics currently under development.

.. _toplevel-getting-started:

***************
Getting Started
***************

.. toctree::
   :maxdepth: 1

   install
   CONTRIBUTING
   CODE_OF_CONDUCT
   about/citation

* `PlasmaPy's GitHub repository
  <https://github.com/PlasmaPy/plasmapy>`_
* `PlasmaPy website
  <http://www.plasmapy.org/>`_
* `Using astropy.units <http://docs.astropy.org/en/stable/units/>`_

.. _toplevel-user-documentation:

******************
User Documentation
******************

.. _toplevel-plasma-parameters:

Theoretical Analysis
--------------------

.. toctree::
   :maxdepth: 1

   physics/index
   transport/index
   mathematics/index

.. _toplevel-experimental-tools:

Experimental Tools
------------------

The `~plasmapy.diagnostics` package is in the early stages of
development.

.. toctree::
   :maxdepth: 1

   diagnostics/index

.. _toplevel-data-structures:

Data Structures and Simulation
------------------------------

.. toctree::
    :maxdepth: 1

    plasma/index
    species/index

.. _toplevel-physical-data:

Physical Data
-------------

.. toctree::
    :maxdepth: 1

    constants/index
    atomic/index

.. _toplevel-utilities:

Utilities
---------

.. toctree::
   :maxdepth: 1

   utils/index

.. _toplevel-examples:

Examples
--------
.. toctree::
    :maxdepth: 1

    auto_examples/index

.. _toplevel-development-guide:

*****************
Development Guide
*****************

The :ref:`plasmapy-development-guide` contains information on how to
contribute to PlasmaPy, along with guidelines for code, testing, and
documentation.

.. toctree::
    :maxdepth: 1
    :hidden:

    development/index

.. _toplevel-project-details:

***************
Project Details
***************

.. The about PlasmaPy section has some important information that would
   be helpful to have more readily accessible from the main doc index
   page.

.. toctree::
   :maxdepth: 2

   about/index

.. _toplevel-index:

*****
Index
*****

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. TODO: Add feedback link: .. _feedback@plasmapy.org: mailto:feedback@plasmapy.org
