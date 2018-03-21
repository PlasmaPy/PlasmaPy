.. PlasmaPy documentation master file, created by
   sphinx-quickstart on Wed May 31 18:16:46 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

######################
PlasmaPy Documentation
######################

`PlasmaPy <http://www.plasmapy.org/`_ is an open source
community-developed core `Python <https://www.python.org/>`_ 3.6+
package for plasma physics in the early stages of development. The
documentation contained on this page is in the process of being written.

***************
Getting Started
***************

* `Installing PlasmaPy
  <https://github.com/PlasmaPy/PlasmaPy/blob/master/INSTALL.md>`_
* `Contributing to PlasmaPy
  <https://github.com/PlasmaPy/PlasmaPy/blob/master/CONTRIBUTING.md>`_
* `Code of conduct
  <https://github.com/PlasmaPy/PlasmaPy/blob/master/CODE_OF_CONDUCT.md>`_
* `PlasmaPy on GitHub
  <https://github.com/PlasmaPy/plasmapy>`_
* `PlasmaPy website
  <http://www.plasmapy.org/>`_

.. user-documentation

******************
User Documentation
******************

.. particles

Particles
---------

.. toctree::
   :maxdepth: 2

   atomic/particle_class
   atomic/functional
   atomic/nuclear
   atomic/decorators

.. plasma-parameters

Plasma Parameters
-----------------

.. toctree::
   :maxdepth: 3

   physics/index
   physics/transport/index

.. _diagnostics

Diagnostics
-----------

The `~plasmapy.diagnostics` package is in the early stages of
development.

.. toctree::
   :maxdepth: 1

   diagnostics/index

.. _mathematics

Mathematics
-----------

The `plasmapy.mathematics` package contains mathematical functions
commonly used by plasma physicists.

.. toctree::
   :maxdepth: 3

   mathematics/index

.. _data-structures

Data Structures
---------------

.. toctree::
    :maxdepth: 1

    plasma/index
    species/index

.. _utilities

Utilities
---------

.. toctree::
   :maxdepth: 3

   utils/index

.. _examples

Examples
--------
.. toctree::
    :maxdepth: 3

    auto_examples/index

.. _development-guide

*****************
Development Guide
*****************

The development guide contains information on how to contribute to
PlasmaPy, along with guidelines for code, testing, and documentation.

.. toctree::
    :maxdepth: 2

    development/code_guide
    development/testing_guide
    development/doc_guide
    development/release_guide

.. _project_details:

***************
Project Details
***************

.. toctree::
   :maxdepth: 2

* `Vision Statement
  <https://github.com/PlasmaPy/PlasmaPy/blob/master/VISION.md>`_

.. _index:

*****
Index
*****

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. Did not include one for `search` since that page didn't show up
   satisfactorily

.. TODO: Add feedback link: .. _feedback@plasmapy.org: mailto:feedback@plasmapy.org
