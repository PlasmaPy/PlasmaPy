:tocdepth: 3

.. _plasmapy-documentation:

######################
PlasmaPy Documentation
######################

.. raw:: html

    <div style="width: 30%;" class="sphx-glr-thumbcontainer" tooltip="Let&#x27;s try to look at ITER plasma conditions using the physics subpackage. ">

.. only:: html

    .. figure:: /auto_examples/images/thumb/sphx_glr_plot_physics_thumb.png

        :ref:`sphx_glr_auto_examples_plot_physics.py`

.. raw:: html

    </div>

.. raw:: html

    <div style="width: 30%;" class="sphx-glr-thumbcontainer" tooltip="Let&#x27;s import some basics (and PlasmaPy!) ">

.. only:: html

    .. figure:: /auto_examples/images/thumb/sphx_glr_plot_dispersion_function_thumb.png

        :ref:`sphx_glr_auto_examples_plot_dispersion_function.py`

.. raw:: html

    </div>

.. raw:: html

    <div style="width: 30%;" class="sphx-glr-thumbcontainer" tooltip="Let&#x27;s analyze a few Langmuir probe characteristics using the diagnostics.langmuir subpackage. F...">

.. only:: html

    .. figure:: /auto_examples/images/thumb/sphx_glr_plot_langmuir_analysis_thumb.png

        :ref:`sphx_glr_auto_examples_plot_langmuir_analysis.py`

.. raw:: html

    </div>

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
   :caption: Getting started
   :maxdepth: 1

   install
   COMMUNICATION
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

.. toctree::
   :maxdepth: 1
   :caption: Formulary

   formulary/index

The `~plasmapy.formulary` subpackage aims to cover the `NRL Plasma Physics
Formulary <https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary>`_.


.. toctree::
   :maxdepth: 1
   :caption: Experimental tools

   diagnostics/index

The `~plasmapy.diagnostics` package is in the early stages of
development.

.. toctree::
    :maxdepth: 1
    :caption: Data structures and simulation

    plasma/index
    simulation/particletracker

.. toctree::
    :maxdepth: 1
    :caption: Particle data

    particles/index

.. toctree::
   :maxdepth: 1
   :caption: Utilities

   utils/index

.. toctree::
    :maxdepth: 1
    :caption: Examples

    auto_examples/index


.. toctree::
    :maxdepth: 1
    :caption: Development guide

    development/index


The :ref:`plasmapy-development-guide` contains information on how to
contribute to PlasmaPy, along with guidelines for code, testing, and
documentation.


.. toctree::
   :maxdepth: 2
   :caption: About

   about/index

.. The about PlasmaPy section has some important information that would
   be helpful to have more readily accessible from the main doc index
   page.

*****
Index
*****

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. TODO: Add feedback link: .. _feedback@plasmapy.org: mailto:feedback@plasmapy.org
