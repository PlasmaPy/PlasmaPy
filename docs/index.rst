:tocdepth: 3

.. _plasmapy-documentation:

.. image:: _static/graphic-circular.png
   :alt: PlasmaPy logo
   :align: right
   :scale: 40%

######################
PlasmaPy Documentation
######################

PlasmaPy_ is an open source community-developed core Python_
|minpython|\ + package for plasma physics currently under development.

Example highlights
------------------

.. nbgallery::
   :hidden:

   notebooks/getting_started/particles
   notebooks/diagnostics/charged_particle_radiography_particle_tracing
   notebooks/dispersion/two_fluid_dispersion
   notebooks/diagnostics/thomson
   notebooks/analysis/swept_langmuir/find_floating_potential
   notebooks/formulary/thermal_bremsstrahlung

.. toctree::
   :caption: First Steps
   :maxdepth: 1

   Installing <install>
   getting_started
   examples
   COMMUNICATION
   CONTRIBUTING
   Code of Conduct <CODE_OF_CONDUCT.rst>
   about/citation

.. toctree::
   :maxdepth: 1
   :caption: Package features

   ad/index
   Dispersion <dispersion/index>
   Formulary <formulary/index>
   Particles <particles/index>
   Simulation <simulation/index>
   Plasma Calculator <plasma_calculator/index>
   Plasma Objects <plasma/index>
   Package Utilities <utils/index>

.. toctree::
   :maxdepth: 1
   :caption: Contributor Guide

   Overview <contributing/index>
   contributing/code_guide
   contributing/changelog_guide
   contributing/doc_guide
   contributing/testing_guide
   contributing/release_guide

.. toctree::
   :maxdepth: 1
   :caption: All the Rest

   whatsnew/index
   about/credits
   bibliography
   glossary
   Vision Statement <about/vision_statement>
   PlasmaPy.org <https://www.plasmapy.org>

.. The about PlasmaPy section has some important information that would
   be helpful to have more readily accessible from the main doc index
   page.

.. TODO: Add feedback link: .. _feedback@plasmapy.org: mailto:feedback@plasmapy.org
