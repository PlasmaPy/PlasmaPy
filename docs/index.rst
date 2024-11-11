:tocdepth: 3

.. _plasmapy-documentation:

.. image:: _static/graphic-circular.png
   :alt: PlasmaPy logo
   :align: right
   :scale: 40%

######################
PlasmaPy Documentation
######################

|PlasmaPy| is an open source community-developed |Python| |minpython|\ +
package for plasma research and education. PlasmaPy is a platform by
which the plasma community can share code and collaboratively develop
new software tools for plasma research.

If you are new to PlasmaPy, please check out our
:ref:`getting started notebooks <getting-started-notebooks>` and our
:ref:`example gallery <examples>`. We invite you to share ideas and ask
questions in our |Matrix chat room| or during our weekly virtual
|office hours|.

PlasmaPy is developed openly `on GitHub`_, where you can
`request a new feature`_ or `report a bug`_.

.. important::

   If you use PlasmaPy for work presented in a publication or talk,
   please help the project by following these instructions to
   :ref:`cite or acknowledge <citation>` PlasmaPy.

.. toctree::
   :caption: First steps
   :maxdepth: 1

   Installing <install>
   getting_started
   examples
   COMMUNICATION
   Code of Conduct <CODE_OF_CONDUCT.rst>
   about/citation

Example notebooks
-----------------

.. nbgallery::
   :hidden:

   notebooks/getting_started/particles
   notebooks/diagnostics/charged_particle_radiography_particle_tracing
   notebooks/dispersion/two_fluid_dispersion
   notebooks/diagnostics/thomson
   notebooks/analysis/swept_langmuir/find_floating_potential
   notebooks/formulary/thermal_bremsstrahlung

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
   :caption: Project details

   changelog/index
   about/credits
   bibliography
   glossary
   performance_tips
   PlasmaPy Enhancement Proposals <https://github.com/PlasmaPy/PlasmaPy-PLEPs>
   PlasmaPy website <https://www.plasmapy.org>
   GitHub Repository <https://github.com/PlasmaPy/PlasmaPy>

.. toctree::
   :maxdepth: 1
   :caption: Contributing

   Contributor Guide <contributing/index>
   contributing/many_ways
   contributing/getting_ready
   contributing/workflow
   contributing/coding_guide
   contributing/testing_guide
   contributing/doc_guide
   contributing/changelog_guide
   contributing/pre-commit

.. _new discussion on GitHub: https://github.com/PlasmaPy/PlasmaPy/discussions/new/choose
.. _on GitHub: https://github.com/PlasmaPy/PlasmaPy
.. _report a bug: https://github.com/PlasmaPy/PlasmaPy/issues/new?assignees=&labels=Bug&template=bug_report.yml
.. _request a new feature: https://github.com/PlasmaPy/PlasmaPy/issues/new?assignees=&labels=Feature+request&template=feature_request.yml
