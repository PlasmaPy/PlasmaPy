
.. _notebook guide:

**************
Notebook Guide
**************

PlasmaPy's online documentation includes :ref:`example notebooks
<examples>`. Jupyter notebooks may be added to PlasmaPy's
`example gallery <example-gallery>`_ by putting them in
:file:`docs/notebooks` or one of its subdirectories. The directory or
notebook must be included in :file:`docs/examples.rst` to show up in the
gallery.

These notebooks are and built using the |nbsphinx|_ extension to Sphinx_
and thus show up in the online documentation. The `Jupyter Project
Documentation`_ has more information on how to use Jupyter notebooks.

.. _Jupyter Project Documentation: https://docs.jupyter.org/en/latest/

.. _example Markdown cells: https://nbsphinx.readthedocs.io/en/latest/markdown-cells.html

Need to discuss `example Markdown cells`_...

.. code-block:: markdown

   [example notebooks]: https://docs.plasmapy.org/en/latest/examples.html
   [astropy.units]: https://docs.astropy.org/en/stable/units/index.html

   PlasmaPy's documentation includes [example notebooks]. We can refer
   to code objects by using backticks, such as for [astropy.units].


Markdown cells
==============


Links
-----

Note that links must be redefined in each Markdown_ cell.

External links
~~~~~~~~~~~~~~

.. code-block:: markdown

   [cosmic latte]: https://en.wikipedia.org/wiki/Cosmic_latte

   The average color of the universe is [cosmic latte].

Internal links
~~~~~~~~~~~~~~

Linking to an object in PlasmaPy's API requires including a relative
link to the corresponding file that gets created in :file:`docs/api`
when the documentation is built.

If a notebook is included in :file:`docs/notebooks/getting_started`,
then a link to |Particle| can be created as

.. code-block:: markdown

   [Particle]: ../../api/plasmapy.particles.particle_class.Particle.rst

   This is a link to the [Particle] class.

The

Glossary terms
~~~~~~~~~~~~~~

.. code-block:: markdown

   [particle-like]: ../../glossary.rst#term-particle-like

   An object is [particle-like] if it can be transformed into a
   `Particle`.

Table of contents
~~~~~~~~~~~~~~~~~

.. code-block:: markdown

   ## Contents

   1. [Introduction](#Introduction)
   2. [Second section](#Second section)

   ## Introduction

   ## Second section

Intentional exceptions
----------------------

If a cell is expected to raise an exception, label it with

reStructuredText cells
======================


Pre-executing notebooks
=======================

Some of the notebooks in the example gallery are

Some of the notebooks needed in
