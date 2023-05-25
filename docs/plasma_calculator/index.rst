.. _plasmapy-calculator:

=================
Plasma Calculator
=================

.. currentmodule:: plasmapy.utils.calculator

Overview
--------

The Plasma Calculator is an interactive Jupyter notebook that is packaged
with `plasmapy`, and allows users to input a set of plasma properties and
immediately calculate multiple plasma parameters.

.. note::

   This functionality is still under development and the API may change
   in future releases.

Using Plasma Calculator
-----------------------

To invoke the app use ``plasma-calculator`` in the command line.
By default this opens the app in a browser, bootstrapping
the notebook in a light theme and using port 8866 by default.

``plasma-calculator`` takes optional arguments such as ``--dark``,
``--port``, and ``--no-browser``. Pass flag ``-h`` or ``--help`` to get
full list of supported arguments.
