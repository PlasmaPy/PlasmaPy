.. _plasmapy-calculator:

=================
Plasma Calculator
=================

.. currentmodule:: plasmapy.utils.calculator

Overview
--------

Plasma calculator is an interactive notebook installed alongside
PlasmaPy that allows users to input a set of plasma properties and get
immediate output for multiple functions at once.

.. note::

   This functionality is still under development and the API may change
   in future releases.

Using Plasma Calculator
-----------------------

To invoke the app use ``plasma-calculator`` in the command line.
By default this opens the app in a browser, bootstrapping
the notebook in a light theme and using port 8866 by default.

``plasma-calculator`` takes optional arguments such as :option:`--dark`,
:option:`--port`, and :option:`--no-browser`. Pass flag :option:`-h` or
:option:`--help` to get the full list of supported arguments.
