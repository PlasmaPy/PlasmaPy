.. _plasmapy-calculator:

=================
Plasma Calculator
=================

.. currentmodule:: plasmapy.utils.calculator

.. attention::

   |expect-api-changes|

Overview
--------

The prototype Plasma Calculator is an interactive Jupyter notebook that
allows users to input a set of plasma properties and immediately
calculate multiple plasma parameters.

.. note::

   This functionality is still under development and the API may change
   in future releases.

Installing the Plasma Calculator
--------------------------------

To install the dependencies needed to run the plasma calculator, it is
necessary to specify the ``calculator`` dependency set when installing
`plasmapy`. This step can be done, for example, with the following
command:

.. code-block:: bash

   pip install plasmapy[calculator]

Using Plasma Calculator
-----------------------

To invoke the app use ``plasma-calculator`` in the command line. By
default this opens the app in a browser, bootstrapping the notebook in a
light theme and using port 8866 by default.

``plasma-calculator`` takes optional arguments such as ``--dark``,
``--port``, and ``--no-browser``. Pass flag ``-h`` or ``--help`` to get
full list of supported arguments.
