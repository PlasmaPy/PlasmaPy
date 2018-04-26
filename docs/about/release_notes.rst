.. _release-notes:

=============
Release Notes
=============

This document contains the release notes for each version of PlasmaPy.
A list of the changes for each version are included in the
:ref:`change_log_`.

PlasmaPy uses `semantic versioning <http://www.semver.org/>`_.  Version
numbers are of the form `MAJOR.MINOR.PATCH`.  Development releases have
`MAJOR` equal to `0`.  The `application programming interface (API)
<https://en.wikipedia.org/wiki/Application_programming_interface>`_
during the development phase is unstable and anything may change at
any time.  Starting with version `1.0.0`, `MAJOR` must be incremented
whenever the API changes in a way that is not backwards compatible,
`MINOR` must be incremented whenever the API expands while maintaining
backwards compatibility, and `PATCH` must be incremented for each bug
fix release.

Version `0.1.0`
---------------

We are excited to announce the first development release of `PlasmaPy
<http://www.plasmapy.org/>`_: a community-developed fully open source
core Python package for plasma physics.

Version `0.1.0` is a preview and a prototype.  It is not yet feature
complete or recommended for production work.  Significant changes to the
API are expected to occur between versions `0.1.0` and `0.2.0`.  Rather,
version `0.1.0` serves as an invitation to plasma students and
scientists to collaboratively develop a community-wide shared software
package for our field.

If you have a scientific Python 3.6+ environment already configured,
you may install PlasmaPy with `pip <https://pypi.org/project/pip/>`_ by
running:

.. code-block::
    pip install plasmapy

We recommend installing PlasmaPy in an Anaconda environment. `PlasmaPy's
GitHub repository <https://github.com/PlasmaPy/PlasmaPy>`_ contains
more detailed `installation instructions
<https://github.com/PlasmaPy/PlasmaPy/blob/master/INSTALL.md>`_.

PlasmaPy uses Astropy's `~astropy.units` subpackage to represent
physical quantities with units, and for compatibility with
`Astropy <http://www.astropy.org/>`_ and `SunPy <http://sunpy.org/>`_.
This subpackage handles unit conversions, and raises exceptions for
operations that have incompatible units.  New users may wish to become
familiar with this functionality by reading Astropy's `units subpackage
documentation <http://docs.astropy.org/en/stable/units/>`_.  An example
use case is:

    >>> import astropy.units as u
    >>> 88 * u.imperial.mi / u.hr
    <Quantity 88. mi / h>
    >>> (1.21e9 * u.J / u.s).to(u.GW)
    <Quantity 1.21 GW>

PlasmaPy's `~plasmapy.physics` subpackage contains functions to
calculate a wide variety of plasma parameters, dielectric tensor
components, and relativity/quantum physics parameters used in plasma
physics.  The `~plasmapy.physics.transport` module of
`~plasmapy.physics` contains functionality to calculate collision rates
and transport parameters (including an object-oriented interface to
classical transport coefficients).  The `~plasmapy.atomic` subpackage
includes both functional and object-oriented interfaces to access atomic
parameters and represent particles. The `~plasmapy.mathematics`
subpackage contains analytical functions that are commonly used in
plasma physics (including the plasma dispersion function).  The
`~plasmapy.classes` subpackage includes prototype classes to represent
plasma configurations, including a particle pusher.

PlasmaPy requires Python 3.6+.  The core developers chose to
support Python 3.6+ because Python 2.7 will cease to be supported by
most scientific Python packages within about a year, Python 3.6 will
likely to be the oldest version of Python still in common use by the
time we release PlasmaPy `1.0.0`, and Python 3.6 contains new features
such as formatted string literals that greatly improve readability.

If there is functionality that you would like future versions of
PlasmaPy to include or if you discover a bug, we encourage you to
`raise an issue <https://github.com/PlasmaPy/PlasmaPy/issues/new>`_ with
your ideas or even contribute code directly.

The following resources provide more information on PlasmaPy, including
how to contribute.

* `PlasmaPy's online documentation <docs.plasmapy.org>`_
* `PlasmaPy's GitHub repository <https://github.com/PlasmaPy/PlasmaPy>`_
* A guide on :ref:`contributing-to-plasmapy`
* :ref:`subpackage-stability`
* :ref:`plasmapy-vision-statement`
* `PlasmaPy's website <http://www.plasmapy.org/>`_
* :ref:`plasmapy-code-of-conduct`

This release includes over 1800 commits and 173 merged pull requests,
with contributions from 35 different people to the code base or the
vision statement.

PlasmaPy's Coordinating Committee includes:

* Drew Leonard
* Nick Murphy
* Tulasi Parashar
* Dominik Stańczak

The people who have contributed to this release include:

* Jasper Beckers
* Ludovico Bessi
* Sean Carroll
* Apoorv Choubey
* cclauss
* Leah Einhorn
* Thomas Fan
* Graham Goudeau
* Silvina Guidoni
* Colby Haggerty
* Julien Hillairet
* Poh Zi How
* Yi-Min Huang
* Nabil Humphrey
* Maria Isupova
* Pawel Kozlowski
* Siddharth Kulshrestha
* Piotr Kuszaj
* Samuel Langendorf
* Drew Leonard
* Ritiek Malhotra
* Stuart Mumford
* Joshua Munn
* Nick Murphy
* Nismirno
* nrb1324
* Tulasi Parashar
* Neil Patel
* Roberto Díaz Pérez
* Raajit Raj
* Dawa Nurbu Sherpa
* David Stansby
* Dominik Stańczak
* Antoine Tavant
* Sixue Xu
