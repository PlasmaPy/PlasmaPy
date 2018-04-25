.. _release-notes:

======================
PlasmaPy Release Notes
======================

This document contains the release notes for each version of PlasmaPy.
A list of the changes for each version are included in the :ref:`change_log_`.

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

Version `0.1.0` (2018-04-27)
--------------------------
We are pleased to announce the first development release of
`PlasmaPy <http://www.plasmapy.org/>`_: a community-developed fully open
source core Python package for plasma physics.

Version `0.1.0` is a preview and a prototype.  It is not yet feature
complete or recommended for production work.  Significant changes to the
API are expected to occur between versions `0.1.0` and `0.2.0` and
throughout the development phase. Rather, version `0.1.0` serves as an
invitation to plasma physicists to collaboratively develop a
community-wide shared software package for our field.

.. What needs to be included still?
   Link to vision statement and code of conduct.
   Requirements
   Link to how to install
   Link to doc page for each subpackage

PlasmaPy is compatible with Python 3.6.  The core developers chose to
support Python 3.6+ because Python 2.7 will cease to be supported by
most scientific Python packages within about a year, Python 3.6 is
likely to be the oldest version of Python still in use by the time we
release PlasmaPy `1.0.0`, and (3) Python 3.6 contains new features such
as formatted string literals that greatly improve readability.

PlasmaPy uses the `~astropy.units` package for compatibility with
`Astropy <http://www.astropy.org/>`_ and `SunPy <http://sunpy.org/>`_.

If there is functionality that you would like future versions of
PlasmaPy to include, we encourage you to
`raise an issue <https://github.com/PlasmaPy/PlasmaPy/issues/new>`_ with
your ideas or even contribute code directly.

This release includes roughly 1800 commits with contributions from about
35 different contributors to the code base and vision statement.  The
people who have contributed to this release include:

* Jasper Beckers
* Ludovico Bessi
* Sean Carroll
* Apoorv Choubey
* C. Clauss
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
