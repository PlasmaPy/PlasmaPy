======================
PlasmaPy Release Notes
======================

This document contains the release notes for each version of PlasmaPy.
A list of the changes for each version are included in the [change
log](https://github.com/PlasmaPy/PlasmaPy/blob/master/CHANGE_LOG.md).

PlasmaPy uses [semantic versioning](http://semver.org/), as described
in [PlasmaPy Enhancement Proposal (PLEP)
5](https://github.com/PlasmaPy/PlasmaPy-PLEPs/blob/master/PLEP-0005.md).
Version numbers are of the form `MAJOR.MINOR.PATCH`.  Development
releases have `MAJOR` equal to zero.  The [application programming
interface
(API)](https://en.wikipedia.org/wiki/Application_programming_interface)
during the development phase is unstable and anything may change at
any time.  Starting with version 1.0.0, `MAJOR` must be incremented
whenever the API changes in a way that is not backwards compatible,
`MINOR` must be incremented whenever the API expands while maintaining
backwards compatibility, and `PATCH` must be incremented for each bug
fix release.

Version 0.1.0 (2018-xx-xx)
--------------------------

We are pleased to announce the first development release of PlasmaPy:
a community-developed fully open source core Python package for plasma
physics.

Version 0.1.0 is a preview and a prototype.  It is not yet feature
complete or recommended for production work.  Rather, version 0.1.0
serves as an invitation to plasma scientists to participate in the
development of shared software for our community.

PlasmaPy uses the `astropy.units` package for compatibility with
Astropy and SunPy.

Significant changes to the API are expected to occur between
versions 0.1.0 and 0.2.0 and through the development phase.

.. What needs to be included still?
   Link to vision statement and code of conduct.
   Requirements
   Link to how to install
   Link to doc page for each subpackage

PlasmaPy is being designed to be compatible with Python 3.6. The
required packages include...

If there is functionality that you would like future versions of
PlasmaPy to include, we encourage you to raise an issue on GitHub with
your ideas or even contribute code directly.

This release includes roughly 1400 commits by 31 different
contributors. The people who have contributed to this release are:

 * ...