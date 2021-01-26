<div align="center"><img src="https://raw.githubusercontent.com/PlasmaPy/PlasmaPy-logo/master/exports/with-text-dark.png" width="600"/></div>

# PlasmaPy

[![PyPI version](https://img.shields.io/pypi/v/plasmapy?style=flat&logo=pypi)](https://pypi.org/project/plasmapy/)
[![Conda version](https://img.shields.io/conda/v/conda-forge/plasmapy?style=flat&logo=anaconda)](https://img.shields.io/conda/v/conda-forge/plasmapy)
[![PyPI version](https://img.shields.io/pypi/pyversions/plasmapy?style=flat&logo=python)](https://img.shields.io/pypi/pyversions/plasmapy?style=plastic)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](./LICENSE.md)

[![Matrix](https://img.shields.io/badge/Matrix-join%20chat-blueviolet?style=flat&logo=matrix)](https://app.element.io/#/room/#plasmapy:openastronomy.org)
[![YouTube](https://img.shields.io/badge/Twitter%20-follow-red?style=flat&logo=twitter)](https://www.youtube.com/channel/UCSH6qzslhqIZKTAJmHPxIxw)
[![YouTube](https://img.shields.io/badge/YouTube%20-subscribe-red?style=flat&logo=youtube)](https://www.youtube.com/channel/UCSH6qzslhqIZKTAJmHPxIxw)

[![GitHub Actions — CI](https://github.com/PlasmaPy/PlasmaPy/workflows/CI/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions?query=workflow%3ACI+branch%3Amaster)
[![GitHub Actions — Style linters](https://github.com/PlasmaPy/PlasmaPy/workflows/Style%20linters/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions?query=workflow%3AStyle-linters+branch%3Amaster)
[![codecov](https://codecov.io/gh/PlasmaPy/PlasmaPy/branch/master/graph/badge.svg)](https://codecov.io/gh/PlasmaPy/PlasmaPy)
[![Read the Docs Status](https://readthedocs.org/projects/plasmapy/badge/?version=latest&logo=twitter)](http://plasmapy.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PlasmaPy/PlasmaPy/master?filepath=docs/notebooks)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1436011.svg)](https://doi.org/10.5281/zenodo.1436011)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat&logo=astropy)](http://www.astropy.org/)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Open Source Helpers](https://www.codetriage.com/plasmapy/plasmapy/badges/users.svg)](https://www.codetriage.com/plasmapy/plasmapy)

[PlasmaPy](https://www.plasmapy.org/) is an open source, community-developed
Python 3.7+ package for plasma science. PlasmaPy intends to be for plasma
science what [Astropy](https://github.com/astropy/astropy) is for astronomy
— a collection of functionality commonly used and shared between plasma
scientists and researchers globally, running within and leveraging the
open source scientific Python ecosystem. The goals of this project are
more thoroughly described in [this recent video](https://youtu.be/E8RwQF5wcXM).
Current functionality is described in [PlasmaPy's online
documentation](http://docs.plasmapy.org/en/latest/).

We created a guide on [contributing to PlasmaPy](http://docs.plasmapy.org/en/stable/CONTRIBUTING.html)
and have a [code of conduct](http://docs.plasmapy.org/en/stable/CODE_OF_CONDUCT.html).
New contributors are very welcome!

# Installation

If you have [installed Python](https://wiki.python.org/moin/BeginnersGuide/Download),
you can install PlasmaPy from [pip](https://pypi.org/project/pip/)
via
```Shell
python -m pip install plasmapy
```
If you have
[installed conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html),
then you can also get PlasmaPy from
```Shell
conda install -c conda-forge plasmapy
```
To contribute
to the package, check out [our instructions on installing PlasmaPy from
source](http://docs.plasmapy.org/en/stable/install.html#building-and-installing-from-source-code).

# Community

## [Matrix chat](https://app.element.io/#/room/#plasmapy:openastronomy.org)

If you have any questions, the quickest way to get a response is to ask
on our
[Matrix](https://app.element.io/#/room/#plasmapy:openastronomy.org)/[Gitter](https://gitter.im/PlasmaPy/Lobby)
channel. Both of these are the same chat channel; Gitter uses a bridge to link the two.

## [Weekly](https://calendar.google.com/calendar?cid=bzVsb3ZkcW0zaWxsam00ZTlrMDd2cmw5bWdAZ3JvdXAuY2FsZW5kYXIuZ29vZ2xlLmNvbQ) [community meetings](https://meet.jit.si/plasmapy)

We have weekly community meetings in the
[PlasmaPy room on Jitsi](https://meet.jit.si/plasmapy).
The schedule of our community meetings is on our [calendar](https://calendar.google.com/calendar?cid=bzVsb3ZkcW0zaWxsam00ZTlrMDd2cmw5bWdAZ3JvdXAuY2FsZW5kYXIuZ29vZ2xlLmNvbQ), and you may access the [minutes and
agendas](https://drive.google.com/drive/folders/0ByPG8nie6fTPV1FQUEkzMTgtRTg?usp=sharing).
Any last minute changes will be discussed on
[Matrix](https://app.element.io/app/#/room/#plasmapy:openastronomy.org).
As of January 2021, our meetings are on Tuesdays at
[19:00 UTC](http://time.unitarium.com/utc/6pm).
Come discuss plasma software with us!

## [Weekly office hours](http://www.plasmapy.org/meetings/office_hours/)

PlasmaPy's weekly [office hours](http://www.plasmapy.org/meetings/office_hours/)
on Thursdays at [19:00 UTC](http://time.unitarium.com/utc/6pm)
are an opportunity to chat with active members of the PlasmaPy
community about the package and project.

## [GitHub discussions](https://github.com/PlasmaPy/PlasmaPy/discussions)

We're now trying out GitHub discussions for more varied topics that aren't
exactly issues with the existing code base. It's a great place to suggest
ideas, bring up discussion topics, and ask questions.

## [Mailing list](https://groups.google.com/forum/#!forum/plasmapy)

You can subscribe to our low-volume
[mailing list](https://groups.google.com/forum/#!forum/plasmapy)
to receive PlasmaPy newsletters and other announcements.

## [Suggestion box](https://docs.google.com/forms/d/e/1FAIpQLSdT3O5iHZrLJRuavFyzoR23PGy0Prfzx2SQOcwJGWtvHyT2lw/viewform?usp=sf_link)

We have
[a suggestion box](https://docs.google.com/forms/d/e/1FAIpQLSdT3O5iHZrLJRuavFyzoR23PGy0Prfzx2SQOcwJGWtvHyT2lw/viewform?usp=sf_link)
if you would like to (optionally anonymously) suggest
a feature/topic for consideration. These will be reposted on the mailing list
or directly in GitHub issues, as appropriate, for further discussion.

# License

PlasmaPy is permissively licensed under a
[3-clause BSD license with added protections
against software patents](LICENSE.md).

# Citing PlasmaPy

An [emerging best practice for software
citation](https://doi.org/10.7717/peerj-cs.86) is to cite the _specific
version_ of each software package used in a research project (instead of
only citing a journal article, website, or GitHub repository). The
citation should include a persistent identifier that uniquely identifies
which version of the software was used. We therefore ask that you cite
the specific version of PlasmaPy used in your research project. Releases
of PlasmaPy are available in the
[PlasmaPy community](https://zenodo.org/communities/plasmapy) on
[Zenodo](https://zenodo.org/), along with many other PlasmaPy resources.
Please check our documentation for more detailed [citation
instructions](./docs/about/citation.rst).

# Acknowledgements

Early development on PlasmaPy was supported in part by the U.S.
Department of Energy, the Smithsonian Institution, and Google Summer of
Code. Ongoing PlasmaPy development is being supported through a
collaborative award from the U.S. National Science Foundation's
Cyberinfrastructure for Sustained Scientific Innovation program and a
NASA Heliophysics Data Environment Enhancements award.
