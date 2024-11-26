<div align="center"><img src="https://raw.githubusercontent.com/PlasmaPy/PlasmaPy-logo/main/exports/with-text-dark.png" width="600"/></div>

# PlasmaPy

[![PyPI version](https://img.shields.io/pypi/v/plasmapy?style=flat&logo=pypi)](https://pypi.org/project/plasmapy/)
[![Conda version](https://img.shields.io/conda/v/conda-forge/plasmapy?style=flat&logo=anaconda)](https://img.shields.io/conda/v/conda-forge/plasmapy)
[![PyPI version](https://img.shields.io/pypi/pyversions/plasmapy?style=flat&logo=python)](https://img.shields.io/pypi/pyversions/plasmapy?style=plastic)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](./LICENSE.md)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](https://docs.plasmapy.org/en/latest/CODE_OF_CONDUCT.html)

[![Matrix](https://img.shields.io/badge/Matrix-join%20chat-blueviolet?style=flat&logo=matrix)](https://app.element.io/#/room/#plasmapy:openastronomy.org)
<a rel="me" href="https://fosstodon.org/@plasmapy">![Mastodon](https://img.shields.io/badge/Mastodon-plasmapy%40fosstodon.org-blue?logo=mastodon&style=fla)</a>
[![YouTube](https://img.shields.io/badge/YouTube%20-subscribe-red?style=flat&logo=youtube)](https://www.youtube.com/channel/UCSH6qzslhqIZKTAJmHPxIxw)

[![CI](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci.yml)
[![weekly tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly.yml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/PlasmaPy/PlasmaPy/main.svg)](https://results.pre-commit.ci/latest/github/PlasmaPy/PlasmaPy/main)
[![codecov](https://codecov.io/gh/PlasmaPy/PlasmaPy/branch/main/graph/badge.svg)](https://codecov.io/gh/PlasmaPy/PlasmaPy)
[![Read the Docs Status](https://readthedocs.org/projects/plasmapy/badge/?version=latest)](http://plasmapy.readthedocs.io/en/latest/?badge=latest)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1436011.svg)](https://doi.org/10.5281/zenodo.1436011)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat&logo=astropy)](http://www.astropy.org/)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)

[Astropy]: https://www.astropy.org
[authors and credits]: https://docs.plasmapy.org/en/latest/about/credits.html
[3-clause BSD license]: ./LICENSE.md
[calendar]: https://calendar.google.com/calendar/embed?src=c_sqqq390s24jjfjp3q86pv41pi8%40group.calendar.google.com&ctz=America%2FNew_York
[cite PlasmaPy]: https://docs.plasmapy.org/en/latest/about/citation.html
[code of conduct]: http://docs.plasmapy.org/en/latest/CODE_OF_CONDUCT.html
[community meetings]: https://www.plasmapy.org/meetings/weekly
[contributor guide]: https://docs.plasmapy.org/en/latest/development/index.html
[Department of Energy]: https://www.energy.gov
[example gallery]: https://docs.plasmapy.org/en/stable/examples.html
[GitHub discussions]: https://github.com/PlasmaPy/PlasmaPy/discussions
[Gitter]: https://gitter.im/PlasmaPy/Lobby
[good first issues]: https://github.com/PlasmaPy/PlasmaPy/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22
[Google Summer of Code]: https://summerofcode.withgoogle.com
[**install PlasmaPy**]: https://docs.plasmapy.org/en/stable/install.html
[GitHub repository]: https://github.com/PlasmaPy/PlasmaPy
[mailing list]: https://groups.google.com/forum/#!forum/plasmapy
[Matrix]: https://app.element.io/#/room/#plasmapy:openastronomy.org
[meetings]: https://www.plasmapy.org/meetings/weekly
[NASA]: https://www.nasa.gov/
[National Science Foundation]: https://nsf.gov
[office hours]: http://www.plasmapy.org/meetings/office_hours
[PlasmaPy Community on Zenodo]: https://zenodo.org/communities/plasmapy
[PlasmaPy]: https://www.plasmapy.org
[**documentation**]: https://docs.plasmapy.org
[protections against software patents]: ./PATENT.md
[Python]: https://www.python.org
[Smithsonian Institution]: https://www.si.edu
[submit a bug report]: https://github.com/PlasmaPy/PlasmaPy/issues/new?assignees=&labels=Bug&template=bug_report.yml
[submit a feature request]: https://github.com/PlasmaPy/PlasmaPy/issues/new?assignees=&labels=Feature+request&template=feature_request.yml
[team@plasmapy.org]: mailto:team@plasmapy.org
[this video]: https://youtu.be/E8RwQF5wcXM
[Zoom]: https://zoom.us/j/91633383503?pwd=QWNkdHpWeFhrYW1vQy91ODNTVG5Ndz09

[PlasmaPy] is an open source, community-developed [Python] package for
plasma research and education. PlasmaPy intends to be for plasma
science what [Astropy] is for astronomy — a collection of
functionality commonly needed by plasma scientists and researchers
globally, running within and leveraging the open source scientific
Python ecosystem. The goals of PlasmaPy are more thoroughly described
in [this video]. Many of our recent presentations are available from
the [PlasmaPy Community on Zenodo].

## Documentation

Please check out our online [**documentation**] to learn more about
PlasmaPy's capabilities.

If you would like an idea of what PlasmaPy can do, go to our [example
gallery] of Jupyter notebooks. To learn more about how to contribute,
check out PlasmaPy's [contributor guide].

## Installing PlasmaPy

PlasmaPy's online documentation has detailed instructions on how to
[**install PlasmaPy**].

To install PlasmaPy on macOS or Linux, open a terminal and run:
```Shell
python -m pip install plasmapy
```

On some systems, it might be necessary to specify the Python version
number, for example by using `python3` or `python3.13` instead of
`python`.

To install PlasmaPy on Windows, open a terminal and run
```Shell
py -3.13 -m pip install plasmapy
```
The `3.13` may be replaced by any version of Python that is installed
and supported by PlasmaPy.

## Citing PlasmaPy

If you use PlasmaPy for research resulting in a publication, please
[cite PlasmaPy]. It really helps support the project! Citing software
used in research provides credit to its authors, promotes open science &
scientific reproducibility, and helps open source projects demonstrate
to funding agencies that continued development should be supported.

Please check out the [PlasmaPy community on Zenodo] for prior releases
of PlasmaPy and other resources.

## Requesting features

Please [submit a feature request] in our [GitHub repository] if you
have an idea for new or improved functionality. PlasmaPy is
community-driven, and feature requests really help guide the future of
the project.

## Submitting bug reports

Please [submit a bug report] on PlasmaPy's GitHub repository if you
notice any problems. We really appreciate it!

## Contributing

If you are interested in contributing, please check out our [contributor
guide] and [code of conduct]. There are a number of [good first issues]
in our GitHub repository. New contributors are very welcome!

## Events

PlasmaPy has several [meetings] that are on our [calendar]. Events are
usually held on PlasmaPy's [Zoom] room.

Last-minute changes are usually announced on the [Matrix]/[Gitter]
chat room. The most up-to-date information about these meetings is on
the [meetings] page of PlasmaPy's website.

### Office hours

Our weekly informal [office hours] are an opportunity to chat with
active members of the PlasmaPy community about topics related to
Python and plasma physics. If you'd like to learn more about PlasmaPy,
our office hours are one of the best places to start. As of July 2024,
our office hours are on most Thursdays at 3 pm Eastern. Please feel
free to come by!

### Community meetings

PlasmaPy's weekly [community meetings] are a place to talk about code
development. If you have an idea for a new feature or would like to
make a code contribution, community meetings are a good place to go
to. As of July 2024, our community meetings are on most Tuesdays at 2 pm
Eastern.

<!--
### Project meetings

PlasmaPy's weekly project meetings are a place to discuss education,
outreach, and project coordination. Topics might range from creating
educational notebooks to organizing community events. As of July
2024, project meetings are held on most Wednesdays at 3 pm Eastern.
-->

<!--
### Working group meetings

PlasmaPy has started several working groups, including on diagnostics,
dispersion relations, and simulation. These working groups usually
meet fortnightly, and their meeting times can be found in PlasmaPy's
event [calendar]. If you would like to join a PlasmaPy working group
or even start a new one, please email us at [team@plasmapy.org]!
-->

## Community

### Matrix chat

If you have any questions, the quickest way to get a response is to
ask on our [Matrix]/[Gitter] channel. Both of these are the same chat
channel; Gitter uses a bridge to link the two.

### GitHub discussions

We're trying out [GitHub discussions] as a place to suggest ideas,
bring up discussion topics, and ask questions.

### Mailing list

You can subscribe to PlasmaPy's low-volume [mailing list] to receive
PlasmaPy newsletters and other announcements.

## Contact information

Please feel free to reach out to us at [team@plasmapy.org] or stop by
our [office hours] with any ideas, questions, and/or puns about
computational magnetohydrodynamics.

Please use these links to [submit a feature request] and to [submit a
bug report] on PlasmaPy's GitHub repository.

## License

PlasmaPy is permissively licensed under a [3-clause BSD license] with
added [protections against software patents].

## Acknowledgments

Development of PlasmaPy has been supported in part by the [National
Science Foundation], [Department of Energy], [NASA], and the
[Smithsonian Institution]. For more details, please see PlasmaPy's
documentation page on [authors and credits].
