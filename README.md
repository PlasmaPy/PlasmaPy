# PlasmaPy

[![Build Status](https://travis-ci.org/PlasmaPy/PlasmaPy.svg?branch=master)](https://travis-ci.org/PlasmaPy/PlasmaPy)
[![Coverage Status](https://coveralls.io/repos/github/PlasmaPy/PlasmaPy/badge.svg?branch=master)](https://coveralls.io/github/PlasmaPy/PlasmaPy?branch=master)

A community developed Python 3.6+ package for plasma physics in 
the early stages. PlasmaPy intends to be for plasmas what
[Astropy](https://github.com/astropy/astropy) is for astronomy - a 
collection of functionality commonly used and shared between plasma physicists 
and researchers globally, running within and leveraging the open source 
scientific Python ecosystem. The goals of this project are better described in our
[vision statement](https://github.com/PlasmaPy/PlasmaPy/blob/master/vision_statement.md)
and [an earlier conference poster](http://doi.org/10.5281/zenodo.163752).

We recently created a guide on
[contributing to PlasmaPy](https://github.com/PlasmaPy/PlasmaPy/blob/master/CONTRIBUTING.md).
New contributors are very welcome! 

If you have any questions, please send us a message at our
[Riot channel](https://riot.im/app/#/room/#plasmapy:matrix.org) or contact
Nick Murphy at <namurphy@cfa.harvard.edu>,
Yi-Min Huang at <yiminh@princeton.edu>,
or Drew Leonard at <andy.j.leonard@gmail.com>.

## Installation

There are multiple options to download the source code for PlasmaPy.
The simplest is to select "Clone or Download" on our 
[repository](https://github.com/PlasmaPy/PlasmaPy) page.  This will provide 
an option to download a zip file plus information on how to 
clone the repository.  If you have git installed on your computer and you
would like to use HTTPS (which is the default and easier to set up), then run:

```ShellSession
git clone https://github.com/PlasmaPy/PlasmaPy.git
```

If you have [set up an SSH key](https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/),
an equivalent and more secure command is:

```ShellSession
git clone git@github.com:PlasmaPy/PlasmaPy.git
```

The [contributing to PlasmaPy](https://github.com/PlasmaPy/PlasmaPy/blob/master/CONTRIBUTING.md)
guide has instructions on how to fork a repository so that you may make pull requests.

In the top level directory, run

```ShellSession
pip install .
```
or
```ShellSession
python setup.py install
```

PlasmaPy is [not compatible with Python 2](https://pythonclock.org/).

## License

PlasmaPy is licensed under a 3-clause BSD style license - see
``LICENSE.md`` file in the top-level directory.
