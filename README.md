# PlasmaPy

[![Build Status](https://travis-ci.org/PlasmaPy/PlasmaPy.svg?branch=master)](https://travis-ci.org/PlasmaPy/PlasmaPy) [![Coverage Status](https://coveralls.io/repos/github/PlasmaPy/PlasmaPy/badge.svg?branch=master)](https://coveralls.io/github/PlasmaPy/PlasmaPy?branch=master)

A community developed Python package for plasma physics.

## Project Status

PlasmaPy is in an early stage of development.  The goals of this project are described in our [vision statement](https://github.com/PlasmaPy/PlasmaPy/blob/master/vision_statement.md) and an earlier conference poster:

* Murphy, Nicholas A, Huang, Yi-Min, & PlasmaPy Community. (2016, October). PlasmaPy: beginning a community developed Python package for plasma physics. Zenodo. http://doi.org/10.5281/zenodo.163752

We recently created a guide on [contributing to PlasmaPy](https://github.com/PlasmaPy/PlasmaPy/blob/master/CONTRIBUTE.md), which also contains instructions on how to join our email list.  New contributors are very welcome!  

If you have any questions, please contact Nick Murphy at <namurphy@cfa.harvard.edu>, Yi-Min Huang at <yiminh@princeton.edu>, or Drew Leonard at <andy.j.leonard@gmail.com>.

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

If you have [set up an SSH key](https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/), an equivalent and more secure command is:

```ShellSession
git clone git@github.com:PlasmaPy/PlasmaPy.git
```

The [contributing to PlasmaPy](https://github.com/PlasmaPy/PlasmaPy/blob/master/CONTRIBUTE.md) guide has instructions on how to fork a repository so that you may make pull requests.

In the top level directory, run

```ShellSession
python setup.py install
```

PlasmaPy is presently being designed to be compatible with Python 3.6 and above, and does not guarantee support for Python 3.5 and below.  PlasmaPy is not compatible with Python 2.

## License

PlasmaPy is licensed under a 3-clause BSD style license - see
``LICENSE.md`` file in the top-level directory.
