# Installation

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

The [contributing to PlasmaPy](http://docs.plasmapy.org/en/master/CONTRIBUTING.html)
guide has instructions on how to fork a repository so that you may make pull requests.

In the top level directory, run

```ShellSession
pip install .
```
or
```ShellSession
python setup.py install
```
****

We are officially [on PyPI](https://pypi.org/project/plasmapy/)
- PlasmaPy can be installed via

```ShellSession
pip install plasmapy
```

Though be warned that the version on PyPI is our distribution version and not
suitable for development. If you wish to contribute to the project, please
install from GitHub.

We're not on conda just yet, but we're working on it!
