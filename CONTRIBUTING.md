# Contributing to PlasmaPy

There are numerous ways to contribute to PlasmaPy, including by
providing code and documentation, suggesting and discussing ideas,
submitting issues and bug reports, and engaging the broader plasma
physics community.  

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

*This disclaimer was originally written by
[Adrienne Lowe](https://github.com/adriennefriend) for a
[PyCon talk](https://www.youtube.com/watch?v=6Uj746j9Heo), and was adapted by 
[yt](https://github.com/yt-project/yt) in their README file based on its use 
in the README file for the [MetPy project](https://github.com/Unidata/MetPy).
It was then adapted by PlasmaPy.*

## Sharing ideas

There are several methods of communication that are being used in the
early stages of PlasmaPy development:

* [Signing up for the PlasmaPy email
  list](https://groups.google.com/forum/#!forum/plasmapy) will allow
  you to participate in broader discussions and keep up with the
  latest developments.

* The [PlasmaPy repository on
  GitHub](https://github.com/PlasmaPy/plasmapy) is the best place to
  [submit issues](https://github.com/PlasmaPy/plasmapy/issues) and
  review [pull requests](https://github.com/PlasmaPy/plasmapy/pulls).

* The [PlasmaPy Matrix channel](https://riot.im/app/#/room/#plasmapy:matrix.org)
  is a great place to have informal conversations, coordinate efforts,
  and share ideas.  

* We have biweekly telecons which are announced on the email list.

## Contributing code or documentation to PlasmaPy

### Preliminaries

Before contributing to the PlasmaPy code base, one must [**join
GitHub**](https://github.com/join?source=header-home).  A free account
will suffice for you to have unlimited public repositories.  If you
are new to [git](https://git-scm.com/), helpful resources include
documentation on [git
basics](https://git-scm.com/book/en/v2/Getting-Started-Git-Basics) and
an [interactive git
tutorial](https://try.github.io/levels/1/challenges/1).  You must also
[install
git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
locally on your computer.  We highly recommend getting familiar with
git by going through these tutorials or a [Software
Carpentry](https://software-carpentry.org/) workshop prior to making
code contributions.

### Forking and cloning PlasmaPy

After creating your GitHub account, go to the [main
repository](https://github.com/PlasmaPy/PlasmaPy) and **fork a copy of
PlasmaPy to your account**.

Next you must **clone your fork to your computer**.  Go to the
directory that will host your PlasmaPy directory, and run one of the
following commands (after changing *your-username* to your username).
If you would like to use HTTPS (which is the default and easier to set
up), then run:

```ShellSession
git clone https://github.com/your-username/PlasmaPy.git
```

SSH is a more secure option, but requires you to [set up an SSH
key](https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/) beforehand.  The equivalent SSH command is:

```ShellSession
git clone git@github.com:your-username/PlasmaPy.git
```

After cloning, we must tell git where the development version of
PlasmaPy is by running:

```ShellSession
git remote add upstream git://github.com/PlasmaPy/PlasmaPy.git
```

To check on which remotes exist, run `git remote -v`.  You should get
something like this:

```ShellSession
origin		git@github.com:namurphy/PlasmaPy.git (fetch)
origin		git@github.com:namurphy/PlasmaPy.git (push)
upstream	git@github.com:PlasmaPy/PlasmaPy.git (fetch)
upstream	git@github.com:PlasmaPy/PlasmaPy.git (push)
```

### Branches, commits, and pull requests

Before making any changes, it is prudent to update your local
repository with the most recent changes from the development
repository:

```ShellSession
git fetch upstream
```

Changes to PlasmaPy should be made using branches.  It is usually best
to avoid making changes on your master branch so that it can be kept
consistent with the upstream repository.  Instead we can create a new
branch for the specific feature that you would like to work on:

```ShellSession
git branch *your-new-feature*
``` 

Descriptive branch names such as `grad-shafranov` or
`adding-eigenfunction-poetry` are helpful, while vague names like
`edits` are considered harmful.  After creating your branch locally,
let your fork of PlasmaPy know about it by running:

```ShellSession
git push --set-upstream origin *your-new-feature*
``` 

It is also useful to configure git so that only the branch you are
working on gets pushed to GitHub:

```ShellSession
git config --global push.default simple
```

Once you have set up your fork and created a branch, you are ready to
make edits to PlasmaPy.  Switch to your new branch by running:

```ShellSession
git checkout *your-new-feature*
```

Go ahead and modify files with your favorite text editor.  Be sure to
include tests and documentation with any new functionality.  We also
recommend reading about [best practices for scientific
computing](https://doi.org/10.1371/journal.pbio.1001745).  PlasmaPy
uses the [PEP 8 style guide for Python
code](https://www.python.org/dev/peps/pep-0008/) and the [numpydoc
format for
docstrings](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt)
to maintain consistency and readability.  New contributors should not 
worry too much about precisely matching these styles when first 
submitting a pull request, as the [PEP8 Speaks](http://pep8speaks.com/)
GitHub integration will check pull requests for PEP 8 compatibility, and
further changes to the style can be suggested during code review.

You may periodically commit changes to your branch by running

```ShellSession
git add filename.py
git commit -m "*brief description of changes*"
```

Committed changes may be pushed to the corresponding branch on your
GitHub fork of PlasmaPy using 

```ShellSession
git push origin *your-new-feature* 
```

or, more simply,

```ShellSession
git push
```

Once you have completed your changes and pushed them to the branch on
GitHub, you are ready to make a pull request.  Go to your fork of
PlasmaPy in GitHub.  Select "Compare and pull request".  Add a
descriptive title and some details about your changes.  Then select
"Create pull request".  Other contributors will then have a chance to
review the code and offer contructive suggestions.  You can continue
to edit the pull request by changing the corresponding branch on your
PlasmaPy fork on GitHub.  After a pull request is merged into the
code, you may delete the branch you created for that pull request.
