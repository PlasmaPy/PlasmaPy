# A vision for an open source core Python package for plasma physics

## About PlasmaPy

PlasmaPy is a community-developed and community-driven free and open
source Python package that provides common functionality required for
plasma physics in a single, reliable codebase.

## Motivation

In recent years, researchers in many different scientific disciplines
have worked together to develop core Python packages such as Astropy,
SunPy, and SpacePy. These packages provide core functionality, common
frameworks for data visualization and analysis, and educational tools
for their respective scientific disciplines. We believe that a similar
cooperatively developed package for plasma physics will greatly
benefit our field. In this document, we lay out our vision for
PlasmaPy: a community-developed and community-driven open source core
Python software package for plasma physics.

There is considerable need in plasma physics for open, general purpose
software framework using modern [best practices for scientific
computing](http://dx.doi.org/10.1371/journal.pbio.1001745). As most
scientific programmers are largely self-taught, software often does not
take advantage of these practices and is instead written in a rush to
produce results for the next research paper. The resulting code is
often difficult to read and maintain, the documentation is usually
inadequate, and tests are typically implemented late in the
development process if at all.  Legacy code is often written in low
level languages such as Fortran, which typically makes compiling and
installing packages difficult and frustrating, especially if it
calls external libraries. It is also unusual to share code, and access
to major software is often restricted in some way, resulting in many
different programs and packages which do essentially the same thing but
with little or no interoperability. These factors lead to research
that is difficult to reproduce, and present a significant barrier to
entry for new users.

The plasma physics community is slowly moving in the open source
direction.  Several different types of packages and software have 
been released under open source licences, including the UCLA 
PIC codes, PICCANTE, EPOCH, VPIC, PIConGPU, WARP, the FLASH framework, 
Athena, and PENCIL.  These projects are built as individual packages, 
are written in different programming languages, and often have many 
dependencies  on specific packages.  Python packages such as Astropy, 
SunPy, and SpacePy have had notable success providing
open source alternatives to legacy code in related fields.  We are
grateful to these communities for their hard work, and hope to build
upon their accomplishments for the field of plasma physics.

An end user might not always be interested in a complicated powerpack
to perform one specific task on supercomputers. She might also be
interested in performing some basic plasma physics calculations,
running small desktop scale simulations to test preliminary ideas
(e.g., 1D MHD/PIC or test particles), or even comparing data from two
different sources (simulations vs. spacecraft). Such tasks require a
central platform. This is where PlasmaPy comes in.

## Ensuring a welcoming and inclusive environment

PlasmaPy strives to follow the best practices in open source software
development. New contributors are encouraged to join the team and
contribute to the codebase. We anticipate/encourage a global
participation from people with diverse backgrounds, skills, interests,
and opinions. We believe that such diversity is critical in ensuring a
growth of ideas in our community. We as a community pledge to abide by
the following guidelines:

* We pledge to treat all people with respect and provide a harassment-
  and bullying-free environment, regardless of sex, sexual orientation
  and/or gender identity, disability, physical appearance, body size,
  race, nationality, ethnicity, and religion. In particular, sexual
  language and imagery, sexist, racist, or otherwise exclusionary
  jokes are not appropriate.
* We pledge to respect the work of others by recognizing
  acknowledgment/citation requests of original authors. As authors, we
  pledge to be explicit about how we want our own work to be cited or
  acknowledged.
* We pledge to welcome those interested in joining the community, and
  realize that including people with a variety of opinions and
  backgrounds will only serve to enrich our community. In particular,
  discussions relating to pros/cons of various technologies,
  programming languages, and so on are welcome, but these should be
  done with respect, taking proactive measure to ensure that all
  participants are heard and feel confident that they can freely
  express their opinions.
* We pledge to welcome questions and answer them respectfully, paying
  particular attention to those new to the community. We pledge to
  provide respectful criticisms and feedback in forums, especially in
  discussion threads resulting from code contributions.
* We pledge to be conscientious of the perceptions of the wider
  community and to respond to criticism respectfully. We will strive
  to model behaviors that encourage productive debate and
  disagreement, both within our community and where we are
  criticized. We will treat those outside our community with the same
  respect as people within our community.
* We pledge to work from the very beginning of this project to make
  PlasmaPy accessible to people with disabilities.
* We pledge to help the entire community follow these guidelines, and
  to not remain silent when we see violations of them. We will take
  action when members of our community violate these
  guidelines. Members of the PlasmaPy community may contact any member
  of the Coordinating Committee to report violations. Members of the
  Coordinating Committee will treat these reports in the strictest
  confidence.  The Coordinating Committee will develop formal
  procedures for how to handle reported violations.

Parts of these guidelines have been adapted from the [Astropy code of
conduct](http://www.astropy.org/about.html#codeofconduct) and the
[Python Software Foundation code of
conduct](https://www.python.org/psf/codeofconduct/).

## Organizational structure

The Coordinating Committee (CC) will oversee the PlasmaPy project and
code development.  The CC will ensure that roles are being filled,
facilitate community-wide communication, coordinate and delegate
tasks, manage the project repository, oversee the code review process,
regulate intercompatibility between different subpackages, seek
funding mechanisms, facilitate compromises and cooperation, enforce
the code of conduct, and foster a culture of appreciation.

The Community Engagement Committee (CEC) will be responsible for
organizing conferences, trainings, and workshops; maintaining and
moderating social media groups and accounts; overseeing PlasmaPy’s
website; and communicating with the PlasmaPy and plasma physics
communities. The CEC will facilitate partnerships with groups such as
[Software Carpentry](https://software-carpentry.org/).

Each subpackage will have lead and deputy coordinators who will guide
and oversee the development of that subpackage.

The Accessibility Coordinator will work to ensure that the PlasmaPy
codebase, documentation, and practices are accessible to disabled
students and scientists.  Additional roles include the Webpage
Maintainer, the Release Coordinator, and the Testing Coordinator.

The work undertaken by each of these groups and coordinators should be
done openly and transparently, except where confidentiality is needed.
We will strive to have multiple subfields from plasma physics in each
committee.  Major decisions should ideally be made by general
consensus among the PlasmaPy community, but when consensus is not
possible then the committees may decide via majority vote.  Much of
this section is following the [organizational structure of
Astropy](http://www.astropy.org/team.html).

## Development procedure

The initial developers of PlasmaPy will create a flexible 
development roadmap that outlines and prioritizes subpackages to be
developed.  The developers will survey existing open source Python
packages in plasma physics.  Priority will be given to determining how
data will be stored and structured.  Developers will break up into
small groups to work on different subpackages.  These small groups
will communicate regularly and work towards interoperability and
common coding practices.

Because Python is new to many plasma physicists, community engagement
is vital.  The CEC will arrange occasional informal trainings early in
the project that are director towards the initial developers.

New code and edits should be submitted as a pull request to the
development branch of the PlasmaPy repository on GitHub.  The pull
request will undergo a code review by the subpackage maintainers
and/or the CC, who will provide suggestions on how the contributor may
update the pull request.  Subpackage maintainers will generally be
responsible for deciding on pull requests with minor changes, while
pull requests with major changes should be decided jointly by the
subpackage maintainers and the CC.  The CC and CEC will develop a
friendly guide on how users may contribute new code to PlasmaPy.

New code should conform to the [PEP 8 style guide for Python
code](https://www.python.org/dev/peps/pep-0008/) and the established
coding style within PlasmaPy.  New code should be submitted with
documentation and tests.  Documentation should be written primarily in
docstrings and follow the [numpydoc documentation style
guide](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt).
Every new module, class and function should have an appropriate
docstring.  The documentation should describe the interface and the
purpose for the method, but generally not the implementation.  The
code itself should be readable enough to be able to explain how it
works.  Documentation should be updated when the code is edited.  The
tests should cover new functionality (especially methods with complex
logic), but the tests should also be readable and easy to maintain.
Existing tests should be updated when necessary (e.g., during
the initial development of a new feature when the API is not yet
stable), but with caution since this may imply loss of backwards
compatibility.

Members of the PlasmaPy community may submit PlasmaPy Enhancement
Proposals (PLEPs) to suggest changes such as major reorganization of a
subpackage, creation of a new subpackage, non-backwards compatible
changes to a stable package, or significant changes to policies and
procedures related to the organization of this project.  The issues
list on GitHub will generally be more appropriate for changes that do
not require community discussion.  The CC shall maintain a GitHub
repository of PLEPs.  PLEPs will be made openly available for
community discussion and transparency for a period of at least four
weeks, during which time the proposal may be updated and revised by
the proposers.  The CC shall approve or decline these proposals after
seeking community input.  The rationale behind the decision and a
summary of the community discussion shall be recorded along with the
PLEP.

## Programming guidelines

### Choice of languages

PlasmaPy shall be written using Python 3.  PlasmaPy shall initially
guarantee compatibility with Python 3.6 and above.  Python 3 is
continually growing, so we will proceed on the general principle that
future updates to PlasmaPy remain compatible with releases of Python
that are up to two years old.  Python 2.7 and below will not be
supported as these versions will no longer be updated past 2020.  The
core package will initially be written solely in Python.

Code readability is more important than optimization, except when
performance is critical.  Code should be optimized only after getting
it to work, and primarily for where there is a performance bottleneck.
Performance-critical parts of the core package will preferably be
written using Cython or Numba to achieve compiled speeds while
maintaining the significant advantages of using a high level language.

### Dependencies

Dependencies have the advantage of providing capabilities that will
enhance PlasmaPy and speed up its development, but the disadvantage
that they can make manual installation more difficult and potentially
frustrating.  Package managers such as Anaconda and Homebrew greatly
simplify installation of Python packages, but there will be situations
where manual installation is necessary (e.g., on some supercomputers
without package managers).  The core package should be able to be
imported using a minimal number of packages (e.g., NumPy, SciPy, and
matplotlib) without getting an import error.  Additional packages may
be included as dependencies of the core package if there is a strong
need for it, and if these packages are easily installed with currently
available package managers.  Subpackages may use additional
dependencies when appropriate.

### Affiliated packages

We will follow the practice of Astropy by having a core package and
affiliated packages.  The core package will contain common tools and
base functionality that most plasma physicists will need.  The
affiliated packages contained in separate repositories will include
more specialized functionality that is needed for subfields of plasma
physics.  This approach will reduce the likelihood of scope creep for
the core package while maintaining avenues for broader development.

### Units

Multiple sets of units are used by plasma physicists.  There exist
some peculiarities with how units are used within plasma physics, such
as how an electron volt is typically used as a measurement of
temperature.  Code will be most readable and maintainable if written
assuming a particular set of units, but there should be enough
flexibility for people in different subfields to choose their
preferred set of units. As the generally most common accepted 
international standard, SI base units will be utilized.  We will use an
existing Python module (e.g., astropy.units or pint) to assign units
to variables and allow straightforward conversion between different
systems of units.
