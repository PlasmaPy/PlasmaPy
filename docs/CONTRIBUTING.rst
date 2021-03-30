.. _contributing-to-plasmapy:

How to Contribute
=================

There are numerous ways to contribute to PlasmaPy, including by
providing code and documentation, suggesting and discussing ideas,
submitting issues and bug reports, and engaging the broader plasma
physics community.

.. _impostor-syndrome-disclaimer:

Impostor syndrome disclaimer [1]_
---------------------------------

We want your help. No, really.

There may be a little voice inside your head that is telling you that
you're not ready to be an open source contributor; that your skills
aren't nearly good enough to contribute. What could you possibly offer a
project like this one?

We assure you - the little voice in your head is wrong. If you can write
code at all, you can contribute code to open source. Contributing to
open source projects is a fantastic way to advance one's coding skills.
Writing perfect code isn't the measure of a good developer (that would
disqualify all of us!); it's trying to create something, making
mistakes, and learning from those mistakes. That's how we all improve,
and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either.
You can help out by writing documentation, tests, or even giving
feedback about the project (and yes - that includes giving feedback
about the contribution process). Some of these contributions may be the
most valuable to the project as a whole, because you're coming to the
project with fresh eyes, so you can see the errors and assumptions that
seasoned contributors have glossed over.

Contributing code or documentation to PlasmaPy
----------------------------------------------

If you see something you’d like to work on amongst our
`issues <https://github.com/PlasmaPy/PlasmaPy/issues>`_, start hacking
away on that! However, please announce your intent first in the relevant
issue to make sure there is no work duplication.

Please note that PlasmaPy has a :ref:`plasmapy-code-of-conduct`.

Issues marked by the community as *help wanted* mean just that - either
they’re good contributions for outsiders or there’s an issue in the
ongoing work that requires a second opinion. Please consider these
first!

Work on PlasmaPy is done via GitHub, so you’ll need a `(free)
account <https://github.com/join?source=header-home>`_. If you are new
to `git <https://git-scm.com/>`_, helpful resources include
documentation on `git
basics <https://git-scm.com/book/en/v2/Getting-Started-Git-Basics>`_
and an `interactive git
tutorial <https://try.github.io/levels/1/challenges/1>`_. You must also
`install
git <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_
locally on your computer. We highly recommend getting reasonably
familiar with git by going through these tutorials or a `Software
Carpentry <https://software-carpentry.org/>`_ workshop prior to making
code contributions. Do note that you can usually find help in the
`PlasmaPy Matrix
chatroom <https://app.element.io/#/room/#plasmapy:matrix.org>`_.

For actual guidelines for working on PlasmaPy, please see our
:ref:`plasmapy-development-guide`.

Towncrier changelog entries
---------------------------

Every pull request should include a changelog entry. Please see
`changelog/README.rst` for instructions.

To summarize, put a file like ``<PULL REQUEST>.<TYPE>.rst``, where ``<PULL
REQUEST>`` is a pull request number, and ``<TYPE>`` is one of ``breaking``,
``feature``, ``bugfix``, ``doc``, ``removal``, ``trivial``. If unsure, ask
a maintainer.

Footnotes
^^^^^^^^^

.. [1] The `imposter syndrome disclaimer
       <https://github.com/adriennefriend/imposter-syndrome-disclaimer>`_
       was originally written by `Adrienne Lowe
       <https://github.com/adriennefriend>`_ for a `PyCon talk
       <https://www.youtube.com/watch?v=6Uj746j9Heo>`_.  It was adapted
       in the README files for
       `MetPy <https://github.com/Unidata/MetPy>`_ and `yt
       <https://github.com/yt-project/yt>`_, and was then adapted by
       PlasmaPy.
