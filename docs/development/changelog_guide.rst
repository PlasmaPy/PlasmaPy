***************
Changelog Guide
***************

A changelog tells users and contributors what notable changes have been
made between each release.

Creating a changelog entry
==========================

Each changelog entry should be created after a pull request to the main
branch has been made. If the change is sufficiently minor (e.g., a typo
fix or a tweak to the GitHub Actions configuration), then a package
maintainer will add the "No changelog entry needed" label.

#. In the :file:`changelog/` directory, create a new file entitled,
   :file:`⟨number⟩.⟨type⟩.rst`, where ``⟨number⟩`` is replaced with the pull
   request number and ``⟨type⟩`` is replaced with one of the following
   changelog types:

   * ``breaking``: A change which requires users to change code and is not
     backwards compatible. (Not to be used for removal of deprecated features.)
   * ``bugfix``: Fixes a bug.
   * ``doc``: Documentation addition or improvement, like rewording a
     section or adding missing docs.
   * ``feature``: New user facing features and any new behavior.
   * ``removal``: Feature deprecation and/or feature removal.
   * ``trivial``: A change which has no user-facing effect or is tiny.

   For example, pull request
   `#1198 <https://github.com/PlasmaPy/PlasmaPy/pull/1198>`__
   to PlasmaPy includes an update to the documentation, so the file to be
   created will be entitled :file:`changelog/1198.doc.rst`.

   More than one changelog file
   If a pull request includes the removal of a deprecated feature and an
   update to the documentation, then

#. Open that file and write a short description of the changes that
   were made, using the past tense.

   .. code-block:: rst

      Created a page in the contributor guide that describes how to add
      changelog entries.



   Changelog entries are written using reST_.


For this example, we will assume that we are writing a changelog entry

For this example, we will assume that we are writing a changelog entry
for pull request
to PlasmaPy. Change ``1198`` number to the number of the pull request you
created

Create a news fragment
~~~~~~~~~~~~~~~~~~~~~~



1. In the `changelog/` directory, create a file named

Create a file ``⟨pull request number⟩.⟨pull request type⟩.rst``, where the



To create a changelog entry



Changelog entries should be made for any



Changelog entries will be read by users,

Changelog entries will be read by users and developers of other
packages, so please write each entry to be understandable.

If a change is supplanted by another change, both changelog entries
should be included in the change log. A reference to the first package
may be

***********************


The change log should be


This directory contains "news fragments" which are short files that contain a
small **ReST**-formatted text that will be added to the next ``CHANGELOG``.

The ``CHANGELOG`` will be read by users, so this description should be aimed at
PlasmaPy users instead of describing internal changes which are only relevant
to the developers.

Make sure to use full sentences with correct case and punctuation, for example:

    Add another Braginskii transport coefficient to `plasmapy.transport`.

Please try to use Sphinx intersphinx using backticks.

Each file should be named like ``<PULL REQUEST>.<TYPE>.rst``, where ``<PULL
REQUEST>`` is a pull request number and ``<TYPE>`` is one of:

So for example: ``123.feature.rst``, ``456.bugfix.rst``.

If you are unsure what pull request type to use, don't hesitate to ask in your
PR.

Note that the ``towncrier`` tool will automatically reflow your text, so it
will work best if you stick to a single paragraph, but multiple sentences and
links are OK and encouraged.  You can install ``towncrier`` and then run
``towncrier --draft`` if you want to get a preview of how your change will look
in the final release notes.


Changelog guidelines
====================

* Use the past tense to describe
