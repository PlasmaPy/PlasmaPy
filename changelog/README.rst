***************
Changelog Guide
***************

.. This directory contains "news fragments" which are files that contain
   a short description of the changes that will be added to the
   changelog for the next release.

.. The rendered version of this document is in PlasmaPy's online
   documentation at:
   https://docs.plasmapy.org/en/latest/contributing/changelog_guide.html

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

Introduction
============

A changelog tells users and contributors what notable changes have been
made between each release. Pull requests to PlasmaPy need changelog
entries before they can be merged, except when the changes are very
minor. PlasmaPy uses |towncrier| to convert the changelog entries into
the full changelog. Some example changelog entries are:

.. code-block:: rst

   Added a page in the contributor guide that describes how to add
   changelog entries.

   The ``oldname`` argument to `plasmapy.subpackage.module.function` has
   been deprecated and will be removed in a future release. Use
   ``newname`` instead.

.. _add-changelog:

Adding a changelog entry
========================

Please follow these steps to add a changelog entry after submitting a
pull request to PlasmaPy's ``main`` branch.

#. In the :file:`changelog` directory, create a new file entitled
   :file:`⟨number⟩.⟨type⟩.rst`, where ``⟨number⟩`` is replaced with the
   pull request number and ``⟨type⟩`` is replaced with one of the
   following changelog types:

   * ``breaking``: For backwards incompatible changes that would require
     users to change code. Not to be used for removal of deprecated
     features.
   * ``bugfix``: For changes that fix bugs or problems with the code.
   * ``doc``: For changes to the documentation.
   * ``feature``: For new user-facing features and any new behavior.
   * ``internal``: For refactoring of the internal mechanics of the code
     and tests in ways that do not change the application user interface.
   * ``removal``: For feature deprecation and planned removals.
   * ``trivial``: For minor changes that do not change the application
     programming interface.

   Pull request `#1198 <https://github.com/PlasmaPy/PlasmaPy/pull/1198>`__
   includes an update to the documentation, so the file should be named
   :file:`1198.doc.rst`. If you are unsure of which changelog type to
   use, please feel free to ask in your pull request.

   .. note::

      A ``doc`` changelog entry is not necessary if there is a
      corresponding ``feature`` changelog entry.

   .. tip::

      When a pull request includes multiple changes, use a separate
      changelog file for each distinct change.

      If the changes are in multiple categories, include a separate
      changelog file for each category. For example, pull request
      `#1208 <https://github.com/PlasmaPy/PlasmaPy/pull/1208>`__
      included both a breaking change and a new feature, and thus needed
      both :file:`1208.breaking.rst` and :file:`1208.feature.rst`.

      For multiple changes in a single category, use filenames like
      :file:`1208.trivial.1.rst` and :file:`1208.trivial.2.rst`.

#. Open that file and write a short description of the changes that were
   made. As an example, :file:`1198.doc.rst` might include:

   .. code-block:: rst

      Added a page in the contributor guide that describes how to add
      changelog entries.

#. Commit the file and push the change to branch associated with the
   pull request on GitHub.

Changelog guidelines
====================

* Changelog entries will be read by users and developers of PlasmaPy and
  packages that depend on it, so please write each entry to be
  understandable to someone with limited familiarity of the package.

* Changelog entries are not required for changes that are sufficiently
  minor, such as typo fixes or fixed hyperlinks. When this is the case,
  a package maintainer will add the :guilabel:`no changelog entry needed`
  label to the pull request.

* Use the past tense to describe the change, and the present tense to
  describe how the functionality currently works.

* A changelog entry may include multiple sentences to describe important
  context and consequences of the change. Because |towncrier|
  automatically reflows text, keep entries to a single paragraph.

* Use |intersphinx| links to refer to objects within PlasmaPy, and
  include the full namespace. For example, use
  ```~plasmapy.particles.particle_class.Particle``` to refer to
  |Particle|. The tilde is included to hide all but the name of the
  object.

* Show the full former namespace for objects that have been removed or
  moved, and use double back ticks so that the name is rendered as code
  without attempting to create a link.

  .. code-block:: rst

     Removed the ``plasmapy.physics`` subpackage. The functionality from
     that subpackage is now in `plasmapy.formulary`.

* Substitutions as defined in :file:`docs/_global_substitutions.py` may
  be used in changelog entries.

* The pull request number does not need to be included inside the
  changelog entry because it will be added automatically when the
  individual entries are converted into the full changelog.

* When a changelog entry describes changes to functionality, it is not
  necessary to mention the corresponding changes to the tests.

* If a change is supplanted by another change during the release cycle,
  keep the files for both changelog entries. When the change is
  significant, mention in the earlier entry that the change was
  superseded or reverted and include a link to the appropriate pull
  request.

.. _fixing-obsolete-rest-links:

.. tip::

   When removing or moving an object, |reStructuredText| links that
   follow the original namespace will break, causing the documentation
   build to fail.

   Text in single back ticks is used to link to code objects, while text
   in double back ticks is treated as an `inline literal`_. To remedy
   this problem in old changelog entries, change the broken link into an
   inline literal by surrounding it with double back ticks instead.
   Remove the tilde if present. For example,
   ```~plasmapy.subpackage.module.function``` should be changed to:

   .. code-block:: rst

      ``plasmapy.subpackage.module.function``

   Outside of the changelog, the namespace should be corrected rather
   than changed into an inline literal.

Building the changelog
======================

During the release cycle, |towncrier| is used to build the changelog. To
install |towncrier| and the other packages needed to develop PlasmaPy,
go to the top-level directory of your local clone of PlasmaPy and run:

.. code-block:: shell

   pip install -e .[docs]

Configuration files for |towncrier| are in :file:`pyproject.toml`.

To run |towncrier|, enter the top-level directory of PlasmaPy's
repository. To print out a preview of the changelog, run:

.. code-block:: shell

   towncrier --draft

To convert the changelog entries into a changelog prior to the
``v2024.5.0`` release, run:

.. code-block:: shell

   towncrier build --version=v2024.5.0

This command will ask for permission to delete changelog entries, and
then create :file:`CHANGELOG.rst` in the top-level directory.

.. tip::

   |towncrier| can be used to create a new changelog entry and open it
   for editing using a command like:

   .. code-block:: shell

      towncrier create --edit ⟨number⟩.⟨type⟩.rst

   Replace ``⟨number⟩`` with the pull request number and ``⟨type⟩``
   with the one of the changelog types as described above.

.. _inline literal: https://docutils.sourceforge.io/docs/user/rst/quickref.html#inline-markup
