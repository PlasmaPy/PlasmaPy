***************
Changelog Guide
***************

A changelog tells users and contributors what notable changes have been
made between each release. Pull requests to PlasmaPy need changelog
entries before they can be merged, except when the changes are very
minor. PlasmaPy uses towncrier_ to convert the changelog entries into
the full changelog.

Creating a changelog entry
==========================

Please follow these steps to add a changelog entry after submitting a
pull request to PlasmaPy's ``main`` branch.

#. In the :file:`changelog` directory, create a new file entitled
   :file:`⟨number⟩.⟨type⟩.rst`, where ``⟨number⟩`` is replaced with the
   pull request number and ``⟨type⟩`` is replaced with one of the
   following changelog types:

   * ``breaking``: A change which requires users to change code and is
     not backwards compatible. (Not to be used for removal of deprecated
     features.)
   * ``bugfix``: Fixes a bug.
   * ``doc``: Documentation addition or improvement, like rewording a
     section or adding missing docs.
   * ``feature``: New user facing features and any new behavior.
   * ``removal``: Feature deprecation and/or feature removal.
   * ``trivial``: A change which has no user-facing effect or is tiny.

   Pull request `#1198 <https://github.com/PlasmaPy/PlasmaPy/pull/1198>`__
   includes an update to the documentation, so the file should be named
   :file:`1198.doc.rst`. If you are unsure of which changelog type to
   use, please feel free to ask in your pull request.

   .. tip::

      When a pull request includes multiple changes, use a separate
      changelog entry for each change to improve readability.

      When there are changes in multiple categories, include a separate
      changelog file for each category. For example, pull request
      `#1208 <https://github.com/PlasmaPy/PlasmaPy/pull/1208>`__
      included both a breaking change and a new feature, and thus needed
      both :file:`1208.breaking.rst` and :file:`1208.feature.rst`.

      When there are multiple changes in a single category, use
      filenames like :file:`1206.trivial.1.rst` and
      :file:`1206.trivial.2.rst`.

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
  minor, such as typo fixes. When this is the case, a package maintainer
  will add the *No changelog entry needed* label to the pull request.

* Use the past tense to describe the change, and the present tense to
  describe how the functionality currently works.

* A changelog entry may include multiple sentences to describe important
  context and consequences of the change. Because towncrier_
  automatically reflows text, keep entries to a single paragraph.

* Use intersphinx_ links to refer to objects within PlasmaPy, and
  include the full namespace. For example, use
  ```~plasmapy.particles.particle_class.Particle``` to refer to
  |Particle|. The tilde is used to hide all but the name of the object.

* Show the full namespace for objects that have been removed or moved,
  and use double back ticks so that the name is rendered as code without
  attempting to create a link.

  .. code-block:: rst

     Removed the ``plasmapy.physics`` subpackage. The functionality in
     that subpackage is now in `plasmapy.formulary`.

* Substitutions as defined in :file:`common_links.rst` may be used in
  changelog entries.

* The pull request number does not need to be included inside the
  changelog entry because it will be added automatically when the
  individual entries are converted into the full changelog.

* If a change is supplanted by another change during the release cycle,
  keep the files for both changelog entries. When the change is
  significant, mention in the earlier entry that the change was
  superseded or reverted and include a link to the appropriate pull
  request.

Building the changelog
======================

To install towncrier_ and the other packages needed to develop PlasmaPy,
go to the top-level directory of your local clone of PlasmaPy and run:

.. code-block:: shell

   pip install -r requirements.txt

To run towncrier_, enter the top-level directory of PlasmaPy's
repository. To print out a preview of the changelog, run:

.. code-block:: shell

   towncrier --draft

To convert the changelog entries into a changelog prior to the 0.7.0
release, run:

.. code-block:: shell

   towncrier --version v0.7.0

This will create :file:`CHANGELOG.rst` in the top-level directory. The
full steps to update the changelog are described in the
:ref:`Release Guide`.

Configuration files for towncrier_ are in :file:`pyproject.toml`.
