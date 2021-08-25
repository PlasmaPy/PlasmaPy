***************
Changelog Guide
***************

A changelog tells users and contributors what notable changes have been
made between each release.

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

      If a pull request includes multiple changes, add a separate
      changelog entry for each change. For example, because pull request
      `#1208 <https://github.com/PlasmaPy/PlasmaPy/pull/1208>`__
      includes a breaking change as well as a new feature, the two files
      to be created will be :file:`1208.breaking.rst` and
      :file:`1208.feature.rst`.

      If multiple changes were made in the same category, use filenames
      like :file:`1206.trivial.1.rst` and :file:`1206.trivial.2.rst`.

#. Open that file and write a short description of the changes that were
   made, using the past tense.

   For example, :file:`1198.doc.rst` might include:

   .. code-block:: rst

      Added a page in the contributor guide that describes how to add
      changelog entries.

   Changelog entries are written using reST_.



#. Commit the file and push the change to branch associated with the
   pull request on GitHub.

Changelog guidelines
====================

* Changelog entries will be read by users and developers of PlasmaPy and
  packages that depend on it, so please write each entry to be
  understandable to someone with limited familiarity of the package.

* Use the past tense when describing the changes that were made.

* Use full sentences with correct case and punctuation.

* Use intersphinx_ links to refer to objects within PlasmaPy, and
  include the full namespace. For example, use
  ``` `~plasmapy.particles.particle_class.Particle` ``` to refer to
  |Particle|. The tilde is used to hide all but the name of the object.

* Show the full namespace for features that have been removed or moved,
  and use double back ticks so that the name is rendered as code without
  attempting to create a link. For example, an entry could mention that
  the ``` ``plasmapy.physics`` ``` subpackage was removed with the contents
  being moved into `plasmapy.formulary`.

  .. code-block:: rst

     Removed the ``plasmapy.physics`` subpackage and move the contents
     into `plasmapy.formulary`.

* Changelog entries are not required for changes that are sufficiently
  minor, such as typo fixes. When this is the case, a package maintainer
  will add the *No changelog entry needed* label to the pull request.

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

PlasmaPy uses towncrier_ to convert the changelog entries (called "news
fragments") into the full changelog.

To install towncrier_ and the other packages needed to develop PlasmaPy,
go to the top-level directory of your local clone of PlasmaPy and run:

.. code-block:: shell

   pip install -r requirements.txt

To print out a preview of the changelog, run

.. code-block:: shell

   towncrier --draft


.. Please try to use Sphinx intersphinx using backticks.

.. Each file should be named like ``<PULL REQUEST>.<TYPE>.rst``, where ``<PULL
REQUEST>`` is a pull request number and ``<TYPE>`` is one of:

.. So for example: ``123.feature.rst``, ``456.bugfix.rst``.

.. If you are unsure what pull request type to use, don't hesitate to ask in your
PR.

Note that the ``towncrier`` tool will automatically reflow your text, so it
will work best if you stick to a single paragraph, but multiple sentences and
links are OK and encouraged.  You can install ``towncrier`` and then run
``towncrier --draft`` if you want to get a preview of how your change will look
in the final release notes.
