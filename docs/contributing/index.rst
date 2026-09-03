.. _contributor guide:

*****************
Contributor Guide
*****************

Thank you for your interest in contributing to PlasmaPy! |:sparkles:|
The future of the project depends on people like you, so we deeply
appreciate it! |:seedling:|

Please feel free to reach out to us during one of PlasmaPy's
|community meetings|.  The PlasmaPy community abides by the
:ref:`Contributor Covenant Code of Conduct <plasmapy-code-of-conduct>`.

If you are becoming a first-time contributor, we recommend starting with:

.. toctree::
   :caption: Getting started
   :maxdepth: 1

   many_ways
   getting_ready
   workflow

After gaining some familiarity with the code contribution workflow, we
suggest checking out:

.. toctree::
   :caption: Contributing to PlasmaPy
   :maxdepth: 2

   coding_guide
   testing_guide
   doc_guide
   changelog_guide
   pre-commit

The contributions are made to |PlasmaPy's GitHub repository|. Thank you
again!

.. _ai-contribution-policy:

PlasmaPy AI Contribution Policy
=================================

Overview
--------

PlasmaPy welcomes contributions from all the community members. As AI-assisted development
is becoming increasingly prevalent, we are adopting this policy to ensure contributions
maintain transparency and alignment with the project's goals.

Core Principles
---------------

1. **Transparency**: All AI usage must be clearly disclosed
2. **Human Responsibility**: Contributors retain full responsibility for their submissions
3. **Code Quality**: AI-generated code must be thoroughly reviewed and understood by the contributor
4. **Maintainer Respect**: We value our maintainers' time and require contributors to do due diligence

Policy Details
--------------

What is Allowed
~~~~~~~~~~~~~~~

- **AI-assisted code** with full human review and testing
- **AI for grammar/documentation improvement** (especially for non-native English speakers)
- **AI for code refactoring or optimization** (with human validation)

What is Not Allowed
-------------------

- **Fully automated contributions**: PRs that are generated entirely by AI without any sort of human involvement
- **Unreviewed AI code**: Submitting AI-generated code without manually understanding it
- **Bulk AI-pasted text**: Copying AI-generated text directly into PR descriptions, issue comments, or documentation without review

Requirements for AI-Assisted Contributions
------------------------------------------

Mandatory Disclosures
~~~~~~~~~~~~~~~~~~~~~

You **must** include in your PR description:

1. **Which AI tools were used** (e.g., "Claude Code for code structure", "ChatGPT for docstring improvement")
2. **Which parts of the contribution were AI-assisted** (be specific: "AI helped with function logic in ``process_data()``, docstrings written manually")
3. **A brief description of how AI was used** (e.g., "Used ChatGPT to generate initial algorithm structure, then rewrote and tested thoroughly")

Testing Requirements
--------------------

- AI-generated code must pass all existing tests
- New tests should be written (can be AI-assisted, but must be manually reviewed)
- Do not rely solely on AI-generated tests to validate AI-generated code

Labeling and Review Process
---------------------------

PR Labels
---------

Contributors are advised to add the **``AI assisted``** label if their PR uses AI assistance.
This helps project maintainers:

- Plan review capacity (AI-assisted PRs may face longer review times)
- Understand how AI is being used in the project
- Prioritize based on complexity of the PR

Review Expectations
~~~~~~~~~~~~~~~~~~~

- **AI-assisted PRs may likely face delayed code reviews**

Questions?
----------

If you have any questions/doubts around this, feel free to reach out to the project maintainers.
