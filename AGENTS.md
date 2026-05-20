# AI Agent Instructions

This file provides context, rules, and guidelines for AI coding assistants working on this repository.

## Development environment

- This project uses `uv` to manage the development environment.
- The development environment can be created with `uv venv` and synced with `uv sync`.
- If there is a merge conflict in `uv.lock`, first try running `uv lock`. If that doesn't work, run `uvx nox --session lock`.

## Contribution guidelines

- Every commit made by an AI coding assistant should indicate that it was made by an AI coding assistant, such as by starting each commit message with the name of the AI coding assistant in square brackets (e.g., `[Gemini]`, `[ChatGPT]`, or `[CoPilot]`).
- Do not fix issues labeled with "good first issue" because they are reserved to provide human contributors with a chance to learn how to manually make a contribution.

## Contributor guide

- PlasmaPy's contributor guide is located in `docs/contributing/`.
- Coding guidelines are in `docs/contributing/coding_guide.rst`.
- Documentation guidelines are in `docs/contributing/doc_guide.rst`.
- Testing guidelines are in `docs/contributing/testing_guide.rst`.
- Changelog instructions are in `changelog/README.rst`.

## Style

- Follow the PEP 8 style guide.
- Lint files with `pre-commit`, which uses `ruff`.
- Follow the style conventions in `.editorconfig`.

## Documentation

- The documentation is written in reStructuredText, built with Sphinx, and hosted on Read the Docs.
- The documentation can be built with `nox --session docs`.
- Use the numpydoc style for docstrings.

## Changelog

- Except for trivial changes, each pull request should have a changelog entry in the `changelog/` directory using the guidelines in `changelog/README.md`.

## Testing instructions

- Tests are located in the `tests/` directory.
- Find the CI plan in `.github/workflows/`, with Nox sessions defined in `noxfile.py`. The list of Nox sessions can be found by running `nox --list`.
- The full test suite is run with `uvx nox --session 'tests-3.14(all)'` in the top-level directory, which invokes the `tests` session defined in `noxfile.py`.
- Use pytest to run small numbers of tests when trying to fix failing tests. For example, use `pytest tests/particles/test_atomic.py::test_half_life` to run the test named `test_half_life` in `tests/particles/test_atomic.py`, or use `pytest tests/formulary/test_lengths.py` to run all tests in `tests/formulary/test_lengths.py`.
- Fix any test errors until all tests are passing.

## Static type checking

- Static type checking is performed with `nox --session ty`.
- Fix static type checking errors when possible, but use in-line ignore comments if necessary.
