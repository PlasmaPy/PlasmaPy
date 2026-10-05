# AI Agent Instructions

This file provides context, rules, and guidelines for AI coding assistants working on this repository.

## Project overview

- `plasmapy` is a Python package that contains software tools for plasma science.
- Source code is located in `src/`.

## Development environment

- Set up the development environment with `uv venv` followed by `uv sync`.
- If there is a merge conflict in `uv.lock`, run `uv lock`. If that doesn't work, run `uvx nox --session lock`.

## Contribution guidelines

- Start every commit message with the name of the coding agent in square brackets (e.g., `[Gemini]`, `[Claude]` or `[Copilot]`).

## Contributor guide

- PlasmaPy's contributor guide is located in `docs/contributing/`.
- Coding guidelines are in `docs/contributing/coding_guide.rst`.
- Documentation guidelines are in `docs/contributing/doc_guide.rst`.
- Testing guidelines are in `docs/contributing/testing_guide.rst`.
- Changelog instructions are in `changelog/README.rst`.

## Style

- Lint files with `uvx pre-commit`.
- Use SI units.

## Units and particles

- Represent physical quantities with `astropy.units` (usually imported as
  `u`; e.g., `5 * u.m`). Keep operations unit-aware and convert with
  `.to(...)`.
- Treat eV as energy. For temperature arguments that accept energy units,
  add a temperature-energy equivalency to that parameter in
  `@validate_quantities`, e.g., `T={"equivalencies": u.temperature_energy()}`.
- Use `Particle` for known species and `CustomParticle` for species with
  specified mass or charge (e.g., an average ion). Use `ParticleList` to
  collect particles.
- Use `@particle_input` when a function should accept particle-like
  values and convert them to particle objects. Annotate each relevant
  parameter with `ParticleLike` or `ParticleListLike`; use `require`,
  `any_of`, or `exclude` for category constraints. Avoid implicit
  default particle assumptions.
- Use `@validate_quantities` to validate units and value constraints and
  convert `Quantity` arguments or configured return values to requested
  units. Prefer type annotations for expected units, and use decorator
  options for additional value constraints or equivalencies.

## Documentation

- The documentation is written in reStructuredText, built with Sphinx, and hosted on Read the Docs.
- The documentation is built with `uvx nox --session docs`.
- Use numpydoc style docstrings.

## Changelog

- When changes are not trivial, include a changelog entry in the `changelog/` directory using the guidelines in `changelog/README.md`.

## Testing instructions

- Tests are located in the `tests/` directory.
- Find the CI plan in `.github/workflows/`, with Nox sessions defined in `noxfile.py`. The list of Nox sessions can be found by running `nox --list`.
- The full test suite is run with `uvx nox --session 'tests-3.14(all)'` in the top-level directory, which invokes the `tests` session defined in `noxfile.py`.
- Use pytest to run small numbers of tests when trying to fix failing tests. For example, use `pytest tests/particles/test_atomic.py::test_half_life` to run the test named `test_half_life` in `tests/particles/test_atomic.py`, or use `pytest tests/formulary/test_lengths.py` to run all tests in `tests/formulary/test_lengths.py`.
- Fix any test errors until all tests are passing.

## Static type checking

- Static type checking is performed with `nox --session ty`.
- Fix static type checking errors when possible, but use in-line ignore comments if necessary.
