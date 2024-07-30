# Pinned requirements files

[`.github/workflows/update-pinned-reqs.yml`]: ../.github/workflows/update-pinned-reqs.yml
[Nox]: https://nox.thea.codes/en/stable/
[`pyproject.toml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/pyproject.toml
[workflow to update pinned requirements]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/update-pinned-reqs.yml

This directory contains pinned requirements files for use by the [Nox]
sessions that run tests, build documentation, and perform other
checks. Pinning requirements files ensures that continuous integration
(CI) tests are consistently run on a set of dependencies that are
known to work.

These requirements files are regenerated via weekly automated pull
requests (PRs) created through the GitHub workflow defined in
[`.github/workflows/update-pinned-reqs.yml`]. These PRs help us
discover and address new test failures that result from breaking
changes in upstream dependencies. If we did not pin dependencies,
breaking changes would result in PRs having spontaneous CI failures
that are unrelated to the PR and thus difficult to track down.

## Regenerating requirements

### Locally with Nox

When a PR changes the dependencies as defined in [`pyproject.toml`],
it is necessary to regenerate these requirements files. To regenerate
the requirements locally using [Nox], go to the top-level directory of
your clone of PlasmaPy and run

```console
nox -s requirements
```

before committing the changes and pushing them to GitHub.

### On GitHub

Package maintainers may trigger the [workflow to update pinned
requirements] by going to that page, selecting _Run workflow_, and
choosing the `main` branch. This workflow will run the [Nox] session
described above to regenerate the requirements files in this directory
and create a PR to `main` with those changes (e.g., [#2597]).
