# Pinned requirements files

[`.github/workflows/update-pinned-reqs.yml`]: ../.github/workflows/update-pinned-reqs.yml
[Nox]: https://nox.thea.codes/en/stable/

This directory contains pinned requirements files for use by the [Nox]
sessions that are used in continuous integration testing. Pinning
requirements files ensures that continuous integration tests are
consistently run on a set of dependencies that are known to work.

These requirements files are regenerated via weekly automated pull
requests created through the GitHub workflow defined in
[`.github/workflows/update-pinned-reqs.yml`]. These pull requests help
us find out when an upstream dependency causes new test failures, and
usually prevent spontaneous test failures.

## Regenerating requirements

[workflow to update pinned requirements]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/update-pinned-reqs.yml
[#2597]: https://github.com/PlasmaPy/PlasmaPy/pull/2597
[`pyproject.toml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/pyproject.toml

When the dependencies defined in [`pyproject.toml`] change in a pull
request, it is necessary to regenerate these requirements files. To
regenerate the requirements locally, run:

```console
nox -s requirements
```

Package maintainers may trigger the [workflow to update pinned
requirements] through going to that page, selecting _Run workflow_ and
choosing the `main` branch. This workflow will create a pull request
(e.g., [#2597]) that updates the files in this directory.
