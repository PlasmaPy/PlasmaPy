# Pinned requirements files

[`.github/workflows/update-pinned-reqs.yml`]: ../.github/workflows/update-pinned-reqs.yml

This directory contains pinned requirements files for use by the `tox`
environments that are used in continuous integration testing. Pinning
requirements files makes sure that continuous integration tests are
all run on a set of dependencies that are known to work.

These requirements files are regenerated weekly via an automated pull
request created through the GitHub workflow defined in
[`.github/workflows/update-pinned-reqs.yml`]. These pull requests help
us find out when an upstream dependency causes new test failures, and
usually prevent spontaneous test failures.

## Regenerating requirements

[workflow to update pinned requirements]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/update-pinned-reqs.yml
[#2597]: https://github.com/PlasmaPy/PlasmaPy/pull/2597

Package maintainers may trigger the [workflow to update pinned
requirements] through going to that page, selecting _Run workflow_ and
choosing the `main` branch. This workflow will create a pull request
(e.g., [#2597]) that updates the files in this directory.

To regenerate the requirements locally, run:

```console
tox -e requirements
```
