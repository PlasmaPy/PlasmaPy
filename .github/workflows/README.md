# GitHub workflows

The [`.github/workflows`](.) directory contains [YAML] files that
describe the [GitHub Actions] workflows used when developing PlasmaPy.
Several of the workflows invoke [Nox] sessions that are defined in the
top-level [`noxfile.py`](../../noxfile.py).

Continuous integration (CI) workflows include:

- [`ci.yml`](./ci.yml) — perform CI checks during pull requests (PRs)
- [`weekly.yml`](./weekly.yml) — run weekly tests
- [`check-author-included.yml`](./check-author-included.yml) — verify
  that the author of a PR is included in the top-level
  [`CITATION.cff`](../../CITATION.cff) metadata file

Workflows associated with the release process include:

- [`create-release-issue.yml`](./create-release-issue.yml) — create an
  issue containing the release checklist (triggered manually)
- [`mint-release.yml`](./mint-release.yml) — preparing the repository
  for a release (triggered manually)
- [`publish-to-pypi.yml`](./publish-to-pypi.yml) — perform the official
  release to the Python Package Index (triggered by performing a release
  on GitHub)

[github actions]: https://docs.github.com/en/actions
[nox]: https://nox.thea.codes
[yaml]: https://en.wikipedia.org/wiki/YAML
