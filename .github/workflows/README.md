# GitHub workflows

The [`.github/workflows`](.) directory contains [YAML] files that describe the [GitHub Actions] workflows used for during PlasmaPy development, including for continuous integration (CI) checks on pull requests (PRs).

## Workflows

### CI

- [`ci.yml`](./ci.yml) — perform standard continuous integration (CI) checks on PRs
- [`ci-comprehensive.yml`](./ci-comprehensive.yml) — run comprehensive tests
- [`upstream-tests.yml`](./upstream-tests.yml) — test against unreleased versions of upstream dependencies to find breaking changes before they make their way into a release
- [`upstream-docs.yml`](./upstream-docs.yml) — build documentation against unreleased versions of upstream dependencies
- [`linkcheck.yml`](./linkcheck.yml) — check the documentation for broken hyperlinks

### Maintenance and triage

- [`upgrade-uv-lock.yml`](./upgrade-uv-lock.yml) — update the locked Python environments used in CI

### Release process

- [`create-release-issue.yml`](./create-release-issue.yml) — trigger manually to create an issue containing the [release checklist](../content/release-checklist.md)
- [`prepare-release-pr.yml`](./prepare-release-pr.yml) — trigger manually to prepare the repository for a release
- [`publish-to-pypi.yml`](./publish-to-pypi.yml) — automatically perform the official release to the Python Package Index (PyPI) after performing a release on GitHub

### Quality assurance

- [`changelog.yml`](./changelog.yml) — verify that a valid changelog entry exists, unless labeled with "no changelog entry needed"
- [`check-author-included.yml`](./check-author-included.yml) — verify that the author of a PR is included in [`CITATION.cff`](../../CITATION.cff)
- [`installability.yml`](./installability.yml) – test package installation from official channels
- [`pyhc-actions.yml`](./pyhc-actions.yml) – check for interoperability with the Python in Heliophysics Community (PyHC) environment

### Triage

- [`comment-on-pr.yml`](./comment-on-pr.yml) — comment on PRs with information on how to contribute
- [`labeler.yml`](./labeler.yml) — add labels to PRs
- [`stale.yml`](./stale.yml) — close issues and PRs that have been inactive for a very long time
- [`unlabel-pr-after-merge.yml`](./unlabel-pr-after-merge.yml) — remove certain labels from PRs after merging

[github actions]: https://docs.github.com/en/actions
[yaml]: https://en.wikipedia.org/wiki/YAML
