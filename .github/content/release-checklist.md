This issue is created from the [release checklist] for releasing a new version of PlasmaPy.

### Planning the release

- [x] [Create an issue for the release]. 📝
- [ ] Update [milestones] for issues & pull requests (PRs). 🛣️

### Quality assurance checks (low priority)

- [ ] Run the [GitHub workflow for checking hyperlinks], and update
  broken links. Use `linkcheck_allowed_redirects` in [`docs/conf.py`] to
  allow redirects (e.g., from `doi.org`). Update or delete the `alias`
  field for authors in [`CITATION.cff`] who have changed their GitHub
  username. [![linkcheck](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/linkcheck.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/linkcheck.yml)

- [ ] Run `git log --format="%aN <%aE>" | sort -u` in a bash terminal,
  and update [`.mailmap`] if there are any duplicate contributors in the
  output ([gitmailmap documentation]). 📫

- [ ] Build the [upstream docs] and fix failures, as appropriate. These checks build the documentation against unreleased versions of core dependencies. [![upstream docs](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-docs.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-docs.yml)

### Quality assurance checks (high priority)

- [ ] [Upgrade `uv.lock`], and fix new test failures. 🔒

- [ ] Remove any remaining upper limits in the `dependencies` field in [
  `pyproject.toml`][`pyproject.toml`]. ⏩

- [ ] Run `pre-commit autoupdate` followed by `pre-commit run --all-files`. Fix new errors and commit the changes. [![pre-commit.ci status](https://results.pre-commit.ci/badge/github/PlasmaPy/PlasmaPy/main.svg)](https://results.pre-commit.ci/latest/github/PlasmaPy/PlasmaPy/main)

  - Occasionally, certain hooks may need to be manually downgraded after running `pre-commit autoupdate` because of problems with the latest versions of these hooks. Look for comments in [`.pre-commit-config.yaml`] for guidance. 🪝

- [ ] Address issues labeled as [needed for release](https://github.com/PlasmaPy/PlasmaPy/issues?q=state%3Aopen%20label%3A%22needed%20for%20release%22).

- [ ] Run the [upstream tests] and fix failures, as appropriate. These checks run tests against unreleased versions of core dependencies. [![upstream tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-tests.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-tests.yml)

- [ ] Fix any remaining deprecation warnings, including any ignored under `tool.pytest.filterwarnings` in [`pyproject.toml`].

### Perform the release

- [ ] Begin an upload to [Zenodo] for a new version of [this record], using the `team@plasmapy.org` login. Reserve a DOI. 🔢

- [ ] Run the GitHub workflow to [prepare a release], specifying the version (i.e., `2026.1.0`) and copying the reserved DOI from Zenodo. This workflow will create a PR that builds the changelog and updates package metadata. 🤖 [![prepare release PR](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/prepare-release-pr.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/prepare-release-pr.yml)

  - [ ] Revise changelog entries to make sure that they are understandable, necessary, and correctly categorized. 📜
  - [ ] Make sure that all tests are passing in the PR. ✅
  - [ ] Merge the PR. 📦

- [ ] Run the [tests]. [![tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci.yml)

- [ ] Run the [comprehensive tests]. [![comprehensive tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-comprehensive.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-comprehensive.yml) 🔍

- [ ] [Create a release on GitHub]. 🚀

  - [ ] Choose the newly created tag (e.g., `v2026.1.0`), and use it as the title. (The release will be performed from the tag, so it is not necessary to select the branch.) 🏷️
  - [ ] Set the tag for the previous release, and select the option to automatically generate release notes. 📜
  - [ ] Select the option to create a discussion for the release under the _General_ category. 📣
  - [ ] For official releases, choose _Set as the latest release_. For beta releases or release candidates (e.g., `v2026.1.0rc1`), specify it as a pre-release. 🆕
  - [ ] Click on <kbd>Publish release</kbd>, which will create the GitHub release and trigger the GitHub workflow to [publish to PyPI]. [![publish](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/publish-to-pypi.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/publish-to-pypi.yml)
  - [ ] Check the [release history] on PyPI to make sure that the release was successful. [![PyPI version](https://img.shields.io/pypi/v/plasmapy?style=flat&logo=pypi)](https://pypi.org/project/plasmapy/) [![PyPI version](https://img.shields.io/pypi/pyversions/plasmapy?style=flat&logo=python)](https://img.shields.io/pypi/pyversions/plasmapy?style=plastic)

### Upload release to Zenodo

- [ ] Download a `.tar.gz` file of the tagged release from the [list of tagged versions] on GitHub, and upload it to [Zenodo]. 📤
  - [ ] Update the version number and release date in the record. 📅
  - [ ] Update the author list with new authors from the automatically generated release notes or [`CITATION.cff`]. 👥
  - [ ] Update the bibliography. 📖
  - [ ] Publish the record. 🏛️

> [!TIP]
> To compare two files across different tags, use commands like:
>
> ```shell
> git diff v2024.7.0:CITATION.cff v2024.10.0:CITATION.cff
> git diff v2024.7.0:docs/bibliography.bib v2024.10.0:docs/bibliography.bib
> ```

### Documentation

- [ ] Delete the [`stable`] branch on GitHub if it exists. 🗑️

- [ ] Verify that the current and previous release are activated on
  the [versions page on RTD].
  📚 [![Read the Docs Status](https://readthedocs.org/projects/plasmapy/badge/?version=stable)](http://plasmapy.readthedocs.io/en/latest/?badge=stable)

- [ ] Verify that the [citation page] is up-to-date and the DOI link points to the most recent release. 🧾

## Availability on conda-forge

Within a day, an automated PR will be made to PlasmaPy's [conda-forge feedstock], which will require updates if any requirements have changed.

- [ ] Update [`recipe/meta.yaml`] to match [`pyproject.toml`] in the release. 🔧

- [ ] Verify that `python_min` near the beginning of [`recipe/meta.yaml`] is consistent with `requires-python` [`pyproject.toml`]. 🐍

- [ ] Merge the PR to the conda-forge feedstock. 🚀

- [ ] Verify that the new version shows up on conda-forge. 📦 [![Conda version](https://img.shields.io/conda/v/conda-forge/plasmapy?style=flat&logo=anaconda)](https://img.shields.io/conda/v/conda-forge/plasmapy)

> [!TIP]
> If the documentation build fails, fix any problems with the RTD build and create a new [`stable`] branch from the release. The [stable documentation build] will point to the `stable` branch on GitHub if it exists. Otherwise, it will point to the most recent release on GitHub. An alternative would be to fix the problem on `main` and perform a patch release. Because the release process has changed since we needed a `stable` branch, please update the [release checklist] accordingly.

## Test releases

- [ ] Run the [installability tests] to make sure that the new version of PlasmaPy can be installed via pip, uv, and conda-forge. [![installability tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/installability.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/installability.yml)

## After the release

- [ ] Update the [release checklist], if necessary. 📋

- [ ] [Create an issue for the next release], and then [pin the issue]. ⏳

- [ ] Close the issue and celebrate. 🎆

[citation page]: https://docs.plasmapy.org/en/stable/about/citation.html
[comprehensive tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-comprehensive.yml
[conda-forge feedstock]: https://github.com/conda-forge/plasmapy-feedstock
[create a release on github]: https://github.com/PlasmaPy/PlasmaPy/releases/new
[create an issue for the next release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/create-release-issue.yml
[create an issue for the release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/create-release-issue.yml
[github workflow for checking hyperlinks]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/linkcheck.yml
[gitmailmap documentation]: https://git-scm.com/docs/gitmailmap
[installability tests]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.github/workflows/installability.yml
[list of tagged versions]: https://github.com/PlasmaPy/PlasmaPy/tags
[milestones]: https://github.com/PlasmaPy/PlasmaPy/milestones
[pin the issue]: https://docs.github.com/en/issues/tracking-your-work-with-issues/administering-issues/pinning-an-issue-to-your-repository
[prepare a release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/prepare-release-pr.yml
[publish to pypi]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.github/workflows/publish-to-pypi.yml
[release checklist]: https://github.com/PlasmaPy/PlasmaPy/tree/main/.github/content/release-checklist.md
[release history]: https://pypi.org/project/plasmapy/#history
[stable documentation build]: https://docs.plasmapy.org/en/stable
[tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci.yml
[this record]: https://zenodo.org/doi/10.5281/zenodo.6774349
[upgrade `uv.lock`]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upgrade-uv-lock.yml
[upstream tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-upstream.yml
[versions page on rtd]: https://readthedocs.org/projects/plasmapy/versions/
[zenodo]: https://zenodo.org/me/uploads
[`.mailmap`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.mailmap
[`.pre-commit-config.yaml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.pre-commit-config.yaml
[`citation.cff`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/CITATION.cff
[`docs/conf.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/conf.py
[`pyproject.toml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/pyproject.toml
[`recipe/meta.yaml`]: https://github.com/conda-forge/plasmapy-feedstock/blob/main/recipe/meta.yaml
[`stable`]: https://github.com/PlasmaPy/PlasmaPy/tree/stable
