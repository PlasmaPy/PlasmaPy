[citation page]: https://docs.plasmapy.org/en/stable/about/citation.html[community meeting]: https://www.plasmapy.org/meetings/weekly
[conda-forge feedstock]: https://github.com/conda-forge/plasmapy-feedstock
[Create an issue for the release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/create-release-issue.yml
[Create a release on GitHub]: https://github.com/PlasmaPy/PlasmaPy/releases/new
[`docs/about/citation.rst`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/about/citation.rst
[`docs/conf.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/conf.py
[Draft a new release]: https://github.com/PlasmaPy/PlasmaPy/releases/new
[.mailmap]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.mailmap
[mint a release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/mint-release.yml
[publish to PyPI]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.github/workflows/publish-to-pypi.yml
[release checklist]: https://github.com/PlasmaPy/PlasmaPy/tree/main/.github/markdown/release-checklist.md
[tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/tests.yml
[Update pinned requirements]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/update-pinned-reqs.yml
[weekly tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly-tests.yml
[Zenodo]: https://zenodo.org/me/uploads

This issue contains the procedure for releasing a new version of PlasmaPy.

## Ahead of release

### Planning

 - [x] [Create an issue for the release].
 - [ ] Discuss the release timeline at a [community meeting].
 - [ ] Update milestones to prioritize and postpone issues and PRs.

<!-- We have had less need of a feature freeze as the package has become more mature, but we may wish to add this back in the future.
 - [ ] About three weeks before a minor or major release, announce that a feature freeze will occur one week before the anticipated release date. Only pull requests with a limited scope that do not significantly change functionality should be merged during the feature freeze.
 - [ ] Begin a code freeze about two weekdays before a release. Only bugfixes and pull requests that are directly related to the release should be merged during the code freeze.
-->

### Code quality updates

> [!IMPORTANT]
> Split up these changes into multiple small PRs rather than a single monolithic PR.

 - [ ] Revise changelog entries to make sure that they are understandable, necessary, and correctly categorized.
   - Add the `no changelog entry needed` label to skip doing changelog checks.
 - [ ] Build the docs using `make linkcheck` or via [weekly tests], and update any broken links. Use `linkcheck_allowed_redirects` in [`docs/conf.py`] to allow redirects (e.g., from `doi.org`), as necessary.
 - [ ] Open a pull request to re-execute pre-executed notebooks, such as those for charged particle radiography <!-- automate this step? -->
 - [ ] Update [`.mailmap`] with new contributors <!-- delete this file? -->

## Day of release

## Make sure that all tests are passing on `main`

 - [ ] [Update pinned requirements] and perform `pre-commit autoupdate`.
 - [ ] Run the [tests]. [![CI](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/tests.yml)
 - [ ] Run the [weekly tests]. [![weekly tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly-tests.yml/badge.svg?branch=main)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly-tests.yml)
 - [ ] Enjoy life for 15 minutes.
 - [ ] Fix any failures and repeat.

## Perform the release

 - [ ] Begin an upload to [Zenodo] for the new release using the `team@plasmapy.org` login, and reserve a DOI.
 - [ ] Run the GitHub Action to [mint a release]. Specify the version (i.e., `2024.1.0` or `2024.1.0rc1` for a release candidate) and copy/paste the reserved DOI from Zenodo.  This action will update the DOI, build the changelog, and tag the release.
 - [ ] [Create a release on GitHub].
   - Choose the newly created tag (e.g., `v2024.1.0`) and use it as the title.
   - Automatically generate release notes.
   - Create a discussion for the release.
   - For release candidates (e.g., `v2024.1.0rc1`), specify it as a pre-release.
 - [ ] Create a pull request to merge the release branch back into `main`. <!-- Automate this step? -->
 - [ ] Download a `.tar.gz` file of the release from the [GitHub releases page] and upload it to [Zenodo].
   - [ ] Update the author list with new authors from the automatically generated release notes.
   - [ ] Update the bibliography, and publish the release to Zenodo.
 - [ ] After fixing any broken tests, merge the pull request from the release branch back into `main`.
 - [ ] Fix any problems with the automated pull request to [conda-forge feedstock], if necessary. This step should be automatic, but may take a while.

## Update documentation

[stable documentation build]: https://docs.plasmapy.org/en/stable
[`stable`]: https://github.com/PlasmaPy/PlasmaPy/tree/stable
[Read the Docs]: https://readthedocs.org/projects/plasmapy
[versions page on RTD]: https://readthedocs.org/projects/plasmapy/versions/

<!-- This section of the checklist may need revision. -->

> [!NOTE]
> The [stable documentation build] on [Read the Docs] (RTD) will point to the [`stable`] branch on GitHub, if it exists. Otherwise, it will point to the most recent release on GitHub. The [`stable`] branch is needed if the documentation build for the release fails or if we wish to make any changes to documentation between releases.

 - [ ] Delete the [`stable`] branch on GitHub, if it exists.
 - [ ] Activate the new release on the [versions page on RTD]
 - [ ] Check that the documentation builds correctly for the release branch, and that the [citation page] is up-to-date.
 - [ ] If the documentation build fails, create a new [`stable`] branch from the release branch (e.g., `2023.10.x`) and fix any problems with the documentation build.

## Test the release

 - [ ] After activating a new virtual or conda environment,  make sure that the released version installs correctly.
 - [ ] Open Python and run `import plasmapy`, `dir(plasmapy)`, and `plasmapy.__version__`.
 - [ ] Run `plasma-calculator` from the command line.
 - [ ] Verify that the new version can be installed with conda.

## After the release

 - [ ] Announce the release
 - [ ] Discuss the release at the [community meeting]
 - [ ] Update the [release checklist] and raise issues about the release process could be improved
