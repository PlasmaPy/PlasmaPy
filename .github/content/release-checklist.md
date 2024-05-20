[citation page]: https://docs.plasmapy.org/en/stable/about/citation.html
[`CITATION.cff`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/CITATION.cff
[community meeting]: https://www.plasmapy.org/meetings/weekly
[conda-forge feedstock]: https://github.com/conda-forge/plasmapy-feedstock
[Create an issue for the release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/create-release-issue.yml
[Create a release on GitHub]: https://github.com/PlasmaPy/PlasmaPy/releases/new
[`docs/about/citation.rst`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/about/citation.rst
[`docs/conf.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/conf.py
[documentation]: https://docs.plasmapy/org/en/stable
[Draft a new release]: https://github.com/PlasmaPy/PlasmaPy/releases/new
[GitHub Action for checking hyperlinks]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/linkcheck.yml
[list of tagged versions]: https://github.com/PlasmaPy/PlasmaPy/tags
[`.mailmap`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.mailmap
[milestones]: https://github.com/PlasmaPy/PlasmaPy/milestones
[mint a release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/mint-release.yml
[publish to PyPI]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.github/workflows/publish-to-pypi.yml
[release checklist]: https://github.com/PlasmaPy/PlasmaPy/tree/main/.github/content/release-checklist.md
[tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/tests.yml
[Update pinned requirements]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/update-pinned-reqs.yml
[weekly tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly-tests.yml
[Zenodo]: https://zenodo.org/me/uploads
[stable documentation build]: https://docs.plasmapy.org/en/stable
[`stable`]: https://github.com/PlasmaPy/PlasmaPy/tree/stable
[Read the Docs]: https://readthedocs.org/projects/plasmapy
[versions page on RTD]: https://readthedocs.org/projects/plasmapy/versions/

This issue contains the procedure for releasing a new version of PlasmaPy.

### Planning the release

 - [x] [Create an issue for the release].
 - [ ] Update [milestones] for issues & pull requests (PRs).

<!-- We have had less need of a feature freeze as the package has become more mature, but we may wish to add this back in the future.
 - [ ] About three weeks before a minor or major release, announce that a feature freeze will occur one week before the anticipated release date. Only pull requests with a limited scope that do not significantly change functionality should be merged during the feature freeze.
 - [ ] Begin a code freeze about two weekdays before a release. Only bugfixes and pull requests that are directly related to the release should be merged during the code freeze.
-->

### Code quality updates

 - [ ] Revise changelog entries to make sure that they are understandable, necessary, and correctly categorized. Add the `no changelog entry needed` label to skip doing changelog checks.
 - [ ] Run the [GitHub Action for checking hyperlinks] and update broken links. Use `linkcheck_allowed_redirects` in [`docs/conf.py`] to allow redirects (e.g., from `doi.org`). Update or delete the `alias` field for authors in [`CITATION.cff`] who have changed their GitHub username.
 - [ ] Update [`.mailmap`] with new contributors <!-- delete this file? -->
 - [ ] [Update pinned requirements] and perform `pre-commit autoupdate`.

### Make sure that all tests are passing on `main`

 - [ ] Run the [tests]. [![CI](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/tests.yml)
 - [ ] Run the [weekly tests]. [![weekly tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly-tests.yml/badge.svg?branch=main)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly-tests.yml)
 - [ ] Fix any failures and repeat.

### Perform the release

 - [ ] Begin an upload to [Zenodo] for the new release using the `team@plasmapy.org` login, and reserve a DOI.
 - [ ] Run the GitHub Action to [mint a release]. Specify the version (i.e., `2024.5.0` or `2024.5.0rc1` for a release candidate) and copy/paste the reserved DOI from Zenodo.  This action will update the DOI, build the changelog, and tag the release.
 - [ ] [Create a release on GitHub]. Choose the newly created tag (e.g., `v2024.5.0`) and use it as the title. Select the options to automatically generate release notes and create a discussion for the release. For beta releases or release candidates (e.g., `v2024.5.0rc1`), specify it as a pre-release.
 - [ ] Download a `.tar.gz` file of the tagged release from the [list of tagged versions] and upload it to [Zenodo].
   - [ ] Update the author list with new authors from the automatically generated release notes or [`CITATION.cff`].
   - [ ] Update the bibliography, and publish the release to Zenodo.
 - [ ] Create and merge a pull request from the release branch back into `main`. <!-- Automate pull request creation? Change it into a commit? -->
 - [ ] Fix any problems with the automated pull request to [conda-forge feedstock], if necessary. This step should be automatic, but may take a while.
 - [ ] Update requirements in the [conda-forge feedstock] in `recipe/meta.yaml`, in particular when there is a new version of Python.

### Update documentation

 - [ ] Delete the [`stable`] branch on GitHub, if it exists.
 - [ ] Activate the current and prior release on the [versions page on RTD], if necessary. If the documentation fails to build for a release, activate the corresponding branch (e.g., activate the `v2023.10.x` branch instead of the `v2023.10.0` tag). <!-- true example! -->
 - [ ] Check that the [documentation] builds correctly for the release branch. If the documentation build fails, create a new [`stable`] branch from the release branch (e.g., `2024.5.x`) and fix any problems with the documentation build. The [`stable`] branch is needed if the documentation build for the release fails or if we make any changes to the documentation between releases. The [stable documentation build] will point to the [`stable`] branch on GitHub if it exists. Otherwise, it will point to the most recent release on GitHub.
 - [ ] Verify that the [citation page] is up-to-date and the DOI link points to the most recent release.

## Test the release

 - [ ] After activating a new virtual or conda environment, make sure that the released version installs correctly with `pip install --upgrade plasmapy`.
 - [ ] Open Python and run `import plasmapy`, `dir(plasmapy)`, and `plasmapy.__version__`.
 - [ ] Run `plasma-calculator` from the command line.
 - [ ] Verify that the new version can be installed with conda.

## After the release

 - [ ] Announce the release at the [community meeting].
 - [ ] [Create an issue for the release] to occur ∼3–4 months after this one.
 - [ ] Update the [release checklist], as needed.
 - [ ] Close this issue.
