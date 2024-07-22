[citation page]: https://docs.plasmapy.org/en/stable/about/citation.html
[`CITATION.cff`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/CITATION.cff
[community meeting]: https://www.plasmapy.org/meetings/weekly
[conda-forge feedstock]: https://github.com/conda-forge/plasmapy-feedstock
[Create an issue for the release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/create-release-issue.yml
[Create a pull request]: https://github.com/PlasmaPy/PlasmaPy/compare
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
[release history]: https://pypi.org/project/plasmapy/#history
[tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/tests.yml
[Update pinned requirements]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/update-pinned-reqs.yml
[weekly tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly.yml
[this record]: https://zenodo.org/doi/10.5281/zenodo.6774349
[Zenodo]: https://zenodo.org/me/uploads
[stable documentation build]: https://docs.plasmapy.org/en/stable
[`stable`]: https://github.com/PlasmaPy/PlasmaPy/tree/stable
[Read the Docs]: https://readthedocs.org/projects/plasmapy
[versions page on RTD]: https://readthedocs.org/projects/plasmapy/versions/

This issue contains the procedure for releasing a new version of PlasmaPy.

### Planning the release

 - [x] [Create an issue for the release].
 - [ ] Update [milestones] for issues & pull requests (PRs). 🛣️

### Code quality updates

 - [ ] Revise changelog entries to make sure that they are understandable, necessary, and correctly categorized. Add the `no changelog entry needed` label to skip doing changelog checks. 📜
 - [ ] Run the [GitHub Action for checking hyperlinks] 🔗, and update broken links. Use `linkcheck_allowed_redirects` in [`docs/conf.py`] to allow redirects (e.g., from `doi.org`). Update or delete the `alias` field for authors in [`CITATION.cff`] who have changed their GitHub username.
 - [ ] Run `git log --format="%aN <%aE>" | sort -u`, and update [`.mailmap`] if there are any duplicate contributors in the output. 📧
 - [ ] [Update pinned requirements]. 📍
 - [ ] Run `pre-commit autoupdate` followed by `pre-commit run --all-files`. Fix any new errors and commit the changes.

### Make sure that all tests are passing on `main`

 - [ ] Run the [tests]. [![CI](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/tests.yml)
 - [ ] Run the [weekly tests]. [![weekly tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly.yml/badge.svg?branch=main)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly-tests.yml)
 - [ ] ✅ Fix any failures and repeat.

### Perform the release

 - [ ] Begin an upload to [Zenodo] for a new version of [this record], using the `team@plasmapy.org` login. Reserve a DOI. 🏷️
 - [ ] Run the GitHub Action to [mint a release] 🍬, specifying the version (i.e., `2024.10.0` or `2024.10.0rc1` for a release candidate) and copying the reserved DOI from Zenodo. This workflow will update metadata, build the changelog, create the release branch (i.e., `2024.10.x`), and tag the release.
 - [ ] [Create a release on GitHub]. 🚀
   - Choose the newly created tag (e.g., `v2024.10.0`), and use it as the title. (The release will be performed from the tag, so it is not necessary to select the branch.) 🏷️
   - Select the option to automatically generate release notes. 📜
   - Select the option to create a discussion for the release under the _General_ category.
   - For official releases, make sure the checkbox is selected for _Set as the latest release_. For beta releases or release candidates (e.g., `v2024.10.0rc1`), specify it as a pre-release.
   - Click on _Publish release_, which will create the GitHub release and trigger the GitHub workflow to [publish to PyPI]. 🚀
   - Check the [release history] on PyPI to make sure that the release was successful. 🆕
 - [ ] [Create a pull request] to merge the `v2024.10.x` branch back into main.
 - [ ] Merge **but do not squash** this PR back into `main`, and **do not delete** the `v2024.10.x` branch.

> [!CAUTION]
> Squashing the pull request back into `main` causes problems by removing the tagged release commit from the history of `main`. For example, `plasmapy.__version__` in development installations is chosen based on the most recent tag in `main`.

### Following the release

 - [ ] Download a `.tar.gz` file of the tagged release from the [list of tagged versions] on GitHub, and upload it to [Zenodo].
   - [ ] Update the author list with new authors from the automatically generated release notes or [`CITATION.cff`].
   - [ ] Update the bibliography, and publish the release to Zenodo.
 - [ ] Fix any problems with the automated pull request to [conda-forge feedstock], if necessary. This step should be automatic, but may take a while.
 - [ ] Update requirements in the [conda-forge feedstock] in `recipe/meta.yaml`, in particular when there is a new version of Python.

### Update documentation

 - [ ] Delete the [`stable`] branch on GitHub, if it exists.
 - [ ] Activate the current and prior release on the [versions page on RTD], if necessary. If the documentation fails to build for a release, activate the corresponding branch (e.g., activate the `v2023.10.x` branch instead of the `v2023.10.0` tag).
 - [ ] Check that the [documentation] builds correctly for the release branch. If the documentation build fails, create a new [`stable`] branch from the release branch (e.g., `2024.10.x`) and fix any problems with the documentation build. The [`stable`] branch is needed if the documentation build for the release fails or if we make any changes to the documentation between releases. The [stable documentation build] will point to the [`stable`] branch on GitHub if it exists. Otherwise, it will point to the most recent release on GitHub. 📚
 - [ ] Verify that the [citation page] is up-to-date and the DOI link points to the most recent release.

## Test the release

 - [ ] After activating a new virtual or conda environment, make sure that the released version installs correctly with `pip install --upgrade plasmapy`.
 - [ ] Open Python and run `import plasmapy`, `dir(plasmapy)`, and `plasmapy.__version__`.
 - [ ] Run `plasma-calculator` from the command line.
 - [ ] Verify that the new version can be installed with conda.

## After the release

 - [ ] Announce the release at the [community meeting].
 - [ ] Update the [release checklist], as needed.
 - [ ] [Create an issue for the release] to occur ∼3–4 months after this one.
 - [ ] Close this issue.
