[community meeting]: https://www.plasmapy.org/meetings/weekly
[Create an issue for the release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/create-release-issue.yml
[`docs/about/citation.rst`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/about/citation.rst
[`docs/conf.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/conf.py
[release checklist]: https://github.com/PlasmaPy/PlasmaPy/tree/main/.github/markdown/release-checklist.md
[Zenodo]: https://zenodo.org
[.mailmap]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.mailmap
[weekly tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly-tests.yml
[tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/tests.yml
[mint a release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/mint-release.yml
[Draft a new release]: https://github.com/PlasmaPy/PlasmaPy/releases/new
[publish to PyPI]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.github/workflows/publish-to-pypi.yml

# Release checklist

This issue contains the procedure for releasing a new minor version of PlasmaPy.

Below, `2023.10.0` represents the version to be released.

> [!IMPORTANT]
> Split up pre-release changes into multiple focused PRs rather than a  single monolothic PR.

## Ahead of release

 - [x] [Create an issue for the release]
 - [ ] Discuss the release timeline at a [community meeting]

<!--
 - [ ] About three weeks before a minor or major release, announce that a feature freeze will occur one week before the anticipated release date. Only pull requests with a limited scope that do not significantly change functionality should be merged during the feature freeze.

 - [ ] Begin a code freeze about three weekdays before a release. Only bugfixes and pull requests that are directly related to the release should be merged during the code freeze.
-->


## Perform code quality checks

 - [ ] Revise changelog entries to make sure that they are understandable, necessary, and correctly categorized.
   - Add the `no changelog entry needed` label to skip doing changelog checks.
 - [ ] Build the docs using `make linkcheck` or via the [weekly tests], and update any broken links.
   - Use `linkcheck_allowed_redirects` in [`docs/conf.py`] to specify allowed redirects.
 - [ ] Open a pull request to re-execute pre-executed notebooks, as those for charged particle radiography.
 - [ ] Update [`.mailmap`] with new contributors.


## Update requirements

 - [ ] Go to the [Actions] page
 - [ ] ...

## Make sure that all tests are passing

 - [ ] Run the [tests]
 - [ ] Run the [weekly tests]
 - [ ] Enjoy life for 15 minutes
 - [ ] Fix any failures, and repeat these steps

<!-- Add badges for the tests & weekly tests? -->

## Perform the release

 - [ ] Begin an upload to [Zenodo] for the new release using the `team@plasmapy.org` login, and reserve a DOI.
 - [ ] Run the GitHub Action to [mint a release]. Specify the version (i.e., `2024.1.0` or `2024.1.0rc` for a release candidate) and copy/paste the reserve DOI from Zenodo.  This action will:

   - Create or checkout the appropriate branch
   - Update [`CITATION.cff`] and [`docs/about/citation.rst`] with the new DOI
   - Build the changelog
   - Tag the release
   - Push the tag, branch state <!-- Add more detail? -->

This action will then trigger the action to [publish to PyPI].

## Test the release

 - [ ] After activating a new virtual or conda environment, run `pip install plasmapy==2024.1.0` to make sure that the new version installs correctly.
 - [ ] Open Python and run `import plasmapy` and `dir(plasmapy)`.
 - [ ] Run `plasma-calculator` from the command line to make sure that the plasma calculator activates.


## Post-release

 - [ ] Announce the release
 - [ ] Discuss the release at the [community meeting]
 - [ ] Modify the [release checklist], as necessary
