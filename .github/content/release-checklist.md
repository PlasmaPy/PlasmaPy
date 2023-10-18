[community meeting]: https://www.plasmapy.org/meetings/weekly
[`docs/about/citation.rst`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/about/citation.rst
[release checklist]: https://github.com/PlasmaPy/PlasmaPy/tree/main/.github/markdown/release-checklist.md
[Zenodo]: https://zenodo.org

# Release checklist

This issue contains the procedure for releasing a new minor version of
PlasmaPy. Below, `2023.10.0` represents the version to be released.

> [!IMPORTANT]
> Split up pre-release changes into multiple focused PRs rather than a
> single monolothic PR.

## Ahead of release

 - [x] [Create an issue with release checklist]()
 - [ ] Announce the forthcoming release at a [community meeting]

### Update metadata

 - [ ] Begin an upload to [Zenodo] for the new release using the
   `team@plasmapy.org` login, and reserve a DOI
 - [ ] Open a pull request to update the DOI in [`docs/about/citation.rst`]
 - [ ] Update `.mailmap`

### Perform code quality checks

 - [ ] Build the docs using `make linkcheck`, and update any broken links
 - [ ] Open a pull request to re-execute pre-executed notebooks, as
       those for charged particle radiography, if there have been any
       substantive changes in the relevant functionality

## Day of release

### Update requirements

 - [ ] Go to the [Actions] page
 - [ ] ...

### Make sure that all tests are passing

 - [ ] Go to the [Actions] page
 - [ ] Click on the "CI" tab and select "Run workflow"
 - [ ] Click on the "weekly tests" tab and select "Run workflow"
 - [ ] Enjoy life for 15 minutes
 - [ ] Fix any failures, and repeat these steps

### Create the release branch

### Publish the release

Go to the GitHub page to [draft a new release].


## After release

 - [ ] Discuss the release at the [community meeting]
 - [ ] Modify the [release checklist], if necessary
