This issue contains the procedure for releasing a new version of PlasmaPy.

### Planning the release

- [x] [Create an issue for the release]. 📝
- [ ] Update [milestones] for issues & pull requests (PRs). 🛣️

### Code quality updates (optional, but recommended)

- [ ] Revise changelog entries to make sure that they are understandable, necessary, and correctly categorized. Add the https://github.com/PlasmaPy/PlasmaPy/labels/skip%20changelog%checks label to skip doing changelog checks. 📜
- [ ] Run the [GitHub Action for checking hyperlinks], and update broken links. Use `linkcheck_allowed_redirects` in [`docs/conf.py`] to allow redirects (e.g., from `doi.org`). Update or delete the `alias` field for authors in [`CITATION.cff`] who have changed their GitHub username. 🔗
- [ ] Run `git log --format="%aN <%aE>" | sort -u`, and update [`.mailmap`] if there are any duplicate contributors in the output. 📫
- [ ] [Update pinned requirements] in `uv.lock`. 📍
- [ ] Run `pre-commit autoupdate` followed by `pre-commit run --all-files`. Fix new errors and commit the changes. 🧹

### Perform the release

- [ ] Begin an upload to [Zenodo] for a new version of [this record], using the `team@plasmapy.org` login. Reserve a DOI. 🔢

- [ ] Run the GitHub workflow to [prepare a release], specifying the version (i.e., `2025.8.0`) and copying the reserved DOI from Zenodo. This workflow will create a pull request that builds the changelog and updates package metadata. 🤖

  - [ ] Revise changelog entries to make sure that they are understandable, necessary, and correctly categorized. 📜
  - [ ] Make sure that all tests are passing in the pull request. ✅
  - [ ] Merge the pull request. 📦

- [ ] Run the [tests]. [![CI](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci.yml) 🧪

- [ ] Run the [comprehensive tests]. [![comprehensive tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-comprehensive.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-comprehensive.yml) 🔍

- [ ] Run the [upstream tests]. [![upstream tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-upstream.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-upstream.yml) 🔮

- [ ] [Create a release on GitHub]. 🚀

  - Choose the newly created tag (e.g., `v2025.10.0`), and use it as the title. (The release will be performed from the tag, so it is not necessary to select the branch.) 🏷️
  - Set the tag for the previous release, and select the option to automatically generate release notes. 📜
  - Select the option to create a discussion for the release under the _General_ category. 📣
  - For official releases, make sure the checkbox is selected for _Set as the latest release_. For beta releases or release candidates (e.g., `v2025.10.0rc1`), specify it as a pre-release. 🆕
  - Click on _Publish release_, which will create the GitHub release and trigger the GitHub workflow to [publish to PyPI]. 🚀
  - Check the [release history] on PyPI to make sure that the release was successful. 🗓️

### Following the release

- [ ] Download a `.tar.gz` file of the tagged release from the [list of tagged versions] on GitHub, and upload it to [Zenodo]. 📤
  - [ ] Update the version number and release date in the record. 📅
  - [ ] Update the author list with new authors from the automatically generated release notes or [`CITATION.cff`]. 👥
  - [ ] Update the bibliography. 📖
  - [ ] Publish the record. 🏛️
- [ ] Fix any problems with the automated pull request to [conda-forge feedstock], if necessary. This step should be automatic, but may take a while. 🔧
- [ ] Update requirements in the [conda-forge feedstock] in `recipe/meta.yaml`, in particular when there is a new version of Python. 🔄

> [!TIP]
> To compare two files across different tags, use commands like:
>
> ```shell
> git diff v2024.7.0:CITATION.cff v2024.10.0:CITATION.cff
> git diff v2024.7.0:docs/bibliography.bib v2024.10.0:docs/bibliography.bib
> ```

### Update documentation

- [ ] Delete the [`stable`] branch on GitHub if it exists. 🗑️
- [ ] Activate the current and prior release on the [versions page on RTD], if necessary. If the documentation fails to build for a release, activate the corresponding branch (e.g., activate the `v2025.10.x` branch instead of the `v2025.10.0` tag). ⚙️
- [ ] Verify that the [citation page] is up-to-date and the DOI link points to the most recent release. 🧾
- [ ] Check that the [documentation] builds correctly for the release branch. 📘

> [!TIP]
> If the documentation build fails, create a new [`stable`] branch from the release branch (e.g., `2025.10.x`) and fix any problems with the documentation build. The [`stable`] branch is needed if the documentation build for the release fails or if we make any changes to the documentation between releases. The [stable documentation build] will point to the [`stable`] branch on GitHub if it exists. Otherwise, it will point to the most recent release on GitHub. 📚

## Test the release

- [ ] After activating a new virtual or conda environment, make sure that the released version installs correctly with `pip install --upgrade plasmapy`. 💻
- [ ] Open Python and run `import plasmapy`, `dir(plasmapy)`, and `plasmapy.__version__`. 🐍
- [ ] Verify that the new version can be installed with conda. 🧩
- [ ] Verify that the new version can be installed with uv. 🌈

## After the release

- [ ] Update the [release checklist], as needed. 📋
  - An example changelog entry is: "Updated the release checklist following the `v2025.8.0` release." ✍️
- [ ] [Create an issue for the next release]. ⏳
- [ ] Close this issue. 🏁

[citation page]: https://docs.plasmapy.org/en/stable/about/citation.html
[comprehensive tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-comprehensive.yml
[conda-forge feedstock]: https://github.com/conda-forge/plasmapy-feedstock
[create a release on github]: https://github.com/PlasmaPy/PlasmaPy/releases/new
[create an issue for the next release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/create-release-issue.yml
[documentation]: https://docs.plasmapy/org/en/stable
[github action for checking hyperlinks]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/linkcheck.yml
[list of tagged versions]: https://github.com/PlasmaPy/PlasmaPy/tags
[milestones]: https://github.com/PlasmaPy/PlasmaPy/milestones
[prepare a release]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/prepare-release-pr.yml
[publish to pypi]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.github/workflows/publish-to-pypi.yml
[release checklist]: https://github.com/PlasmaPy/PlasmaPy/tree/main/.github/content/release-checklist.md
[release history]: https://pypi.org/project/plasmapy/#history
[stable documentation build]: https://docs.plasmapy.org/en/stable
[tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/tests.yml
[this record]: https://zenodo.org/doi/10.5281/zenodo.6774349
[update pinned requirements]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/update-pinned-reqs.yml
[upstream tests]: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-comprehensive.yml
[versions page on rtd]: https://readthedocs.org/projects/plasmapy/versions/
[zenodo]: https://zenodo.org/me/uploads
[`.mailmap`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.mailmap
[`citation.cff`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/CITATION.cff
[`docs/conf.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/conf.py
[`stable`]: https://github.com/PlasmaPy/PlasmaPy/tree/stable
