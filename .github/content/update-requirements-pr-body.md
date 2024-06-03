[`ci_requirements/`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/ci_requirements
[`ci_requirements/README.md`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/ci_requirements/README.md
[create an issue]: https://github.com/PlasmaPy/PlasmaPy/issues/new?title=Remove+upper+limit+on+version+of
[`.github/workflows/update-pinned-reqs.yml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.github/workflows/update-pinned-reqs.yml
[`noxfile.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/noxfile.py
[`pyproject.toml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/pyproject.toml

This pull request (PR) regenerates the requirements files in [`ci_requirements/`] that are used for running tests, building documentation, and performing other continuous integration (CI) checks. Pinning requirements files reduces the probability that a CI check will spontaneously start failing in a PR due to a breaking change in a dependency. This PR lets us find unexpected breaking changes in newly released dependencies before they get used in CI. 🛡

**If all checks pass ✅, please merge this PR.** If any checks fail due to a breaking change in a dependency 🚨, please fix the problem and get the checks to pass again.

> [!CAUTION]
> When it is necessary to put a *temporary* upper limit on the allowed versions of a dependency in [`pyproject.toml`], please [create an issue] that this upper limit should be removed.

> [!NOTE]
> For more information, see [`ci_requirements/README.md`]. This workflow is defined in [`.github/workflows/update-pinned-reqs.yml`], and makes use of the `requirements` session defined in [`noxfile.py`].
