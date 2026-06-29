## Instructions

**If all checks pass ✅, please merge this pull request.** If any checks fail due to a breaking change in a dependency 🚨, please address the problems before merging. Thank you!

## Description

This pull request upgrades [`uv.lock`]: the cross-platform lockfile that specifies the Python environments used when running tests, building documentation, and performing continuous integration (CI) checks. Locking and periodically updating the Python environment lets us quarantine breaking changes before they start causing spontaneous failures on unrelated pull requests, while ensuring that everyone is using the same versions of dependencies to perform checks.

## Notes

- This upgrade was performed via `nox --session lock` (see [`noxfile.py`]), which uses [`uv lock --upgrade`] to upgrade `uv.lock` with a [dependency cooldown].

- If it is necessary to temporarily place an upper limit on a dependency in [`pyproject.toml`], please [create an issue] to remove this upper limit before the next release.

- Because the [ty] static type checker rule set is still evolving, this workflow adds `# ty: ignore` comments when new errors arise.

## CI snapshot

[![CI](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci.yml) [![comprehensive tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-comprehensive.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/ci-comprehensive.yml) [![Read the Docs](https://readthedocs.org/projects/plasmapy/badge/?version=latest)](http://plasmapy.readthedocs.io/en/latest/?badge=latest) [![PyHC Actions](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/pyhc-actions.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/pyhc-actions.yml) [![Installability](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/installability.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/installability.yml) [![pre-commit.ci status](https://results.pre-commit.ci/badge/github/PlasmaPy/PlasmaPy/main.svg)](https://results.pre-commit.ci/latest/github/PlasmaPy/PlasmaPy/main)
 [![upgrade uv.lock](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upgrade-uv-lock.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upgrade-uv-lock.yml)

These checks do not need to be fixed here, but should be fixed before releases: [![upstream tests](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-tests.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-tests.yml) [![upstream docs](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-docs.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upstream-docs.yml) [![linkcheck](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/linkcheck.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/linkcheck.yml)

## Dependencies updated by the workflow

[create an issue]: https://github.com/PlasmaPy/PlasmaPy/issues/new?title=Remove+upper+limit+on+version+of
[dependency cooldown]: https://blog.yossarian.net/2025/11/21/We-should-all-be-using-dependency-cooldowns
[ty]: https://docs.astral.sh/ty
[`noxfile.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/noxfile.py#:~:text=def%20lock
[`pyproject.toml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/pyproject.toml
[`uv lock --upgrade`]: https://docs.astral.sh/uv/reference/cli/#uv-lock
[`uv.lock`]: https://docs.astral.sh/uv/guides/projects/#uvlock
