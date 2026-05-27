## Instructions

**If all checks pass ✅, please merge this pull request.** If any checks fail due to a breaking change in a dependency 🚨, please address the problems before merging.

> [!IMPORTANT]
> When it is necessary to put an upper limit on a dependency in [`pyproject.toml`], please [create an issue] to remove this upper limit before the next release.

## Description

This pull request upgrades [`uv.lock`]: the cross-platform lockfile that specifies the Python environments used when running tests, building documentation, and performing continuous integration (CI) checks.

This upgrade was performed by running `nox --session lock` in the top-level directory of the repository, using the `lock` session defined in [`noxfile.py`]. This [Nox] session uses [`uv lock --upgrade`] to upgrade `uv.lock` with a [dependency cooldown]. [![upgrade uv.lock](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upgrade-uv-lock.yml/badge.svg)](https://github.com/PlasmaPy/PlasmaPy/actions/workflows/upgrade-uv-lock.yml)

> [!NOTE]
> Because the [ty] static type checker rule set is still growing, this workflow adds `# ty: ignore` comments when new errors arise.

## Motivation

Locking requirements lets us prevent spontaneous test failures due to breaking changes in upstream dependencies. This PR lets us quarantine breaking changes in newly released dependencies before they start being used in CI for other pull requests.

## Updated dependencies

[create an issue]: https://github.com/PlasmaPy/PlasmaPy/issues/new?title=Remove+upper+limit+on+version+of
[dependency cooldown]: https://blog.yossarian.net/2025/11/21/We-should-all-be-using-dependency-cooldowns
[nox]: https://nox.thea.codes/en/stable
[ty]: https://docs.astral.sh/ty
[`noxfile.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/noxfile.py
[`pyproject.toml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/pyproject.toml
[`uv lock --upgrade`]: https://docs.astral.sh/uv/reference/cli/#uv-lock
[`uv.lock`]: https://docs.astral.sh/uv/guides/projects/#uvlock
