# PlasmaPy Security Policy

## Supported versions

Security updates will only be provided for the most recent version of PlasmaPy.

| Version | Supported |
| ------- | --------- |
| latest  | ✅        |
| older   | ❌        |

## Finding vulnerabilities via continuous integration

PlasmaPy has an extensive suite of continuous integration checks,
including several that identify and flag common security vulnerabilities.

- [`zizmor`] finds security vulnerabilities in the GitHub workflows.
- [`ruff`] has several [rule sets] that find security vulnerabilities,
  including the [`flake8-bandit` rule set].

The checks are usually either run as [`pre-commit`] hooks
(defined in [`.pre-commit-config.yaml`])
or as [Nox] sessions (defined in [`noxfile.py`])
invoked during GitHub workflows (located in [`.github/workflows/`]).
The configurations for most of these tools are located in [`pyproject.toml`],
but may have a dedicated configuration file
(for `zizmor`, this would be at [`.github/zizmor.yml`])

[`zizmor`]: https://github.com/woodruffw/zizmor
[`.github/workflows/`]: https://github.com/PlasmaPy/PlasmaPy/tree/main/.github/workflows
[`ruff`]: https://docs.astral.sh/ruff
[rule sets]: https://docs.astral.sh/ruff/rules
[`flake8-bandit` rule set]: https://docs.astral.sh/ruff/rules/#flake8-bandit-s
[`pre-commit`]: https://pre-commit.com/
[Nox]: https://nox.thea.codes
[`noxfile.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/noxfile.py
[`.pre-commit-config.yaml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.pre-commit-config.yaml
[`pyproject.toml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/pyproject.toml
[`.github/zizmor.yml`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/.github/zizmor.yml

## Reporting a vulnerability

If you believe you have found a security vulnerability, please [raise an issue] on [PlasmaPy's GitHub repository].

> [!IMPORTANT]
> Please report critical vulnerabilities to team@plasmapy.org instead of creating a GitHub issue.
> When doing so, please include:
>
> - A description of the vulnerability,
> - Steps to reproduce the issue (if applicable), and
> - Any relevant additional information (if needed).

[plasmapy's github repository]: https://github.com/PlasmaPy/PlasmaPy
[raise an issue]: https://github.com/PlasmaPy/PlasmaPy/issues/new/choose
