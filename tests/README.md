# Tests

[contributor guide]: https://docs.plasmapy.org/en/latest/contributing
[**testing guide**]: https://docs.plasmapy.org/en/latest/contributing/testing_guide.html
[`src/plasmapy/formulary/speeds.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/src/plasmapy/formulary/speeds.py
[`tests/formulary/test_speeds.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/tests/formulary/test_speeds.py
[`tests`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/tests
[`src/plasmapy/`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/src/plasmapy

The [`tests`] directory contains PlasmaPy's tests.

> [!TIP]
> To learn more about software testing, please check out the
> [**testing guide**] in PlasmaPy's [contributor guide].

## Locating tests

The directory structure and organization of [`tests`] largely mirrors
that of [`src/plasmapy/`], which contains PlasmaPy's source code.

Using the example of `plasmapy.formulary.speeds`:

 - The source code is located at [`src/plasmapy/formulary/speeds.py`]
 - Tests are located at [`tests/formulary/test_speeds.py`]

## Running tests

[Nox]: https://nox.thea.codes
[pytest]: https://docs.pytest.org

PlasmaPy uses the [pytest] framework with [Nox] as its test runner. The
test suite (excluding tests marked as slow) can be run locally with the
following steps:

1. Install [Nox] and its dependencies with:
   ```shell
   python -m pip install nox uv
   ```
2. Enter the top-level directory of your clone of PlasmaPy and run
   ```shell
   nox
   ```

The [**testing guide**] describes the different Nox sessions more
thoroughly.
