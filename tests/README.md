# Tests

[**testing guide**]: https://docs.plasmapy.org/en/latest/contributing/testing_guide.html
[`src/plasmapy/formulary/speeds.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/src/plasmapy/formulary/speeds.py
[`tests/formulary/test_speeds.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/tests/formulary/test_speeds.py
[`tests/`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/tests
[`src/plasmapy/`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/src/plasmapy
[Nox]: https://nox.thea.codes
[pytest]: https://docs.pytest.org
[running tests]: https://docs.plasmapy.org/en/latest/contributing/testing_guide.html#running-tests

> [!TIP]
> To learn more about software testing in PlasmaPy, please check out the
> [**testing guide**].

## Locating tests

The [`tests/`] directory contains PlasmaPy's tests. The directory
structure and organization of [`tests/`] largely mirrors that of
[`src/plasmapy/`], which contains PlasmaPy's source code. For example,
the source code of `plasmapy.formulary.speeds` is located at
[`src/plasmapy/formulary/speeds.py`], while its tests are located at
[`tests/formulary/test_speeds.py`]

## Running tests

PlasmaPy's tests are written using the [pytest] framework, with [Nox] as
the test runner. The test suite is automatically run whenever a pull
request to PlasmaPy's GitHub repository is submitted or updated.

To run tests locally, first install Nox and its dependencies with:

```shell
python -m pip install nox uv
```

To run all but the slowest tests, enter the top-level directory of your
clone of PlasmaPy and run:

```shell
nox
```

For more details, please see the section of the testing guide on
[running tests].
