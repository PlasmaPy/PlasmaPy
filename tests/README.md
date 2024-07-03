# Tests

[contributor guide]: https://docs.plasmapy.org/en/latest/contributing
[**testing guide**]: https://docs.plasmapy.org/en/latest/contributing/testing_guide.html
[`src/plasmapy/formulary/speeds.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/src/plasmapy/formulary/speeds.py
[`tests/formulary/test_speeds.py`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/tests/formulary/test_speeds.py
[`tests/`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/tests
[`src/plasmapy/`]: https://github.com/PlasmaPy/PlasmaPy/blob/main/src/plasmapy
[Nox]: https://nox.thea.codes
[noxfile.py]: https://github.com/PlasmaPy/PlasmaPy/blob/main/noxfile.py
[pytest]: https://docs.pytest.org

ℹ️To learn more about software testing in PlasmaPy, please check out
the [**testing guide**] in the online [contributor guide].

## Locating tests

The [`tests/`] directory contains PlasmaPy's tests. The directory
structure and organization of [`tests/`] largely mirrors that of
[`src/plasmapy/`], which contains PlasmaPy's source code. For example,
the source code of `plasmapy.formulary.speeds` is located at
[`src/plasmapy/formulary/speeds.py`], while its tests are located at
[`tests/formulary/test_speeds.py`]

## Running tests

PlasmaPy's tests are written using the [pytest] framework, with [Nox] as
the test runner. The test suite (excluding tests marked as slow) can be
run locally with the following steps:

1. Install [Nox] and its dependencies with:
   ```shell
   python -m pip install nox uv
   ```
2. Enter the top-level directory of your clone of PlasmaPy and run
   ```shell
   nox
   ```

Additional Nox sessions are defined in [`noxfile.py`], and can be viewed
with `nox -l`. For example, to run all tests in Python 3.12, run
`nox -s 'tests-3.12(all)'`.
