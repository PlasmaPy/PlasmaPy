git mv src/plasmapy/analysis/tests tests/analysis
git mv src/plasmapy/analysis/swept_langmuir/tests tests/analysis/swept_langmuir
git mv src/plasmapy/analysis/time_series/tests tests/analysis/time_series

git mv src/plasmapy/diagnostics/tests tests/diagnostics
git mv src/plasmapy/diagnostics/charged_particle_radiography/tests tests/diagnostics/charged_particle_radiography


git mv src/plasmapy/dispersion/tests tests/dispersion
git mv src/plasmapy/dispersion/analytical/tests tests/dispersion/analytical
git mv src/plasmapy/dispersion/numerical/tests  tests/dispersion/numerical


git mv src/plasmapy/formulary/tests tests/formulary
git mv src/plasmapy/formulary/collisions/tests tests/formulary/collisions
git mv src/plasmapy/formulary/collisions/helio/tests tests/formulary/collisions/helio

git mv src/plasmapy/particles/tests tests/particles

git mv src/plasmapy/plasma/tests tests/plasma
git mv src/plasmapy/plasma/sources/tests tests/plasma/sources

git mv src/plasmapy/simulation/tests tests/simulation

git mv src/plasmapy/tests tests/tests
git mv src/plasmapy/tests/_helpers/test tests/tests/_helpers

git mv src/plasmapy/utils/tests tests/utils
git mv src/plasmapy/utils/data/tests tests/utils/data
git mv src/plasmapy/utils/decorators/tests tests/utils/decorators
git mv src/plasmapy/utils/_pytest_helpers/tests tests/utils/_pytest_helpers

git add .
git commit -m "Move tests from src/plasmapy/**/tests to tests/" -n
