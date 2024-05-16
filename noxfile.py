"""
Running `nox` without arguments will run tests with the version of
Python that `nox` is installed under, skipping slow tests. To invoke a
nox session, enter the top-level directory of this repository and run
`nox -s "<session>"`, where <session> is replaced with the name of the
session. To list available sessions, run `nox -l`.

The tests can be run with the following options:

* "all": run all tests
* "skipslow": run tests, except tests decorated with `@pytest.mark.slow`
* "cov": run all tests with code coverage checks
* "lowest-direct" : run all tests with lowest version of direct dependencies

Doctests are run only for the most recent versions of Python and
PlasmaPy dependencies, and not when code coverage checks are performed.
"""

import os
import sys

import nox

supported_python_versions: tuple[str, ...] = ("3.10", "3.11", "3.12")

maxpython = max(supported_python_versions)
minpython = min(supported_python_versions)

current_python = f"{sys.version_info.major}.{sys.version_info.minor}"

nox.options.sessions: list[str] = [f"tests-{current_python}(skipslow)"]
nox.options.default_venv_backend = "uv|virtualenv"


def get_requirements_filepath(
    category: str,
    version: str,
    resolution: str = "highest",
) -> str:
    """
    Return the file path to the requirements file.

    Parameters
    ----------
    category : "docs" | "tests" | "all"
    version : "3.10" | "3.11" | "3.12"
    resolution : "highest" (default) | "lowest-direct" | "lowest"
    """
    requirements_directory = "ci_requirements"
    specifiers = [category, version]
    if resolution != "highest":
        specifiers.append(resolution)
    return f"{requirements_directory}/{'-'.join(specifiers)}.txt"


@nox.session
def requirements(session):
    """Regenerate pinned requirements files."""

    session.install("uv >= 0.1.44")

    category_version_resolution: list[tuple[str, str, str]] = [
        ("tests", version, resolution)
        for version in supported_python_versions
        for resolution in ("highest", "lowest-direct")
    ]

    category_version_resolution += [
        ("docs", maxpython, "highest"),
        ("all", maxpython, "highest"),
    ]

    category_flags: dict[str, tuple[str, ...]] = {
        "all": ("--all-extras",),
        "docs": ("--extra", "docs"),
        "tests": ("--extra", "tests"),
    }

    command: tuple[str, ...] = (
        "python",
        "-m",
        "uv",
        "pip",
        "compile",
        "pyproject.toml",
        "--upgrade",
        "--quiet",
        "--custom-compile-command",  # defines command to be included in file header
        "nox -s requirements",
    )

    for category, version, resolution in category_version_resolution:
        filename = get_requirements_filepath(category, version, resolution)
        session.run(
            *command,
            "--python-version",
            version,
            *category_flags[category],
            "--output-file",
            filename,
            "--resolution",
            resolution,
        )


pytest_command: tuple[str, ...] = (
    "pytest",
    "--pyargs",
    "--durations=5",
    "--tb=short",
    "-n=auto",
    "--dist=loadfile",
)

with_doctests: tuple[str, ...] = ("--doctest-modules", "--doctest-continue-on-failure")

with_coverage: tuple[str, ...] = (
    "--cov=plasmapy",
    "--cov-report=xml",
    "--cov-config=pyproject.toml",
    "--cov-append",
    "--cov-report",
    "xml:coverage.xml",
)

skipslow: tuple[str, ...] = ("-m", "not slow")

test_specifiers: list = [
    nox.param("run all tests", id="all"),
    nox.param("with code coverage", id="cov"),
    nox.param("skip slow tests", id="skipslow"),
    nox.param("lowest-direct", id="lowest-direct"),
]


@nox.session(python=supported_python_versions)
@nox.parametrize("test_specifier", test_specifiers)
def tests(session: nox.Session, test_specifier: nox._parametrize.Param):
    """Run tests with pytest."""

    resolution = "lowest-direct" if test_specifier == "lowest-direct" else "highest"

    requirements = get_requirements_filepath(
        category="tests",
        version=session.python,
        resolution=resolution,
    )

    options: list[str] = []

    if test_specifier == "skip slow tests":
        options += skipslow

    if test_specifier == "with code coverage":
        options += with_coverage

    # Doctests are only run with the most recent versions of Python and
    # other dependencies because there may be subtle differences in the
    # output between different versions of Python, NumPy, and Astropy.
    if session.python == maxpython and test_specifier in {"all", "skipslow"}:
        options += with_doctests

    if gh_token := os.getenv("GH_TOKEN"):
        session.env["GH_TOKEN"] = gh_token

    session.install("-r", requirements, ".[tests]")
    session.run(*pytest_command, *options, *session.posargs)


@nox.session(python=maxpython)
@nox.parametrize(
    ["site", "repository"],
    [
        nox.param("github", "numpy/numpy", id="numpy"),
        nox.param("github", "astropy/astropy", id="astropy"),
        nox.param("github", "pydata/xarray", id="xarray"),
        nox.param("github", "lmfit/lmfit-py", id="lmfit"),
        nox.param("github", "pandas-dev/pandas", id="pandas"),
        nox.param("github", "h5py/h5py", id="h5py"),
    ],
)
def run_tests_with_dev_version_of(session: nox.Session, site: str, repository: str):
    """
    Run tests against the development branch of a dependency.

    The purpose of this session is to catch bugs and breaking changes
    so that they can be fixed or updated earlier rather than later.
    """
    session.install(f"git+https://{site}.com/{repository}", ".[tests]")
    session.run(*pytest_command, *session.posargs)


sphinx_commands: tuple[str, ...] = (
    "sphinx-build",
    "docs/",
    "docs/build/html",
    "--nitpicky",
    "--fail-on-warning",
    "--keep-going",
    "-q",
)

build_html: tuple[str, ...] = ("-b", "html")
check_hyperlinks: tuple[str, ...] = ("-b", "linkcheck")
docs_requirements = get_requirements_filepath(category="docs", version=maxpython)

doc_troubleshooting_message = """

To learn how to address common build failures, please check out
PlasmaPy's documentation troubleshooting guide at:

https://docs.plasmapy.org/en/latest/contributing/doc_guide.html#troubleshooting
"""


@nox.session(python=maxpython)
def docs(session: nox.Session):
    """Build the docs."""
    session.debug(doc_troubleshooting_message)
    session.install("-r", docs_requirements, ".")
    session.run(*sphinx_commands, *build_html, *session.posargs)


@nox.session(python=maxpython)
@nox.parametrize(
    ["site", "repository"],
    [
        nox.param("github", "sphinx-doc/sphinx", id="sphinx"),
        nox.param("github", "readthedocs/sphinx_rtd_theme", id="sphinx_rtd_theme"),
        nox.param("github", "spatialaudio/nbsphinx", id="nbsphinx"),
    ],
)
def build_docs_with_dev_version_of(session: nox.Session, site: str, repository: str):
    """
    Build docs against the development branch of a dependency.

    The purpose of this session is to catch bugs and breaking changes
    so that they can be fixed or updated earlier rather than later.
    """
    session.install(f"git+https://{site}.com/{repository}", ".[docs]")
    session.run(*sphinx_commands, *build_html, *session.posargs)


@nox.session(python=maxpython)
def linkcheck(session: nox.Session):
    """
    Check hyperlinks in documentation.

    Use ``linkcheck_ignore`` and ``linkcheck_allowed_redirects`` in
    :file:`docs/conf.py` to specify hyperlink patterns that should be
    ignored.
    """
    session.install("-r", docs_requirements)
    session.install(".")
    session.run(*sphinx_commands, *check_hyperlinks, *session.posargs)


@nox.session(python=maxpython)
def mypy(session: nox.Session):
    """Perform static type checking."""
    mypy_command: tuple[str, ...] = (
        "mypy",
        ".",
        "--install-types",
        "--non-interactive",
        "--show-error-context",
        "--show-error-code-links",
        "--pretty",
    )
    requirements = get_requirements_filepath(
        category="tests",
        version=session.python,
        resolution="highest",
    )
    session.install("-r", requirements, ".[tests]")
    session.run(*mypy_command, *session.posargs)


@nox.session(name="import")
def try_import(session: nox.Session):
    """Install PlasmaPy and import it."""
    session.install(".")
    session.run("python", "-c", "import plasmapy")


@nox.session
def packaging(session: nox.Session):
    """Build and verify a source distribution and wheel."""
    session.install("twine", "build")
    build_command = ("python", "-m", "build")
    session.run(*build_command, "--sdist")
    session.run(*build_command, "--wheel")
    session.run("twine", "check", "dist/*")


@nox.session
def cff(session: nox.Session):
    """Validate CITATION.cff."""
    session.install("cffconvert")
    session.run("cffconvert", "--validate")
