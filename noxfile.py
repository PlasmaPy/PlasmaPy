"""
Nox is an automation tool used by PlasmaPy to run tests, build
documentation, and perform other checks. Nox sessions are defined in
noxfile.py.

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
Some of the checks require the most recent supported version of Python
to be installed.
"""

# Documentation: https://nox.thea.codes
import os
import pathlib
import platform
import re
import sys
from typing import Literal

import nox

supported_python_versions: tuple[str, ...] = ("3.10", "3.11", "3.12", "3.13")
supported_operating_systems: tuple[str, ...] = ("linux", "macos", "windows")

maxpython = max(supported_python_versions)
minpython = min(supported_python_versions)

current_python = f"{sys.version_info.major}.{sys.version_info.minor}"

# The documentation should be build always using the same version of
# Python, which should be the latest version of Python supported by Read
# the Docs. Because Read the Docs takes some time to support new
# releases of Python, we should not link docpython to maxpython.

docpython = "3.12"

nox.options.sessions: list[str] = [f"tests-{current_python}(skipslow)"]
nox.options.default_venv_backend = "uv|virtualenv"

running_on_ci = os.getenv("CI")


def _get_requirements_filepath(
    category: Literal["docs", "tests", "all"],
    version: Literal["3.10", "3.11", "3.12", "3.13", "3.14", "3.15"],
    resolution: Literal["highest", "lowest-direct", "lowest"] = "highest",
    os_platform: Literal["linux", "macos", "windows"] | None = None,
) -> str:
    """
    Return the file path to the requirements file.

    Parameters
    ----------
    category : str
        The name of the optional dependency set, as defined in
        :file:`pyproject.toml`.

    version : str
        The supported version of Python.

    resolution : str
        The resolution strategy used by uv.

    os_platform : str, optional
        The name of the target platform. By default, it will attempt to find the
        requirement file associated with the current platform.
    """

    if os_platform is None:
        current_platform = platform.system().lower()
        os_platform = (
            current_platform
            if current_platform in supported_operating_systems
            else "linux"
        )

    requirements_directory = "ci_requirements"
    specifiers = [category, version, os_platform]
    if resolution != "highest":
        specifiers.append(resolution)
    return f"{requirements_directory}/{'-'.join(specifiers)}.txt"


@nox.session
def requirements(session) -> None:
    """
    Regenerate the pinned requirements files used in CI.

    This session uses `uv pip compile` to regenerate the pinned
    requirements files in `ci_requirements/` for use by the Nox sessions
    for running tests, building documentation, and performing other
    continuous integration checks.
    """

    session.install("uv")

    category_version_resolution: list[tuple[str, str, str]] = [
        ("tests", version, "highest") for version in supported_python_versions
    ]

    category_version_resolution += [
        ("tests", minpython, "lowest-direct"),
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

    for os_platform in supported_operating_systems:
        for category, version, resolution in category_version_resolution:
            filename = _get_requirements_filepath(
                category, version, resolution, os_platform
            )
            session.run(
                *command,
                "--python-version",
                version,
                *category_flags[category],
                "--output-file",
                filename,
                "--resolution",
                resolution,
                *session.posargs,
                "--python-platform",
                os_platform,
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
    nox.param("skip slow tests", id="skipslow"),
    nox.param("with code coverage", id="cov"),
    nox.param("lowest-direct", id="lowest-direct"),
]


@nox.session(python=supported_python_versions)
@nox.parametrize("test_specifier", test_specifiers)
def tests(session: nox.Session, test_specifier: nox._parametrize.Param) -> None:
    """Run tests with pytest."""

    resolution = "lowest-direct" if test_specifier == "lowest-direct" else "highest"

    requirements = _get_requirements_filepath(
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
    if session.python == maxpython and test_specifier not in {"lowest-direct", "cov"}:
        options += with_doctests

    if gh_token := os.getenv("GH_TOKEN"):
        session.env["GH_TOKEN"] = gh_token

    session.install("-r", requirements, ".[tests]")
    session.run(*pytest_command, *options, *session.posargs)


@nox.session(python=maxpython)
@nox.parametrize(
    ["repository"],
    [
        nox.param("numpy", id="numpy"),
        nox.param("https://github.com/astropy/astropy", id="astropy"),
        nox.param("https://github.com/pydata/xarray", id="xarray"),
        nox.param("https://github.com/lmfit/lmfit-py", id="lmfit"),
        nox.param("https://github.com/pandas-dev/pandas", id="pandas"),
    ],
)
def run_tests_with_dev_version_of(session: nox.Session, repository: str) -> None:
    """
    Run tests against the development branch of a dependency.

    Running this session helps us catch problems resulting from breaking
    changes in an upstream dependency before its official release.
    """
    if repository != "numpy":
        session.install(f"git+{repository}")
    else:
        # From: https://numpy.org/doc/1.26/dev/depending_on_numpy.html
        session.run_install(
            "uv",
            "pip",
            "install",
            "-U",
            "--pre",
            "--only-binary",
            ":all:",
            "-i",
            "https://pypi.anaconda.org/scientific-python-nightly-wheels/simple",
            "numpy",
        )
    session.install(".[tests]")
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

build_html: tuple[str, ...] = ("--builder", "html")
check_hyperlinks: tuple[str, ...] = ("--builder", "linkcheck")
docs_requirements = _get_requirements_filepath(category="docs", version=maxpython)

doc_troubleshooting_message = """

📘 Tips for troubleshooting common documentation build failures are in
PlasmaPy's documentation guide at:

🔗 https://docs.plasmapy.org/en/latest/contributing/doc_guide.html#troubleshooting
"""


@nox.session(python=docpython)
def docs(session: nox.Session) -> None:
    """
    Build documentation with Sphinx.

    This session may require installation of pandoc and graphviz.
    """
    if running_on_ci:
        session.debug(doc_troubleshooting_message)
    session.install("-r", docs_requirements, ".[docs]")
    session.run(*sphinx_commands, *build_html, *session.posargs)
    landing_page = (
        pathlib.Path(session.invoked_from) / "docs" / "build" / "html" / "index.html"
    )

    if not running_on_ci and landing_page.exists():
        session.debug(f"The documentation may be previewed at {landing_page}")
    elif not running_on_ci:
        session.debug(f"Documentation preview landing page not found: {landing_page}")


@nox.session(python=docpython)
@nox.parametrize(
    ["site", "repository"],
    [
        nox.param("github", "sphinx-doc/sphinx", id="sphinx"),
        nox.param("github", "readthedocs/sphinx_rtd_theme", id="sphinx_rtd_theme"),
        nox.param("github", "spatialaudio/nbsphinx", id="nbsphinx"),
    ],
)
def build_docs_with_dev_version_of(
    session: nox.Session, site: str, repository: str
) -> None:
    """
    Build documentation against the development branch of a dependency.

    The purpose of this session is to catch bugs and breaking changes
    so that they can be fixed or updated earlier rather than later.
    """
    session.install(f"git+https://{site}.com/{repository}", ".[docs]")
    session.run(*sphinx_commands, *build_html, *session.posargs)


LINKCHECK_TROUBLESHOOTING = """
The Sphinx configuration variables `linkcheck_ignore` and
`linkcheck_allowed_redirects` in `docs/conf.py` can be used to specify
hyperlink patterns to be ignored along with allowed redirects. For more
information, see:

🔗 https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-linkcheck_ignore
🔗 https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-linkcheck_allowed_redirects

These variables are in the form of Python regular expressions:

🔗 https://docs.python.org/3/howto/regex.html
"""


@nox.session(python=docpython)
def linkcheck(session: nox.Session) -> None:
    """Check hyperlinks in documentation."""
    if running_on_ci:
        session.debug(LINKCHECK_TROUBLESHOOTING)
    session.install("-r", docs_requirements, ".[docs]")
    session.run(*sphinx_commands, *check_hyperlinks, *session.posargs)


MYPY_TROUBLESHOOTING = """
🛡 To learn more about type hints, check out mypy's cheat sheet at:
  https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html

For more details about specific mypy errors, go to:
🔗 https://mypy.readthedocs.io/en/stable/error_codes.html

🪧 Especially difficult errors can be ignored with an inline comment of
the form: `# type: ignore[error]`, where `error` is replaced with the
mypy error code. Please use sparingly!

🛠 To automatically add type hints for common patterns, run:
  nox -s 'autotyping(safe)'
"""


@nox.session(python=maxpython)
def mypy(session: nox.Session) -> None:
    """Perform static type checking."""
    if running_on_ci:
        session.debug(MYPY_TROUBLESHOOTING)
    MYPY_COMMAND: tuple[str, ...] = (
        "mypy",
        ".",
        "--install-types",
        "--non-interactive",
        "--show-error-context",
        "--show-error-code-links",
        "--pretty",
    )

    requirements = _get_requirements_filepath(
        category="tests",
        version=session.python,
        resolution="highest",
    )
    session.install("pip")
    session.install("-r", requirements, ".[tests]")
    session.run(*MYPY_COMMAND, *session.posargs)


@nox.session(name="import")
def try_import(session: nox.Session) -> None:
    """Install PlasmaPy and import it."""
    session.install(".")
    session.run("python", "-c", "import plasmapy", *session.posargs)


@nox.session
def validate_requirements(session: nox.Session) -> None:
    """Verify that the pinned requirements are consistent with pyproject.toml."""
    requirements_file = _get_requirements_filepath(
        category="all",
        version=maxpython,
        resolution="highest",
    )
    session.install("uv")
    session.debug(
        "🛡 If this check fails, regenerate the pinned requirements files "
        "with `nox -s requirements` (see `ci_requirements/README.md`)."
    )
    session.run(
        "uv",
        "pip",
        "install",
        "-r",
        requirements_file,
        ".[docs,tests]",
        "--dry-run",
    )


@nox.session
def build(session: nox.Session) -> None:
    """Build & verify the source distribution and wheel."""
    session.install("twine", "build")
    build_command = ("python", "-m", "build")
    session.run(*build_command, "--sdist")
    session.run(*build_command, "--wheel")
    session.run("twine", "check", "dist/*", *session.posargs)


AUTOTYPING_SAFE: tuple[str, ...] = (
    "--none-return",
    "--scalar-return",
    "--annotate-magics",
)
AUTOTYPING_RISKY: tuple[str, ...] = (
    *AUTOTYPING_SAFE,
    "--bool-param",
    "--int-param",
    "--float-param",
    "--str-param",
    "--bytes-param",
    "--annotate-imprecise-magics",
)


@nox.session
@nox.parametrize("draft", [nox.param(False, id="draft"), nox.param(True, id="final")])
def changelog(session: nox.Session, final: str) -> None:
    """
    Build the changelog with towncrier.

     - 'final': build the combined changelog for the release, and delete
       the individual changelog entries in `changelog`.
     - 'draft': print the draft changelog to standard output, without
       writing to files

    When executing this session, provide the version of the release, as
    in this example:

       nox -s 'changelog(final)' -- 2024.7.0
    """

    if len(session.posargs) != 1:
        raise TypeError(
            "Please provide the version of PlasmaPy to be released "
            "(i.e., `nox -s changelog -- 2024.9.0`"
        )

    source_directory = pathlib.Path("./changelog")

    extraneous_files = source_directory.glob("changelog/*[0-9]*.*.rst?*")
    if final and extraneous_files:
        session.error(
            "Please delete the following extraneous files before "
            "proceeding, as the presence of these files may cause "
            f"towncrier errors: {extraneous_files}"
        )

    version = session.posargs[0]

    year_pattern = r"(202[4-9]|20[3-9][0-9]|2[1-9][0-9]{2}|[3-9][0-9]{3,})"
    month_pattern = r"(1[0-2]|[1-9])"
    patch_pattern = r"(0?[0-9]|[1-9][0-9])"
    version_pattern = rf"^{year_pattern}\.{month_pattern}\.{patch_pattern}$"

    if not re.match(version_pattern, version):
        raise ValueError(
            "Please provide a version of the form YYYY.M.PATCH, where "
            "YYYY is the year past 2024, M is the one or two digit month, "
            "and PATCH is a non-negative integer."
        )

    session.install(".", "towncrier")

    options = ("--yes",) if final else ("--draft", "--keep")

    session.run(
        "towncrier",
        "build",
        "--config",
        "pyproject.toml",
        "--dir",
        ".",
        "--version",
        version,
        *options,
    )

    if final:
        original_file = pathlib.Path("./CHANGELOG.rst")
        destination = pathlib.Path(f"./docs/changelog/{version}.rst")
        original_file.rename(destination)


@nox.session
@nox.parametrize(
    "options",
    [
        nox.param(AUTOTYPING_SAFE, id="safe"),
        nox.param(AUTOTYPING_RISKY, id="aggressive"),
    ],
)
def autotyping(session: nox.Session, options: tuple[str, ...]) -> None:
    """
    Automatically add type hints with autotyping.

    The `safe` option generates very few incorrect type hints, and can
    be used in CI. The `aggressive` option may add type hints that are
    incorrect, so please perform a careful code review when using this
    option.

    To check specific files, pass them after a `--`, such as:

        nox -s 'autotyping(safe)' -- noxfile.py
    """
    session.install(".[tests,docs]", "autotyping", "typing_extensions")
    DEFAULT_PATHS = ("src", "tests", "tools", "*.py", ".github", "docs/*.py")
    paths = session.posargs or DEFAULT_PATHS
    session.run("python", "-m", "autotyping", *options, *paths)


@nox.session
def monkeytype(session: nox.Session) -> None:
    """
    Add type hints to a module based on variable types from running pytest.

    Examples
    --------
    nox -s monkeytype -- plasmapy.particles.atomic
    """

    if not session.posargs:
        session.error(
            "Please add at least one module using a command like: "
            "`nox -s monkeytype -- plasmapy.particles.atomic`"
        )

    session.install(".[tests]")
    session.install("MonkeyType", "pytest-monkeytype", "pre-commit")

    database = pathlib.Path("./monkeytype.sqlite3")

    if not database.exists():
        session.log(f"File {database.absolute()} not found. Running MonkeyType.")
        session.run("pytest", f"--monkeytype-output={database.absolute()}")
    else:
        session.log(f"File {database.absolute()} found.")

    for module in session.posargs:
        session.run("monkeytype", "apply", module)

    session.run("pre-commit", "run", "ruff", "--all-files")
    session.run("pre-commit", "run", "ruff-format", "--all-files")

    session.log("Please inspect newly added type hints for correctness.")
    session.log("Check new type hints with `nox -s mypy`.")


@nox.session
def cff(session: nox.Session) -> None:
    """Validate CITATION.cff against the metadata standard."""
    session.install("cffconvert")
    session.run("cffconvert", "--validate", *session.posargs)


@nox.session
def manifest(session: nox.Session) -> None:
    """
    Check for missing files in MANIFEST.in.

    When run outside of CI, this check may report files that were
    locally created but not included in version control. These false
    positives can be ignored by adding file patterns and paths to
    `ignore` under `[tool.check-manifest]` in `pyproject.toml`.
    """
    session.install("check-manifest")
    session.run("check-manifest", *session.posargs)


@nox.session
def lint(session: nox.Session) -> None:
    """Run all pre-commit hooks on all files."""
    session.install("pre-commit")
    session.run(
        "pre-commit",
        "run",
        "--all-files",
        "--show-diff-on-failure",
        *session.posargs,
    )
