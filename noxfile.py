"""
Running `nox` without arguments will run tests with the version of
Python that `nox` is installed under, skipping slow tests. To invoke a
nox session, run `nox -s "<session>"`, where <session> is replaced with
the name of the session. To list available sessions, run `nox -l`.

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

    session.install("-r", requirements)
    session.install(".")
    session.run(*pytest_command, *options, *session.posargs)


sphinx_commands: tuple[str, ...] = (
    "sphinx-build",
    "docs/",
    "docs/build/html",
    "--nitpicky",
    "--fail-on-warning",
    "--keep-going",
)

html: tuple[str, ...] = ("-b", "html")
check_hyperlinks: tuple[str, ...] = ("-b", "linkcheck", "-q")
docs_requirements = get_requirements_filepath(category="docs", version=maxpython)


@nox.session(python=maxpython)
def docs(session: nox.Session):
    """Build documentation with most recent supported version of Python."""
    session.install("-r", docs_requirements)
    session.install(".")
    session.run(*sphinx_commands, *html, *session.posargs)


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


MYPY_TROUBLESHOOTING = """

To learn more about type hints 🦺, check out mypy's cheat sheet at:

  https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html

For more details about specific mypy errors, go to:

  https://mypy.readthedocs.io/en/stable/error_codes.html

Errors can be ignored on individual lines by adding a comment such as:

  # type: ignore[assignment]

where `assignment` is the mypy error code. Please use `# type: ignore`
judiciously!

To automatically add type hints for common patterns, run:

  nox -s 'autotyping(safe)'
"""


@nox.session
def mypy(session: nox.Session):
    """Perform static type checking."""
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
    session.install("mypy >= 1.10.0", "pip")
    session.install("-r", "requirements.txt")
    session.run(*MYPY_COMMAND, *session.posargs)


@nox.session(name="import")
def try_import(session: nox.Session):
    """Install PlasmaPy and import it."""
    session.install(".")
    session.run("python", "-c", "import plasmapy")


@nox.session
def build(session: nox.Session):
    """Build and verify a source distribution and wheel."""
    session.install("twine", "build")
    build_command = ("python", "-m", "build")
    session.run(*build_command, "--sdist")
    session.run(*build_command, "--wheel")
    session.run("twine", "check", "dist/*")


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
def cff(session: nox.Session):
    """Validate CITATION.cff."""
    session.install("cffconvert")
    session.run("cffconvert", "--validate")
