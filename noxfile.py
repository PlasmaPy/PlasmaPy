"""
Running `nox` without arguments will run tests with the version of
Python that `nox` is installed under, skipping slow tests. To invoke a
nox session, run `nox -s <session>`, where <session> is replaced with
the name of the session. To list available sessions, run `nox -l`.

The tests can be run with the following options:

* "all": run all tests
* "skipslow": run tests, except tests decorated with `@pytest.mark.slow`
* "cov": run all tests with code coverage checks
* "lowest-direct" : run all tests with lowest version of direct dependencies

Doctests are run only for the most recent versions of Python and
PlasmaPy dependencies, and not when code coverage checks are performed.
"""

# nox documentation: https://nox.thea.codes

import sys

import nox

nox.options.default_venv_backend = "uv|virtualenv"


supported_python_versions = ("3.10", "3.11", "3.12")

maxpython = max(supported_python_versions)
minpython = min(supported_python_versions)

current_python = f"{sys.version_info.major}.{sys.version_info.minor}"
nox.options.sessions = [f"tests-{current_python}(skipslow)"]


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

    session.install("uv >= 0.1.39")

    category_version_resolution = [
        ("tests", version, resolution)
        for version in supported_python_versions
        for resolution in ("highest", "lowest-direct")
    ]

    category_flags = {
        "all": ("--all-extras",),
        "docs": ("--extra", "docs"),
        "tests": ("--extra", "tests"),
    }

    command = (
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


pytest_command = (
    "pytest",
    "--pyargs",
    "--durations=5",
    "--tb=short",
    "-n=auto",
    "--dist=loadfile",
)

with_doctests = ("--doctest-modules", "--doctest-continue-on-failure")

with_coverage = (
    "--cov=plasmapy",
    "--cov-config=pyproject.toml",
    "--cov-append",
    "--cov-report xml:coverage.xml",
)

skipslow = ("-m", "not slow")

test_selections = [
    nox.param("all", id="all"),
    nox.param("with code coverage", id="cov"),
    nox.param("skip slow tests", id="skipslow"),
    nox.param("lowest-direct", id="lowest-direct"),
]


@nox.session(python=supported_python_versions)
@nox.parametrize("test_selection", test_selections)
# @nox.parametrize(
#    ["test_selection", "coverage"],
# )
def tests(session, test_selection):
    """Run tests with pytest."""

    resolution = "lowest-direct" if test_selection == "lowest-direct" else "highest"

    requirements = get_requirements_filepath(
        category="tests",
        version=session.python,
        resolution=resolution,
    )

    options = []

    if test_selection == "skipslow":
        options += skipslow

    if test_selection == "cov":
        options += with_coverage

    # Doctests are only run with the most recent versions of Python and
    # other dependencies because there may be subtle differences in the
    # output between different versions of Python, NumPy, and Astropy.
    if session.python == maxpython and test_selection in {"all", "skipslow"}:
        options += with_doctests

    session.install("-r", requirements)
    session.install(".")
    session.run(*pytest_command, *options, *session.posargs)


sphinx_commands = (
    "sphinx-build",
    "docs/",
    "docs/build/html",
    "--nitpicky",
    "--fail-on-warning",
    "--keep-going",
)

html = ("-b", "html")
check_hyperlinks = ("-b", "linkcheck", "-q")
docs_requirements = get_requirements_filepath(category="docs", version=maxpython)


@nox.session(python=maxpython)
def docs(session):
    """Build documentation with most recent supported version of Python."""
    session.install("-r", docs_requirements)
    session.install(".")
    session.run(*sphinx_commands, *html, *session.posargs)


@nox.session(python=maxpython)
def linkcheck(session):
    """
    Check hyperlinks in documentation.

    Use ``linkcheck_ignore`` and ``linkcheck_allowed_redirects`` in
    :file:`docs/conf.py` to specify hyperlink patterns that should be
    ignored.
    """
    session.install("-r", docs_requirements)
    session.install(".")
    session.run(*sphinx_commands, *check_hyperlinks, *session.posargs)


@nox.session
def mypy(session):
    """Perform static type checking."""
    mypy_command = (
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
    session.run(*mypy_command, *session.posargs)
