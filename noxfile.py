"""
Configuration file for nox.

nox is an automation tool that allows us to configure and perform tasks
using programmable Python sessions. Each nox session is defined via a
function decorated with ``@nox.session``. The name of a session is given
by the name of the function.

To list available sessions, run `nox -l`. To invoke a nox session, run
`nox -s <session>`, where <session> is replaced with the name of the
session.

Documentation for nox is at: https://nox.thea.codes
"""

import nox

nox.options.default_venv_backend = "uv|virtualenv"
nox.options.sessions = []

supported_python_versions = ("3.10", "3.11", "3.12")

maxpython = max(supported_python_versions)
minpython = min(supported_python_versions)


def get_requirements_filepath(
    category: str,
    version: str,
    resolution: str = "highest",
) -> str:
    """
    Return the file path to the requirements file.

    Parameters
    ----------
    category : str
        The category for determining requirements: "docs", "tests", or
        "all".

    version : str
        The version of Python to get the requirements for.

    resolution : str, default: "highest"
        The resolution strategy to be used by ``uv pip compile``. Other
        options include ``"lowest-direct"`` or ``"lowest"``.

    Returns
    -------
    str
        The path to the requirements file.
    """
    requirements_directory = "ci_requirements"
    specifiers = [category, version]
    if resolution != "highest":
        specifiers.append(resolution)
    return f"{requirements_directory}/{'-'.join(specifiers)}.txt"


@nox.session
def requirements(session):
    """Regenerate pinned requirements files used during CI."""

    session.install("uv >= 0.1.39")

    category_version_resolution = [
        ("all", maxpython, "highest"),
        ("docs", maxpython, "highest"),
        ("tests", minpython, "lowest-direct"),
    ]

    category_version_resolution += [
        ("tests", version, "highest") for version in supported_python_versions
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
    )

    flags = (
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
            *flags,
        )


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
    """Build documentation with Sphinx."""
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
    """Perform static type checking with mypy."""
    mypy_command = ("mypy", ".")
    mypy_options = (
        "--install-types",
        "--non-interactive",
        "--show-error-context",
        "--show-error-code-links",
        "--pretty",
    )
    session.install("mypy >= 1.10.0", "pip")
    session.install("-r", "requirements.txt")
    session.run(*mypy_command, *mypy_options, *session.posargs)


@nox.session(name="import")
def try_import(session):
    """Install PlasmaPy and import it."""
    session.install(".")
    session.run("python", "-c", "import plasmapy")


@nox.session
def cff(session):
    """Validate CITATION.cff."""
    session.install("cffconvert")
    session.run("cffconvert", "--validate")


@nox.session
def build(session):
    """Build and verify a source distribution and wheel."""
    session.install("twine", "build")
    build_command = ("python", "-m", "build")
    session.run(*build_command, "--sdist")
    session.run(*build_command, "--wheel")
    session.run("twine", "check", "dist/*")
