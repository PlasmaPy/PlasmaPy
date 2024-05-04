"""Experimental nox configuration file."""

import nox

nox.options.default_venv_backend = "uv"

supported_python_versions = ["3.10", "3.11", "3.12"]

maxpython = max(supported_python_versions)
minpython = min(supported_python_versions)

requirements_directory = "ci_requirements"


def get_requirements_file(
    category: str, version: str, resolution: str = "highest"
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

    resolution : str
        The resolution strategy to be used by ``uv pip compile``
        ``uv pip compile``.

    Returns
    -------
    str
        The path to the requirements file.
    """
    specifiers = [category, version, resolution]
    return f"{requirements_directory}/{'-'.join(specifiers)}.txt"


@nox.session
def requirements(session):
    """Regenerate pinned requirements files used during CI."""

    session.install("uv >= 0.1.37")

    command = (
        "python",
        "-m",
        "uv",
        "pip",
        "compile",
        "pyproject.toml",
        "--upgrade",
        "--quiet",
    )

    # Generate documentation requirements file for the most recent
    # version of Python and the newest versions of dependencies.

    doc_requirements_file = get_requirements_file(category="docs", version=maxpython)
    session.run(
        *command, "-p", maxpython, "--extra", "docs", "-o", doc_requirements_file
    )

    # Generate testing requirements files for all versions of Python with
    # the newest versions of dependencies.

    for version in supported_python_versions:
        requirements_file = get_requirements_file(category="tests", version=version)

        session.run(
            *command, "-p", version, "--extra", "tests", "-o", requirements_file
        )

    # Generate testing requirements using the lowest-direct resolution strategy

    minimal_requirements_file = get_requirements_file(
        category="tests",
        version=minpython,
        resolution="lowest-direct",
    )

    session.run(
        *command,
        "-p",
        minpython,
        "-o",
        minimal_requirements_file,
        "--resolution=lowest-direct",
    )


# Environments for building documentation

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
documentation_requirements = get_requirements_file(category="docs", version=maxpython)


@nox.session(python=maxpython)
def docs(session):
    """Build documentation with Sphinx."""

    session.install("-r", documentation_requirements)
    session.install(".")
    session.run(*sphinx_commands, *html, *session.posargs)


@nox.session(python=maxpython)
def linkcheck(session):
    """Check hyperlinks in documentation."""
    session.install("-r", documentation_requirements)
    session.install(".")
    session.run(*sphinx_commands, *check_hyperlinks, *session.posargs)


# Environments for static type checking


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
    session.install("mypy >= 1.9.0", "pip")
    session.install("-r", "requirements.txt")
    session.run(*mypy_command, *mypy_options, *session.posargs)
