"""
Configuration file for nox.

To see the list of available environments, run:

   nox --list

To run an environment, use:

   nox -e <environment>

where `<environment>` is the name of the environment to run.
"""

import nox

supported_python_versions = ("3.10", "3.11", "3.12")

minpython = min(supported_python_versions)
maxpython = max(supported_python_versions)

doc_requirements = ("-r", "ci_requirements/requirements_docs_py312.txt")


nox.options.default_venv_backend = "uv|virtualenv"


@nox.session(python=supported_python_versions)
def tests(session):
    """Run tests with pytest."""


# Documentation builds

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


@nox.session(python=supported_python_versions[-1])
def docs(session):
    """Build documentation with Sphinx."""
    session.install(*doc_requirements)
    session.install(".")
    session.run(*sphinx_commands, *html, *session.posargs)


@nox.session(python=supported_python_versions[-1])
def linkcheck(session):
    """Check hyperlinks in documentation."""
    session.install(*doc_requirements)
    session.install(".")
    session.run(*sphinx_commands, *check_hyperlinks, *session.posargs)


# Static type checking


@nox.session(python=supported_python_versions[-1])
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
