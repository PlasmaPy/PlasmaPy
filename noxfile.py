"""Experimental nox configuration file."""

import nox

nox.options.default_venv_backend = "uv"


maxpython = "3.12"

doc_requirements = ("-r", "requirements/requirements_docs_py312.txt")

sphinx_commands = (
    "sphinx-build",
    "docs/",
    "docs/build/html",
    "--nitpicky",
    "--fail-on-warning",
    "--keep-going",
)

html = ("-b", "html")
check_hyperlinks = ("-b", "linkcheck")


@nox.session(python=maxpython)
def docs(session):
    """Build documentation with Sphinx."""
    session.install(*doc_requirements)
    session.install(".")
    session.run(*sphinx_commands, *html, *session.posargs)


@nox.session(python=maxpython)
def linkcheck(session):
    """Check hyperlinks in documentation."""
    session.install(*doc_requirements)
    session.install(".")
    session.run(*sphinx_commands, *check_hyperlinks, *session.posargs)


@nox.session(python=maxpython)
def mypy(session):
    """Perform static type checking."""
    session.install("mypy >= 1.9.0", "pip")
    session.install("-r", "requirements.txt")
    session.run(
        "mypy",
        ".",
        "--install-types",
        "--non-interactive",
        "--show-error-context",
        "--show-error-code-links",
        "--pretty",
        "--exclude",
        "build",
        "--exclude",
        "_version.py",
        *session.posargs,
    )
