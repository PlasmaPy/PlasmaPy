"""Experimental nox configuration file."""

import nox

nox.options.default_venv_backend = "uv"


maxpython = "3.12"

doc_requirements = ("-r", "ci_requirements/requirements_docs_py312.txt")

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


pytest = (
    "pytest",
    "--pyargs",
    "--durations=5",
    "--tb=short",
    "-n=auto",
    "--dist=loadfile",
)

slow = ("-m", "slow")
skipslow = ("-m", "not slow")
all = ()

coverage = (
    "--cov=plasmapy",
    "--cov-config=pyproject.toml",
    "--cov-append",
    "--cov-report xml:coverage.xml",
)


@nox.session
@nox.parametrize(
    "test_selection",
    [
        nox.param(all, id="all"),
        nox.param(skipslow, id="skipslow"),
        nox.param(slow, id="slow"),
    ],
)
@nox.parametrize(
    "doctest_option",
    [
        nox.param((), id="no-doctest"),
    ],
)
def tests(session, test_selection, doctest_option):
    """Run pytest for PlasmaPy."""
    session.install("-r", "requirements.txt")
    session.run(*pytest, *session.posargs, *test_selection, *doctest_option)

    # Make it so that doctests are run if and only if python=maxpython
    # and the requirements files are
