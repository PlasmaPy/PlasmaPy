import nox

nox.options.sessions = ["tests", "linters"]

python_versions = ("3.8", "3.9", "3.10")

sphinx_paths = ["docs", "docs/_build/html"]
sphinx_fail_on_warnings = ["-W", "--keep-going"]
sphinx_builder = ["-b", "html"]
sphinx_opts = sphinx_paths + sphinx_fail_on_warnings + sphinx_builder
sphinx_no_notebooks = ["-D", "nbsphinx_execute=never"]
sphinx_nitpicky = ["-n"]

pytest_options = []


@nox.session(python=python_versions)
def tests(session):
    session.install("-r", "requirements/tests.txt")
    session.install(".")
    session.run("pytest", *pytest_options)


@nox.session
def linters(session):
    session.install("-r", "requirements/tests.txt")
    flake8_options = ["--count", "--show-source", "--statistics"]
    session.run("flake8", "plasmapy", *flake8_options, *session.posargs)


@nox.session
def codespell(session):
    session.install("codespell")
    session.run("codespell", ".")


@nox.session
def import_package(session):
    session.install(".")
    session.run("python", "-c", 'import plasmapy')  # fmt: skip


@nox.session
def build_docs(session):
    session.install("-r", "requirements/docs.txt")
    session.install(".")
    session.run(
        "sphinx-build",
        *sphinx_opts,
        *session.posargs,
    )


@nox.session
def build_docs_nitpicky(session):
    session.install("-r", "requirements/docs.txt")
    session.install(".")
    session.run(
        "sphinx-build",
        *sphinx_opts,
        *sphinx_nitpicky,
        *session.posargs,
    )


@nox.session
def build_docs_no_examples(session):
    session.install("-r", "requirements/docs.txt")
    session.install(".")
    session.run(
        "sphinx-build",
        *sphinx_opts,
        *sphinx_no_notebooks,
        *session.posargs,
    )
