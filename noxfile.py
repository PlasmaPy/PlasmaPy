import nox

nox.options.sessions = ["tests", "linters", "codespell"]

python_versions = ("3.9", "3.10")
newest_python_version = python_versions[-1]

pytest_options = ["--showlocals"]


@nox.session(python=python_versions)
def tests(session):
    session.install(".[tests]")
    session.run("pytest", *pytest_options)


@nox.session
def linters(session):
    session.install(".[tests]")
    flake8_options = ["--count", "--show-source", "--statistics"]
    session.run("flake8", "plasmapy", *flake8_options, *session.posargs)


@nox.session
def import_package(session):
    session.install(".")
    session.run("python", "-c", "import plasmapy")


sphinx_paths = ["docs", "docs/_build/html"]
sphinx_fail_on_warnings = ["-W", "--keep-going"]
sphinx_builder = ["-b", "html", "-n"]
sphinx_opts = sphinx_paths + sphinx_fail_on_warnings + sphinx_builder
sphinx_skip_notebooks = ["-D", "nbsphinx_execute=never"]

post_doc_build_comments = """
For troubleshooting documentation builds, check out:
https://docs.plasmapy.org/en/latest/contributing/doc_guide.html#troubleshooting
"""


@nox.session
def build_docs(session):
    session.install(".[docs]")
    session.run(
        "sphinx-build",
        *sphinx_opts,
        *session.posargs,
    )
    print(post_doc_build_comments)


@nox.session
def build_docs_no_examples(session):
    session.install(".[docs]")
    session.run(
        "sphinx-build",
        *sphinx_opts,
        *sphinx_skip_notebooks,
        *session.posargs,
    )
    print(post_doc_build_comments)
