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
* "lowest-direct" : run all tests with lowest versions of direct dependencies
* "lowest-direct-skipslow" : run non-slow tests with lowest versions of direct dependencies

Doctests are run only for the most recent versions of Python and
PlasmaPy dependencies, and not when code coverage checks are performed.
Some of the checks require the most recent supported version of Python
to be installed.

Nox documentation: https://nox.thea.codes
"""

import os
import pathlib
import re
import shutil
import sys
import tomllib

import nox
from packaging.requirements import Requirement

# SPEC 0 indicates that scientific Python packages should support
# versions of Python that have been released in the last 3 years, or
# equivalently the most three recently released versions of Python.
# The minimum version of Python should be incremented immediately
# following the first release after October of each year.

supported_python_versions: tuple[str, ...] = ("3.11", "3.12", "3.13")
supported_operating_systems: tuple[str, ...] = ("linux", "macos", "windows")

maxpython = max(supported_python_versions)
minpython = min(supported_python_versions)

_HERE = pathlib.Path(__file__).parent

# The documentation should be build always using the same version of
# Python, which should be the latest version of Python supported by Read
# the Docs. Because Read the Docs takes some time to support new
# releases of Python, we should not link docpython to maxpython.

docpython = "3.12"

current_python = f"{sys.version_info.major}.{sys.version_info.minor}"
nox.options.sessions = [f"tests-{current_python}(skipslow)"]

nox.options.default_venv_backend = "uv"

uv_sync = ("uv", "sync", "--no-progress", "--frozen")

running_on_ci = os.getenv("CI")
running_on_rtd = os.environ.get("READTHEDOCS") == "True"


def _create_requirements_pr_message(uv_output: str, session: nox.Session) -> None:
    """
    Create the pull request message during requirements updates.

    This function copies a GitHub flavored Markdown template to a new
    file and appends a table containing the updated requirements, with
    links to the corresponding PyPI pages. This file is then used as the
    body of the pull request message used in the workflow for updating
    requirements.

    Parameters
    ----------
    uv_output : str
        The multi-line output of ``session.run(..., silent=True)``.
    """

    pr_template = pathlib.Path("./.github/content/update-requirements-pr-template.md")
    pr_message = pathlib.Path("./.github/content/update-requirements-pr-body.md")

    shutil.copy(pr_template, pr_message)

    lines = [
        "",
        "| package | old version | new version |",
        "| :-----: | :---------: | :---------: |",
    ]

    for package_update in uv_output.splitlines():
        if not package_update.startswith("Updated"):
            session.debug(f"Line not added to table: {package_update}")
            continue

        try:
            # An example line is "Updated nbsphinx v0.9.6 -> v0.9.7"
            _, package_, old_version_, _, new_version_ = package_update.split()
        except ValueError:
            session.debug(f"Line not added to table: {package_update}:")
            continue

        old_version = f"{old_version_.removeprefix('v')}"
        new_version = f"{new_version_.removeprefix('v')}"

        pypi_link = f"https://pypi.org/project/{package_}/{new_version}"
        package = f"[`{package_}`]({pypi_link})"

        lines.append(f"| {package} | `{old_version}` | `{new_version}` |")

    with pr_message.open(mode="a") as file:
        file.write("\n".join(lines))


def _get_dependencies_from_pyproject_toml(extras: str | None = None):
    _PYTPROJECT_TOML = (_HERE / "pyproject.toml").resolve()
    with _PYTPROJECT_TOML.open(mode="rb") as file:
        data = tomllib.load(file)
        config = data["project"]

    dependencies = {Requirement(item).name: item for item in config["dependencies"]}

    if (
        extras is None
        or "optional-dependencies" not in config
        or not isinstance(extras, str)
        or (extras not in config["optional-dependencies"] and extras != "all")
    ):
        return dependencies

    extras = [extras] if extras != "all" else list(config["optional-dependencies"])
    op_deps = {}
    for extra in extras:
        for dep in config["optional-dependencies"][extra]:
            name = Requirement(dep).name
            op_deps[name] = dep

    return {**dependencies, **op_deps}


@nox.session
def requirements(session: nox.Session) -> None:
    """
    Regenerate the pinned requirements for running tests and building
    documentation.

    This workflow updates :file:`uv.lock` to contain pinned requirements
    for different versions of Python, different operating systems, and
    different dependency sets (i.e., `docs` or `tests`).

    When run in CI, this session will create a file that contains the
    pull request message for the GitHub workflow that updates the pinned
    requirements (:file:`.github/workflows/update-pinned-reqs.yml`).
    """
    uv_lock_upgrade = ["uv", "lock", "--upgrade", "--no-progress"]

    # When silent is `True`, `session.run()` returns a multi-line string
    # with the standard output and standard error.

    uv_output: str | bool = session.run(
        *uv_lock_upgrade,
        *session.posargs,
        silent=running_on_ci,
    )

    if running_on_ci:
        session.log(uv_output)
        _create_requirements_pr_message(uv_output=uv_output, session=session)


@nox.session
def validate_requirements(session: nox.Session) -> None:
    """
    Verify that the requirements in :file:`uv.lock` are compatible
    with the requirements in `pyproject.toml`.
    """
    session.log(
        "🛡 If this check fails, regenerate the pinned requirements in "
        "`uv.lock` with `nox -s requirements`."
    )

    # Generate the cache without updating uv.lock by syncing the
    # current environment. If there ends up being a `--dry-run` option
    # for `uv sync`, we could probably use it here.

    session.run("uv", "sync", "--frozen", "--all-extras", "--no-progress")

    # Verify that uv.lock will be unchanged. Using --offline makes it
    # so that only the information from the cache is used.

    session.run("uv", "lock", "--check", "--offline", "--no-progress")


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
    nox.param("lowest-direct-skipslow", id="lowest-direct-skipslow"),
]


@nox.session(python=supported_python_versions)
@nox.parametrize("test_specifier", test_specifiers)
def tests(session: nox.Session, test_specifier: nox._parametrize.Param) -> None:
    """Run tests with pytest."""

    options: list[str] = []

    if test_specifier in {"skip slow tests", "lowest-direct-skipslow"}:
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

    match test_specifier:
        case "lowest-direct" | "lowest-direct-skipslow":
            session.install(
                ".[tests]",
                "--resolution=lowest-direct",
                "--no-all-extras",
                f"--python={session.virtualenv.location}",
                env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
            )
        case _:
            # From https://nox.thea.codes/en/stable/cookbook.html#using-a-lockfile
            session.run_install(
                *uv_sync,
                "--extra=tests",
                "--no-default-groups",
                f"--python={session.virtualenv.location}",
                env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
            )

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

    if repository == "numpy":
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
    else:
        session.install(f"git+{repository}")

    session.install(".[tests]")
    session.run(*pytest_command, *session.posargs)


if running_on_rtd:
    rtd_output_path = pathlib.Path(os.environ.get("READTHEDOCS_OUTPUT")) / "html"
    rtd_output_path.mkdir(parents=True, exist_ok=True)
    doc_build_dir = str(rtd_output_path)
else:
    doc_build_dir = "docs/_build/html"

sphinx_base_command: list[str] = [
    "sphinx-build",
    "docs/",
    doc_build_dir,
    "--nitpicky",
    "--keep-going",
]

if not running_on_rtd:
    sphinx_base_command.extend(
        [
            "--fail-on-warning",
            "--quiet",
        ]
    )

build_html: tuple[str, ...] = ("--builder", "html")
check_hyperlinks: tuple[str, ...] = ("--builder", "linkcheck")

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

    Configuration file: docs/conf.py
    """

    if running_on_ci:
        session.log(doc_troubleshooting_message)

    session.run_install(
        *uv_sync,
        "--extra=docs",
        "--no-default-groups",
        f"--python={session.virtualenv.location}",
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )
    session.run(*sphinx_base_command, *build_html, *session.posargs)

    landing_page = pathlib.Path(doc_build_dir) / "index.html"

    if landing_page.exists():
        session.log(f"The documentation may be previewed at {landing_page}")
    else:
        session.error(f"Documentation preview landing page not found: {landing_page}")


@nox.session(python=docpython, reuse_venv=True)
def docs_bundle_htmlzip(session: nox.Session) -> None:
    """
    Convert html built docs to a bundle html zip file.
    """

    if not running_on_rtd:
        session.log(
            "Process is NOT being run on Read the Docs.  Will not html ZIP file."
        )
        return None

    html_build_dir = pathlib.Path(doc_build_dir)
    html_landing_page = (html_build_dir / "index.html").resolve()
    READTHEDOCS_OUTPUT = html_build_dir.parent
    if not html_landing_page.exists():
        session.error(
            f"No documentation build found at: {html_landing_page}\n"
            f"It appears the documentation has not been built."
        )

    command = [
        "sphinx-build",
        "--show-traceback",
        "--doctree-dir",
        f"{html_build_dir / '.doctrees'}",
        "--builder",
        "singlehtml",
        "--define",
        "language=en",
        "./docs/",  # source directory
        f"{READTHEDOCS_OUTPUT / 'htmlzip'}",  # output directory
    ]
    session.run(*command)

    # now build the zip file
    READTHEDOCS_PROJECT = os.environ.get("READTHEDOCS_PROJECT")
    READTHEDOCS_LANGUAGE = os.environ.get("READTHEDOCS_LANGUAGE")
    READTHEDOCS_VERSION = os.environ.get("READTHEDOCS_VERSION")
    zip_name = f"{READTHEDOCS_PROJECT}-{READTHEDOCS_LANGUAGE}-{READTHEDOCS_VERSION}.zip"
    # ^ this name mimics how RTD does it by default

    cwd = pathlib.Path.cwd()
    session.chdir(f"{READTHEDOCS_OUTPUT / 'htmlzip'}")
    session.run("zip", "-r", "-m", f"{zip_name}", ".")
    session.chdir(f"{cwd}")

    session.log(f"The htmlzip was placed in: {READTHEDOCS_OUTPUT / 'htmlzip'}")


@nox.session(python=docpython)
@nox.parametrize(
    ["site", "repository"],
    [
        nox.param("github", "sphinx-doc/sphinx", id="sphinx"),
        nox.param("github", "readthedocs/sphinx_rtd_theme", id="sphinx_rtd_theme"),
        nox.param("github", "spatialaudio/nbsphinx", id="nbsphinx"),
        nox.param("github", "plasmapy/plasmapy_sphinx", id="plasmapy_sphinx"),
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
    # Note: Individual dependencies are install in this fashion to
    #       avoid resolution conflicts if an upper dependency limit
    #       had been put on the target package.
    pkg_name = repository.split("/")[-1]
    deps = _get_dependencies_from_pyproject_toml(extras="docs")
    deps.pop(pkg_name, None)

    session.install(
        f"git+https://{site}.com/{repository}",
        *list(deps.values()),
        silent=False,
    )
    session.install("--no-deps", ".")
    session.run(*sphinx_base_command, *build_html, *session.posargs)


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
        session.log(LINKCHECK_TROUBLESHOOTING)

    session.run_install(
        *uv_sync,
        "--extra=docs",
        "--no-default-groups",
        f"--python={session.virtualenv.location}",
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )

    session.run(*sphinx_base_command, *check_hyperlinks, *session.posargs)


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
    """
    Perform static type checking.

    Configuration file: mypy.ini
    """

    session.run_install(
        *uv_sync,
        "--extra=tests",
        "--no-default-groups",
        f"--python={session.virtualenv.location}",
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )

    if running_on_ci:
        session.log(MYPY_TROUBLESHOOTING)

    MYPY_COMMAND: tuple[str, ...] = (
        "mypy",
        ".",
        "--install-types",
        "--non-interactive",
        "--show-error-context",
        "--show-error-code-links",
        "--pretty",
    )

    session.run(*MYPY_COMMAND, *session.posargs)


@nox.session(name="import")
def try_import(session: nox.Session) -> None:
    """Install PlasmaPy and import it."""
    session.install(".")
    session.run("python", "-c", "import plasmapy", *session.posargs)


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
@nox.parametrize("final", [nox.param(False, id="draft"), nox.param(True, id="final")])
def changelog(session: nox.Session, final: str) -> None:
    """
    Build the changelog with towncrier.

     - 'final': build the combined changelog for the release, delete
       the individual changelog entries in `changelog`, and replace
       `CHANGELOG.rst`. Be sure to commit changes before running this
       session.
     - 'draft': print the draft changelog to standard output, without
       writing to files

    When executing this session, provide the version of the release, as
    in this example:

       nox -s 'changelog(final)' -- 2024.7.0
    """

    if len(session.posargs) != 1:
        raise TypeError(
            "Please provide the version of PlasmaPy to be released "
            "(i.e., `nox -s changelog -- 2025.10.0`)"
        )

    version = session.posargs[0]
    year_pattern = r"(202[4-9]|20[3-9][0-9]|2[1-9][0-9]{2}|[3-9][0-9]{3,})"
    month_pattern = r"(1[0-2]|[1-9])"
    patch_pattern = r"(0?[0-9]|[1-9][0-9])"
    version_pattern = rf"^{year_pattern}\.{month_pattern}\.{patch_pattern}$"
    if not re.match(version_pattern, version):
        raise ValueError(
            "Please provide a version of the form YYYY.M.PATCH, where "
            "YYYY is he year, M is the one or two digit month, "
            "and PATCH is a non-negative integer."
        )

    session.install(".", "towncrier")

    towncrier = ["towncrier", "build", "--version", version]

    if not final:
        session.run(*towncrier, "--draft", "--keep")
        return

    original_file = pathlib.Path("./CHANGELOG.rst")
    original_file.unlink()

    session.run(*towncrier, "--yes")

    destination = pathlib.Path(f"./docs/changelog/{version}.rst")
    shutil.copy(original_file, destination)


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
    session.run("python", "-m", "autotyping", *options, *paths, *session.posargs)


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
    """
    Run all pre-commit hooks on all files.

    Configuration file: .pre-commit-config.yaml
    """
    session.install("pre-commit")
    session.run(
        "pre-commit",
        "run",
        "--all-files",
        "--show-diff-on-failure",
        *session.posargs,
    )


zizmor_troubleshooting_message = """

🪧 Run this check locally with `nox -s zizmor` to find potential
security vulnerabilities in GitHub workflows.

📜 Audit rules: https://woodruffw.github.io/zizmor/audits

🔗 If a reported potential vulnerability does not necessitate a fix,
then either append a comment like `# zizmor: ignore[unpinned-uses]` to
the reported line (replacing `unpinned-uses` with the audit rule code),
or add the appropriate configuration settings to: .github/zizmor.yml
"""


@nox.session
def zizmor(session: nox.Session) -> None:
    """
    Find common security issues in GitHub Actions.

    Because some of the zizmor audit rules require a GitHub token,
    running this check locally may produce different results than
    running it in CI.

    Configuration file: .github/zizmor.yml
    """
    if running_on_ci:
        session.log(zizmor_troubleshooting_message)

    session.install("zizmor")
    session.run("zizmor", ".github", "--no-progress", "--color=auto", *session.posargs)


# /// script
# dependencies = ["nox"]
# ///

if __name__ == "__main__":
    nox.main()
