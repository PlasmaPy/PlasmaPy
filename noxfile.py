# /// script
# requires-python = ">=3.12"
# dependencies = ["nox", "nox-uv", "uv"]
# ///

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

import nox
import nox_uv

# SPEC 0 indicates that scientific Python packages should support
# versions of Python that have been released in the last 3 years, or
# equivalently the most three recently released versions of Python.
# The minimum version of Python should be incremented immediately
# following the first release after October of each year.

SUPPORTED_PYTHON_VERSIONS: tuple[str, ...] = ("3.12", "3.13", "3.14")
SUPPORTED_OPERATING_SYSTEMS: tuple[str, ...] = ("linux", "macos", "windows")

MAXPYTHON = max(SUPPORTED_PYTHON_VERSIONS)
MINPYTHON = min(SUPPORTED_PYTHON_VERSIONS)

ROOT_DIR = pathlib.Path(__file__).parent

CURRENT_PYTHON = f"{sys.version_info.major}.{sys.version_info.minor}"

# Define what sessions get run when running `nox` without a session

nox.options.sessions = [f"tests-{CURRENT_PYTHON}(all)"]

# The documentation should be build always using the same version of
# Python, which should be the latest version of Python supported by Read
# the Docs. Because Read the Docs takes some time to support new
# releases of Python, DOCPYTHON should stay independent of MAXPYTHON.
# Changing DOCPYTHON also requires updating .readthedocs.yml and the
# GitHub workflows for building the documentation.

DOCPYTHON = "3.14"

nox.options.default_venv_backend = "uv"

UV_SYNC = ("uv", "sync", "--no-progress", "--frozen")

RUNNING_ON_CI: bool = os.getenv("CI") is not None
RUNNING_ON_RTD: bool = os.getenv("READTHEDOCS") is not None


def _create_lockfile_pr_message(uv_output: str, session: nox.Session) -> None:
    """
    Create the pull request message during requirements updates.

    This function copies a GitHub flavored Markdown template to a new
    file and appends a table containing the updated requirements, with
    links to the corresponding PyPI pages. This file is then used as the
    body of the pull request message used in the workflow for updating
    requirements.

    ⚠️ This function requires that `uv.lock` existed before
    `uv lock --upgrade` was run.

    Parameters
    ----------
    uv_output : str
        The multi-line output of ``session.run(..., silent=True)``.
    """

    pr_template = pathlib.Path(
        ROOT_DIR / "./.github/content/upgrade-uv-lock-pr-template.md"
    )
    pr_message = pathlib.Path(ROOT_DIR / "./.github/content/upgrade-uv-lock-pr-body.md")

    shutil.copy(pr_template, pr_message)

    preamble = [
        "",
        "| package | old version | new version |",
        "| :-----: | :---------: | :---------: |",
    ]

    lines = []
    for line in uv_output.splitlines():
        if line.startswith("Resolved") or not line:
            continue

        if not line.startswith("Updated"):
            session.warn(f"Line not added to table: {line}")
            continue

        try:
            # An example line is "Updated nbsphinx v0.9.6 -> v0.9.7"
            _, package_name_, old_version_, _, new_version_ = line.strip().split()
        except ValueError:
            session.warn(f"Line not added to table: {line}:")
            continue

        old_version = f"{old_version_.removeprefix('v')}"
        new_version = f"{new_version_.removeprefix('v')}"

        pypi_link = f"https://pypi.org/project/{package_name_}/{new_version}"
        package_name = f"[`{package_name_}`]({pypi_link})"

        lines.append(f"| {package_name} | `{old_version}` | `{new_version}` |")

    table_of_package_upgrades = "\n".join(preamble + lines)
    with pr_message.open(mode="a") as file:
        file.write(table_of_package_upgrades)


@nox.session
def lock(session: nox.Session) -> None:
    """
    Upgrade Python environments used in CI, with a dependency cooldown.

    This session upgrades uv.lock: the cross-platform lockfile that
    defines the exact Python environments used when running tests,
    building documentation, and performing continuous integration (CI)
    checks.

    When run in CI, this session generates a file that contains the pull
    request message for the GitHub workflow that uses this session
    (:file:`.github/workflows/upgrade-uv-lock.yml`).
    """

    uv_lock = (
        "uv",
        "lock",
        "--upgrade",
        "--no-progress",
        *session.posargs,
    )
    try:
        # Use session.run() with silent=True to return the command output
        uv_output: str | bool = session.run(*uv_lock, silent=RUNNING_ON_CI)
    except nox.command.CommandFailed:
        session.warn("⚠️ uv.lock is invalid, likely due to a git merge conflict.")
        session.log(
            "📥 Checking out uv.lock from the branch being merged into this one."
        )
        session.log(
            "🪧 If this next attempt is unsuccessful, delete uv.lock and try again."
        )
        session.run("git", "checkout", "--theirs", "--", "uv.lock", external=True)
        uv_output: str | bool = session.run(*uv_lock, silent=RUNNING_ON_CI)

    if RUNNING_ON_CI:
        session.log(uv_output)
        _create_lockfile_pr_message(uv_output=uv_output, session=session)


@nox.session
def validate_lockfile(session: nox.Session) -> None:
    """
    Ensure that uv.lock is consistent with pyproject.toml.

    This check is normally performed locally when running pre-commit or
    prek. Because pre-commit.ci blocks network access, this check is
    instead done in CI via a GitHub workflow that calls this session.
    """
    if RUNNING_ON_CI:
        errmsg = (
            "The Python environments in file 'uv.lock' are inconsistent "
            "with the requirements defined in 'pyproject.toml'. "
            "After installing Nox, this problem can be fixed by running "
            "`nox -s validate_lockfile` in the top-level directory of "
            "your clone of PlasmaPy, and then pushing the updated "
            "'uv.lock' to GitHub. "
        )
    else:
        errmsg = (
            "File 'uv.lock' has been updated for consistency with the "
            "requirements defined in 'pyproject.toml'."
        )

    try:
        session.run("uv", "lock", "--no-progress")
    except nox.command.CommandFailed:
        session.error(errmsg)


# Define pytest flags that are only sometimes used. Define flags that
# are always used in tool.pytest.addopts in pyproject.toml instead of
# noxfile.py, since that allows users to run `pytest` and get

WITH_DOCTESTS: tuple[str, ...] = ("--doctest-modules", "--doctest-continue-on-failure")

WITH_COVERAGE: tuple[str, ...] = (
    "--cov=plasmapy",
    "--cov-report=xml",
    "--cov-config=pyproject.toml",
    "--cov-append",
    "--cov-report",
    "xml:coverage.xml",
)

REPORT_WARNINGS_ONCE = ("-W", "once")

SKIPSLOW: tuple[str, ...] = ("-m", "not slow")

test_specifiers: list[nox._parametrize.Param] = [
    nox.param("run all tests", id="all"),
    nox.param("skip slow tests", id="skipslow"),
    nox.param("with code coverage", id="cov"),
    nox.param("lowest-direct", id="lowest-direct"),
    nox.param("lowest-direct-skipslow", id="lowest-direct-skipslow"),
]


@nox.session(python=SUPPORTED_PYTHON_VERSIONS)
@nox.parametrize("test_specifier", test_specifiers)
def tests(session: nox.Session, test_specifier: nox._parametrize.Param) -> None:
    """Run tests with pytest."""

    options: list[str] = []

    if test_specifier in {"skip slow tests", "lowest-direct-skipslow"}:
        options += SKIPSLOW

    if test_specifier == "with code coverage":
        options += WITH_COVERAGE

    # There may be some warnings that got fixed later with the oldest
    # allowed versions of dependencies
    if test_specifier in {"lowest-direct", "lowest-direct-skipslow"}:
        options += REPORT_WARNINGS_ONCE

    # Doctests are only run with the most recent versions of Python and
    # other dependencies because there may be subtle differences in the
    # output between different versions of Python, NumPy, and Astropy.
    if session.python == MAXPYTHON and test_specifier not in {"lowest-direct", "cov"}:
        options += WITH_DOCTESTS

    if gh_token := os.getenv("GH_TOKEN"):
        session.env["GH_TOKEN"] = gh_token

    match test_specifier:
        case "lowest-direct" | "lowest-direct-skipslow":
            session.install(
                ".",
                "--resolution=lowest-direct",
                "--group=test",
            )
        case _:
            # From https://nox.thea.codes/en/stable/cookbook.html#using-a-lockfile
            # If we separate out the lowest-direct tests, then we can use
            # @nox_uv.session for tests too.
            session.run_install(
                *UV_SYNC,
                "--group=test",
                f"--python={session.virtualenv.location}",
                env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
            )

    session.run("pytest", *options, *session.posargs)


@nox_uv.session(python=MAXPYTHON, uv_groups=["test"])
@nox.parametrize(
    ["package"],
    [
        nox.param("numpy", id="numpy"),
        nox.param("astropy", id="astropy"),
        nox.param("pandas", id="pandas"),
        nox.param("xarray", id="xarray"),
        nox.param("https://github.com/lmfit/lmfit-py", id="lmfit"),
    ],
)
def test_upstream(session: nox.Session, package: str) -> None:
    """
    Run tests against the development branch of an upstream dependency.

    Testing against unreleased versions of upstream dependencies helps
    us catch problems before they make it into an official release.
    """
    if package.startswith("https"):
        session.install(f"git+{package}")
    else:
        session.run_install(
            "uv",
            "pip",
            "install",
            "--upgrade",
            package,
            env={
                "UV_INDEX": "https://pypi.anaconda.org/scientific-python-nightly-wheels/simple",
                "UV_INDEX_STRATEGY": "unsafe-best-match",
                "UV_PRERELEASE": "allow",
            },
        )
        session.run("uv", "pip", "show", package)

    session.run("pytest", *session.posargs)


if RUNNING_ON_RTD:
    rtd_output_path = pathlib.Path(os.environ.get("READTHEDOCS_OUTPUT")) / "html"
    rtd_output_path.mkdir(parents=True, exist_ok=True)
    doc_build_dir = str(rtd_output_path)
else:
    doc_build_dir = "docs/_build/html"

SPHINX_BASE_COMMAND: list[str] = [
    "sphinx-build",
    "docs/",
    doc_build_dir,
    "--nitpicky",
    "--quiet",
    "--keep-going",
]

if not RUNNING_ON_RTD:
    SPHINX_BASE_COMMAND.extend(["--fail-on-warning"])

BUILD_HTML: tuple[str, ...] = ("--builder", "html")
CHECK_HYPERLINKS: tuple[str, ...] = ("--builder", "linkcheck")

DOC_TROUBLESHOOTING_MESSAGE = """

📘 Tips for troubleshooting common documentation build failures are in
PlasmaPy's documentation guide at:

🔗 https://docs.plasmapy.org/en/latest/contributing/doc_guide.html#troubleshooting
"""


@nox_uv.session(python=DOCPYTHON, uv_groups=["docs"])
def docs(session: nox.Session) -> None:
    """
    Build documentation with Sphinx.

    This session may require installation of pandoc and graphviz.

    Configuration file: docs/conf.py
    """

    if RUNNING_ON_CI:
        session.log(DOC_TROUBLESHOOTING_MESSAGE)

    # Can we use pixi or conda to install graphviz and pandoc if they
    # are not installed?

    session.run_install("dot", "-V", external=True)
    session.run_install("pandoc", "--version", external=True)

    session.run(*SPHINX_BASE_COMMAND, *BUILD_HTML, *session.posargs)

    landing_page = pathlib.Path(doc_build_dir) / "index.html"
    if landing_page.exists():
        session.log(f"The documentation may be previewed at {landing_page}")
    else:
        session.error(f"Documentation preview landing page not found: {landing_page}")


@nox_uv.session(python=DOCPYTHON, uv_groups=["docs"])
def htmlzip(session: nox.Session) -> None:
    """Bundle documentation build into a zip file on Read the Docs."""

    if not RUNNING_ON_RTD:
        session.error("This session must be run on Read the Docs.")

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

    # mimic RTD default naming convention
    zip_name = f"{READTHEDOCS_PROJECT}-{READTHEDOCS_LANGUAGE}-{READTHEDOCS_VERSION}.zip"

    cwd = pathlib.Path.cwd()
    session.chdir(f"{READTHEDOCS_OUTPUT / 'htmlzip'}")
    session.run("zip", "-r", "-m", f"{zip_name}", ".", external=True)
    session.chdir(f"{cwd}")

    session.log(f"The htmlzip was placed in: {READTHEDOCS_OUTPUT / 'htmlzip'}")


@nox_uv.session(python=DOCPYTHON, uv_groups=["docs"])
@nox.parametrize(
    ["site", "repository"],
    [
        nox.param("github", "sphinx-doc/sphinx", id="sphinx"),
        nox.param("github", "readthedocs/sphinx_rtd_theme", id="sphinx_rtd_theme"),
        nox.param("github", "spatialaudio/nbsphinx", id="nbsphinx"),
        nox.param("github", "plasmapy/plasmapy_sphinx", id="plasmapy_sphinx"),
    ],
)
def docs_upstream(session: nox.Session, site: str, repository: str) -> None:
    """
    Build documentation against the development branch of an upstream dependency.

    The purpose of this session is to catch bugs and breaking changes
    so that they can be fixed or updated earlier rather than later.
    """
    session.install(f"git+https://{site}.com/{repository}")
    package = repository.split("/")[-1]
    session.run_install("uv", "pip", "show", package)
    session.run(*SPHINX_BASE_COMMAND, *BUILD_HTML, *session.posargs)


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


@nox_uv.session(python=DOCPYTHON, uv_groups="docs")
def linkcheck(session: nox.Session) -> None:
    """Check hyperlinks in documentation."""

    if RUNNING_ON_CI:
        session.log(LINKCHECK_TROUBLESHOOTING)

    session.run(*SPHINX_BASE_COMMAND, *CHECK_HYPERLINKS, *session.posargs)


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


@nox_uv.session(python=MAXPYTHON, uv_groups=["type_check"])
def mypy(session: nox.Session) -> None:
    """
    Perform static type checking.

    Configuration file: mypy.ini
    """
    if RUNNING_ON_CI:
        session.log(MYPY_TROUBLESHOOTING)

    session.run(
        "mypy",
        ".",
        "--install-types",
        "--non-interactive",
        "--show-error-context",
        "--show-error-code-links",
        "--pretty",
        *session.posargs,
    )


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

       nox -s 'changelog(final)' -- 2026.2.0
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

    original_file = pathlib.Path(ROOT_DIR / "CHANGELOG.rst")
    original_file.unlink()

    session.run(*towncrier, "--yes")

    destination = pathlib.Path(ROOT_DIR / f"docs/changelog/{version}.rst")
    shutil.copy(original_file, destination)


@nox_uv.session(python=MINPYTHON, uv_groups=["test"])
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
    session.install("autotyping", "typing_extensions")
    default_paths = ("src", "tests", "tools", "*.py", ".github", "docs/*.py")
    paths = session.posargs or default_paths
    session.run("python", "-m", "autotyping", *options, *paths, *session.posargs)


@nox.session
def manifest(session: nox.Session) -> None:
    """
    Check for missing files in MANIFEST.in.

    When run outside of CI, this check may report files that were
    locally created but not included in version control. These false
    positives can be ignored by adding file patterns and paths to
    `ignore` under `[tool.check-manifest]` in `pyproject.toml`.
    """
    # check-manifest would be suitable as a pre-commit hook, except that
    # it requires ∼10 seconds to build the package, which would triple
    # the time needed to run pre-commit.
    session.install("check-manifest")
    session.run("check-manifest", *session.posargs)


@nox.session
def lint(session: nox.Session) -> None:
    """
    Run all pre-commit hooks on all files with prek.

    Configuration file: .pre-commit-config.yaml
    """
    session.install("prek")
    session.run(
        "prek",
        "run",
        "--all-files",
        *session.posargs,
    )


ZIZMOR_TROUBLESHOOTING_MESSAGE = """

🪧 Run this check locally with `nox -s zizmor` to find potential
security vulnerabilities in GitHub workflows and perform safe fixes.

🧰 Perform safe and unsafe fixes with `nox -s zizmor -- --fix=all`.

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

    If no positional arguments are provided, safe fixes will be applied.
    To perform unsafe fixes, run `nox -s zizmor -- --fix=unsafe-only`.

    Configuration file: .github/zizmor.yml
    """
    session.log(ZIZMOR_TROUBLESHOOTING_MESSAGE)

    args = ["--no-progress", "--color=auto", *session.posargs]
    if not session.posargs:
        args.append("--fix=safe")

    session.install("zizmor")
    session.run(
        "zizmor",
        ".github",
        *args,
    )


if __name__ == "__main__":
    nox.main()
