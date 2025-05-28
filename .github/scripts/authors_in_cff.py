"""Check that the authors of a PR are included in ``CITATION.cff``."""

# /// script
# requires-python = ">=3.13"
# dependencies = ["requests>=2.32"]
# ///

import logging
import os
import pathlib
import sys

import requests

logging.basicConfig(level="INFO")

GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")
PR_NUMBER = os.getenv("PR_NUMBER")
REPO = os.getenv("GITHUB_REPOSITORY")

excluded_authors = {"dependabot[bot]", "pre-commit-ci[bot]", "sourcery-ai"}


def get_pr_authors() -> list[str]:
    """Get the GitHub usernames of the pull request."""
    url = f"https://api.github.com/repos/{REPO}/pulls/{PR_NUMBER}/commits"
    headers = {"Authorization": f"Bearer {GITHUB_TOKEN}"}
    response = requests.get(url, headers=headers, timeout=15)
    response.raise_for_status()

    authors = {
        commit["author"]["login"] for commit in response.json() if commit["author"]
    }

    return sorted(authors - excluded_authors)


def find_missing_github_usernames(authors: list[str]) -> list[str]:
    """Verify that all authors of a PR are included in :file:`CITATION.cff`."""
    with pathlib.Path("CITATION.cff").open() as file:
        lines = file.read()
        return [
            author
            for author in authors
            if f"alias: {author}" not in lines and f"- alias: {author}" not in lines
        ]


def main() -> None:
    """Check that all authors are included in CITATION.cff."""
    authors = get_pr_authors()
    missing_github_usernames = find_missing_github_usernames(authors)

    if not missing_github_usernames:
        msg = (
            f"The authors of pull request {PR_NUMBER} for {REPO} are: "
            f"{', '.join(sorted(authors))}. All authors are included "
            "in CITATION.cff. ✅️"
        )
        logging.info(msg)
        sys.exit(0)

    branch = os.getenv("GITHUB_HEAD_REF")
    username = missing_github_usernames[0]

    instructions_to_add_author = f"""

⚠️ To ensure that you get credit for your contribution, please add the
following authors to CITATION.cff: {", ".join(missing_github_usernames)!r}

This file can be edited at:

🔗️ https://github.com/{username}/PlasmaPy/edit/{branch}/CITATION.cff

Each entry should be of the form:

- given-names: <given names>
  family-names: <family names>
  affiliation: <affiliation>
  orcid: https://orcid.org/<ORCiD number>
  alias: {username}

All fields are optional except "alias", which is the GitHub username.

🏷️ We encourage all contributors to sign up for an ORCID iD: a unique,
persistent identifier used by researchers, authors, and open source
contributors. Sign up at: https://orcid.org/register

Thank you for contributing! 🌱️
"""

    logging.info(instructions_to_add_author)
    sys.exit(1)


if __name__ == "__main__":
    main()
