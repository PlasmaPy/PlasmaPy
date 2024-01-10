"""Check that the authors of a PR are included in ``CITATION.cff``."""

import os
import pathlib
import sys
from typing import Optional

import requests

GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")
PR_NUMBER = os.getenv("PR_NUMBER")
REPO = os.getenv("GITHUB_REPOSITORY")

EXCLUDED_USERS = ["dependabot[bot]", "pre-commit-ci[bot]", "sourcery-ai"]


def get_pr_authors() -> set[str]:
    """Get the GitHub usernames of a pull request."""
    url = f"https://api.github.com/repos/{REPO}/pulls/{PR_NUMBER}/commits"
    headers = {"Authorization": f"Bearer {GITHUB_TOKEN}"}
    response = requests.get(url, headers=headers, timeout=15)
    response.raise_for_status()

    authors = {
        commit["author"]["login"] for commit in response.json() if commit["author"]
    }
    return authors - set(EXCLUDED_USERS)


def check_citation_file(authors: set[str]) -> tuple[bool, Optional[str]]:
    """Verify that all authors of a PR are included in :file:`CITATION.cff`."""
    with pathlib.Path("CITATION.cff").open() as file:
        contents = file.read()
        for author in authors:
            if f"alias: {author}" not in contents:
                return False, author
    return True, None


def main():
    """Check that all authors are included in CITATION.cff."""
    authors = get_pr_authors()
    check_passed, missing_github_username = check_citation_file(authors)
    if not check_passed:
        branch_name = os.getenv("GITHUB_HEAD_REF")

        error_message = f"""
To ensure that you get credit for your contribution to PlasmaPy, please
add {missing_github_username!r} as an author to CITATION.cff.

The entry should be of the form:

- given-names: <given names>
  family-names: <family names>
  affiliation: <affiliation>
  orcid: https://orcid.org/<ORCiD number>
  alias: {missing_github_username}
  email: <email address>

This file can be edited directly on GitHub at:

https://github.com/{REPO}/edit/{branch_name}/CITATION.cff

We encourage all contributors to sign up for an ORCID iD: a unique,
persistent identifier used by researchers, authors, and open source
contributors. Sign up at: https://orcid.org/register

All fields are optional except "alias", which is the GitHub username.
The "affiliation", "orcid", and/or "email" fields are sometimes needed
for conference abstract or journal article submissions about PlasmaPy.

Thank you for contributing!
"""

        print(error_message)  # noqa: T201
        sys.exit(1)


if __name__ == "__main__":
    main()
