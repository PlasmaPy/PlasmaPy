"""Check that the authors of a PR are included in ``CITATION.cff``."""

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


def get_pr_authors() -> set[str]:
    """Get the GitHub usernames of the pull request."""
    url = f"https://api.github.com/repos/{REPO}/pulls/{PR_NUMBER}/commits"
    headers = {"Authorization": f"Bearer {GITHUB_TOKEN}"}
    response = requests.get(url, headers=headers, timeout=15)
    response.raise_for_status()

    authors = {
        commit["author"]["login"] for commit in response.json() if commit["author"]
    }

    return authors - excluded_authors


def find_missing_github_usernames(authors: set[str]) -> set[str]:
    """Verify that all authors of a PR are included in :file:`CITATION.cff`."""
    with pathlib.Path("CITATION.cff").open() as file:
        lines = file.read()
        return {author for author in authors if f"alias: {author}" not in lines}


def main():
    """Check that all authors are included in CITATION.cff."""
    authors = get_pr_authors()
    missing_github_usernames = find_missing_github_usernames(authors)

    if not missing_github_usernames:
        msg = (
            f"The authors of pull request {PR_NUMBER} for {REPO} are: "
            f"{', '.join(sorted(authors))}. No authors need to be "
            "added to CITATION.cff."
        )
        logging.info(msg)
        sys.exit(0)

    branch_name = os.getenv("GITHUB_HEAD_REF")

    error_message = f"""
To ensure that you get credit for your contribution to PlasmaPy, please
add the following authors to CITATION.cff: {missing_github_usernames!r}

The entry should be of the form:

- given-names: <given names>
  family-names: <family names>
  affiliation: <affiliation>
  orcid: https://orcid.org/<ORCiD number>
  alias: <GitHub username>
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

    logging.info(error_message)
    sys.exit(1)


if __name__ == "__main__":
    main()
