"""Check that the author(s) of a PR are included in ``CITATION.cff``."""

import os
import pathlib
import requests
import sys

from typing import NoReturn, Optional

GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")
PR_NUMBER = os.getenv("PR_NUMBER")
REPO = os.getenv("GITHUB_REPOSITORY")

EXCLUDED_USERS = ["dependabot", "pre-commit-ci", "sourcery-ai"]


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


def comment_exists(message: str) -> bool:
    """Check if a comment already exists in the pull request."""
    url = f"https://api.github.com/repos/{REPO}/issues/{PR_NUMBER}/comments"
    headers = {"Authorization": f"Bearer {GITHUB_TOKEN}"}
    response = requests.get(url, headers=headers, timeout=15)
    response.raise_for_status()

    return any(message in comment["body"] for comment in response.json())


def post_comment(message: str) -> NoReturn:
    """Post a comment on a pull request."""
    marker = "To ensure that you get credit for your contribution to PlasmaPy"

    if not comment_exists(marker):
        url = f"https://api.github.com/repos/{REPO}/issues/{PR_NUMBER}/comments"
        headers = {"Authorization": f"Bearer {GITHUB_TOKEN}"}
        response = requests.post(
            url, json={"body": message}, headers=headers, timeout=15
        )
        response.raise_for_status()


def main():
    """Check that all authors are included in CITATION.cff."""
    authors = get_pr_authors()
    check_passed, failed_author = check_citation_file(authors)
    if not check_passed:
        branch_name = os.getenv("GITHUB_HEAD_REF")

        comment = (
            f"To ensure that you get credit for your contribution to "
            f"PlasmaPy, please add yourself as an author to:"
            f"[CITATION.cff](https://github.com/{REPO}/edit/{branch_name}/CITATION.cff)."
            f"The entry should be of the form:\n\n```"
            f"- given-names: <add given names>\n"
            f"  family-names: <family names>\n"
            f"  affiliation: <affiliation>\n"
            f"  orcid: https://orcid.org/<ORCiD ID number>\n"
            f"  alias: <GitHub username>\n```\n"
            f"All fields except `alias` are optional. "
            f"Optionally, [sign up for ORCiD](https://orcid.org/register)."
        )

        post_comment(comment)

        message = (
            f"\nTo ensure that you get credit for your contribution to "
            "PlasmaPy, please add yourself as an author to CITATION.cff "
            f"at: https://github.com/{REPO}/edit/{branch_name}/CITATION.cff\n\n"
            "The entry should be of the form:\n\n"
            "- given-names: <given names>\n"
            "  family-names: <family names>\n"
            "  affiliation: <affiliation>\n"
            "  orcid: https://orcid.org/<ORCiD number>\n"
            "  alias: <GitHub username>\n\n"
            "All fields except `alias` are optional. "
            "Optionally, sign up for ORCiD: https://orcid.org/register\n"
        )

        print(message)

        sys.exit(1)


if __name__ == "__main__":
    main()
