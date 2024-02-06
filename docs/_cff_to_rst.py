"""
Convert author information from :file:`CITATION.cff` into a
reStructuredText-formatted list.
"""
import pathlib
from typing import Union

import yaml
from unidecode import unidecode

# If a contributor changed their contributor name and we don't know
# the new name to update in in CITATION.cff, add their old username to
# the following set to prevent generating a broken link that would get
# detected when doing a `make linkcheck`.

obsolete_github_usernames = {}


def parse_cff(filename: str) -> dict[str, Union[str, list[dict[str, str]]]]:
    """
    Parse a :file:`CITATION.cff` file into a dictionary.

    Parameters
    ----------
    filename : `str`
        The path to :file:`CITATION.cff`.

    Returns
    -------
    `dict`
        A dictionary containing the parsed data from :file:`CITATION.cff`.
    """
    with pathlib.Path(filename).open() as stream:
        return yaml.safe_load(stream)


def author_sorting_key(author: dict[str, str]) -> str:
    """
    Provide sorting key for an author.

    Parameters
    ----------
    author : `dict` of `str` to `str`
        A dictionary containing author information, which should be an
        item in the list returned by `parse_cff`.

    Returns
    -------
    str
        The key for sorting authors.
    """
    if "family-names" in author:
        key = author["family-names"].lower()
    elif "given-names" in author:
        key = author["given-names"].lower()
    elif "alias" in author:
        key = author["alias"].lower()
    else:
        msg = (
            "Need to specify 'family-names', 'given-names', and/or "
            f"'alias' in CITATION.cff for author: {author}"
        )
        raise RuntimeError(msg)

    return unidecode(key)


def get_author_name(author: dict[str, str]) -> str:
    """
    Return the full name of the author if available, and otherwise
    return the author's alias.

    Parameters
    ----------
    author : `dict` of `str` to `str`
        A dictionary containing author information, which should be an
        item in the list returned by `parse_cff`.

    Returns
    -------
    str
        The author name (if available) or the author alias.
    """

    given_names = author.get("given-names", "")
    family_names = author.get("family-names", "")
    alias = author.get("alias", "")

    if given_names or family_names:
        return f"{given_names} {family_names}".strip()

    if alias:
        return alias

    msg = (
        "Invalid author information. Need at least one of "
        "'alias', 'family-names', or 'given-names' in "
        "CITATION.cff. Author information: {author}"
    )

    raise RuntimeError(msg)


def begin_author_line(author: dict[str, str]) -> str:
    """
    Return the beginning of the author line, using the ``:user:`` role
    if the alias is available.

    Parameters
    ----------
    author : `dict` of `str` to `str`
        A dictionary containing author information, which should be an
        item in the list returned by `parse_cff`.

    Returns
    -------
    str
        The beginning of the author line.
    """
    name = get_author_name(author)
    alias = author.get("alias")

    if not alias or alias in obsolete_github_usernames:
        return f"- {name}"

    if alias == name:
        return f"- :user:`{alias}`"

    return f"- :user:`{name} <{alias}>`"


def get_orcid(author: dict[str, str]) -> str:
    """
    Return the string that returns a string containing ORCID information
    using the ``:orcid:`` role.

    Parameters
    ----------
    author : `dict` of `str` to `str`
        A dictionary containing author information, which should be an
        item in the list returned by `parse_cff`.

    Returns
    -------
    str
        The ORCID information using the ``:orcid:`` role.
    """
    orcid = author.get("orcid", "").removeprefix("https://orcid.org/")
    return f" (:orcid:`{orcid}`)" if orcid else ""


def generate_rst_author_list(authors: list[dict[str, str]]) -> str:
    """
    Generate a reStructuredText formatted list of authors.

    Parameters
    ----------
    authors : `list` of `dict` from `str` to `str`
        A list of dictionaries containing information about an author,
        as read in by `parse_cff`.

    Returns
    -------
    str
        The reStructuredText formatted list of authors.
    """
    authors = sorted(authors, key=author_sorting_key)
    author_lines = [begin_author_line(author) + get_orcid(author) for author in authors]
    return "\n".join(author_lines)


def main(cff_file="../CITATION.cff", rst_file="about/_authors.rst", verbose=False):
    """
    Parse :file:`CITATION.cff` file and generate a reStructuredText
    formatted list of authors.

    Parameters
    ----------
    cff_file : str, default: ``"../CITATION.cff"``
        The path to :file:`CITATION.cff`.

    rst_file : str, default: ``about/_authors.rst``
        The path to the output :file:`.rst` file.
    """
    cff_data = parse_cff(cff_file)
    authors_rst = generate_rst_author_list(cff_data["authors"])

    if verbose:
        print(authors_rst)  # noqa: T201

    with pathlib.Path(rst_file).open("w") as file:
        file.write(authors_rst)


if __name__ == "__main__":
    """
    To test the functionality in this file, run:

    .. code-block: bash

        python _cff_to_rst.py
    """
    main(verbose=True)
