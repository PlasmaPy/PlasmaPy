"""
Convert author information from :file:`CITATION.cff` into a
reStructuredText-formatted list.
"""
import pathlib
import yaml

from typing import Union


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

    Raises
    ------
    `yaml.YAMLError`
        If :file:`CITATION.cff` cannot be parsed.
    """
    with pathlib.Path(filename).open() as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as e:
            print(e)


def sorting_key(author: dict[str, str]) -> str:
    """
    Provide sorting key for an author.

    Parameters
    ----------
    author : `dict`
        A dictionary containing author information.

    Returns
    -------
    str
        The key for sorting authors.
    """
    if "family-names" in author:
        return author["family-names"]
    elif "given-names" in author:
        return author["given-names"]
    elif "alias" in author:
        return author["alias"]
    else:
        return ""


def generate_reST(authors: list[dict[str, str]]) -> str:
    """
    Generate a reStructuredText formatted list of authors.

    Parameters
    ----------
    authors : list of dict
        A list of dictionaries, each containing information about an
        author.

    Returns
    -------
    str
        The reStructuredText formatted list of authors.
    """
    authors = sorted(authors, key=sorting_key)
    authors_reST = ""
    for author in authors:
        if "given-names" in author and "family-names" in author:
            name = f'{author["given-names"]} {author["family-names"]}'
        elif "alias" in author:
            name = f'{author["alias"]}'
        else:
            continue

        if "orcid" in author and "alias" in author:
            orcid = author["orcid"].removeprefix("https://orcid.org/")
            authors_reST += f'- :user:`{name} <{author["alias"]}>` (:orcid:`{orcid}`)\n'
        elif "alias" in author:
            authors_reST += f'- :user:`{name} <{author["alias"]}>`\n'
        else:
            authors_reST += f"- :user:`{name}`\n"
    return authors_reST


def main(cff_file="../CITATION.cff", rst_file="about/_authors.rst"):
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
    authors_reST = generate_reST(cff_data["authors"])
    with pathlib.Path(rst_file).open("w") as file:
        file.write(authors_reST)
