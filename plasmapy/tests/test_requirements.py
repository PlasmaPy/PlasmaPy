"""Tests for the consistency of requirements."""

from typing import List


def _read_requirements_txt_file(prefix: str, requirements_directory: str) -> List[str]:
    """
    Read in a text file containing requirements.

    Parameters
    ----------
    prefix : str
        The prefix to a filename, such as `"build"` for :file:`build.txt`.

    requirements_directory : str
        The path to the directory containing the requirements file.

    Returns
    -------
    list of str
        A `list` containing the lines of the requirements file,
        excluding lines that do not start with an alphabetic character
        (like comments).
    """
    filename = f"{requirements_directory}/{prefix}.txt"
    with open(filename) as file:
        lines_of_file = file.readlines()
    return [line.strip() for line in lines_of_file if line[0].isalpha()]
