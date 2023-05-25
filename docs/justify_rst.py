import click
import re
import textwrap

from pathlib import Path
from typing import Set


def should_preserve(line: str) -> bool:
    """
    Check whether a line should be preserved without formatting.

    Parameters
    ----------
    line : str
        Line of text.

    Returns
    -------
    bool
        True if the line should be preserved, False otherwise.
    """
    return line == "" or re.match(r"^\s|^\W(?![\(\[<`])|^\d+\.\s|^\>\>\>", line)


def extract_words(line: str) -> Set[str]:
    """
    Extract all words from a line.

    Parameters
    ----------
    line : str
        Line of text.

    Returns
    -------
    set of str
        Set of all words.
    """
    return set(line.split())


def justify_paragraph(paragraph: str) -> str:
    """
    Justify a paragraph while preserving word count.

    Parameters
    ----------
    paragraph : str
        Paragraph text.

    Returns
    -------
    str
        Justified paragraph, or original paragraph if word count would change.
    """
    justified_paragraph = textwrap.fill(paragraph, width=72)
    justified_paragraph = justified_paragraph.replace(
        ".  ", ". "
    )  # Reduce two spaces after a period to one space
    if extract_words(paragraph) != extract_words(justified_paragraph):
        return paragraph
    return justified_paragraph


def justify_rst_file(filepath: Path) -> None:
    """
    Justify the text in a reStructuredText file.

    Parameters
    ----------
    filepath : pathlib.Path
        Path to the file.

    Returns
    -------
    None
    """
    with filepath.open("r") as f:
        lines = f.readlines()

    justified_text = []
    paragraph = ""
    for line in lines:
        if should_preserve(line):
            if paragraph:
                justified_text.append(justify_paragraph(paragraph).rstrip())
                paragraph = ""
            justified_text.append(
                line.rstrip("\n").rstrip()
            )  # Strip trailing whitespace
        else:
            paragraph += line.strip() + " "

    # Handle last paragraph
    if paragraph:
        justified_text.append(justify_paragraph(paragraph.strip()).rstrip())

    # Ensure file ends with a newline
    if justified_text and not justified_text[-1].endswith("\n"):
        justified_text[-1] += "\n"

    # Overwrite the file with the justified text
    with filepath.open("w") as f:
        f.write("\n".join(justified_text))


@click.command()
@click.argument("dirpath", type=click.Path(exists=True, file_okay=False), default=".")
def main(dirpath: str) -> None:
    """
    Main function. Justify all .rst files in a directory and its subdirectories.

    Parameters
    ----------
    dirpath : str
        Path to the directory.

    Returns
    -------
    None
    """
    # Walk through the directory tree
    for filepath in Path(dirpath).rglob("*.rst"):
        justify_rst_file(filepath)


if __name__ == "__main__":
    main()
