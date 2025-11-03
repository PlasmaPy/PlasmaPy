# /// script
# requires-python = ">=3.13"
# dependencies = [ "ruamel.yaml" ]
# ///

"""Functionality to update citation information prior to making a release."""

import argparse
import datetime
import pathlib
import re

from ruamel.yaml import YAML


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("citation_cff_file")
    parser.add_argument("citation_rst_file")
    parser.add_argument("whatsnew_index_rst_file")
    parser.add_argument("--version")
    parser.add_argument("--doi")
    parser.add_argument(
        "--date_released",
        type=datetime.date.fromisoformat,
        default=datetime.date.today(),  # noqa: DTZ011
        help="Date in iso format (e.g., 2025-12-26)",
    )

    return parser.parse_args()


yaml = YAML()


def update_citation_files(args) -> None:
    citation_cff_file = pathlib.Path(args.citation_cff_file)
    with citation_cff_file.open() as f:
        d = yaml.load(f)

    d["version"] = args.version
    d["identifiers"][0]["value"] = args.doi
    d["date-released"] = args.date_released.isoformat()

    with citation_cff_file.open("w") as f:
        yaml.dump(d, f)

    citation_rst_file = pathlib.Path(args.citation_rst_file)
    citation_rst_text = citation_rst_file.read_text()

    for source_regex, target_value in (
        (
            r"\|version_to_cite\| replace:: (.*)",
            f"|version_to_cite| replace:: {args.version}",
        ),
        (
            r"\|doi_hyperlink\| replace:: (.*)",
            f"|doi_hyperlink| replace:: https://doi.org/{args.doi}",
        ),
        (
            r"\|citation_year\| replace:: (.*)",
            f"|citation_year| replace:: {args.date_released.year}",
        ),
    ):
        citation_rst_text = re.compile(source_regex).sub(
            target_value, citation_rst_text
        )
    with citation_rst_file.open("w") as f:
        f.write(citation_rst_text)


if __name__ == "__main__":
    args = parse_arguments()
    update_citation_files(args)
