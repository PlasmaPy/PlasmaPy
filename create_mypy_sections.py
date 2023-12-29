#!/usr/bin/env python3

import os
import pathlib


def read_mypy_output_file() -> list[str]:
    mypy_output = pathlib.Path("./mypy_output")
    with mypy_output.open() as file:
        return file.readlines()


def extract_error_lines(lines: list[str]) -> list[str]:
    return [line.strip() for line in lines if "error:" in line]


def get_filename(line: str) -> str:
    return line.split(":")[0]


def get_error_code(line: str) -> str:
    return (line.split("[")[-1]).removesuffix("]")


def convert_filename_to_module(filename: str) -> str:
    return filename.removesuffix(".py").replace("/", ".")


def get_filename_error_code_pairs(lines: list[str]) -> list[tuple[str, str]]:
    return [(get_filename(line), get_error_code(line)) for line in lines]


def get_error_codes_per_file(
    filename_error_code_pairs: list[tuple[str, str]],
) -> dict[str, set[str]]:
    error_codes_per_file: dict[str, set[str]] = {
        filename: set() for filename, error_code in filename_error_code_pairs
    }

    for filename, error_code in filename_error_code_pairs:
        error_codes_per_file[filename].add(error_code)

    return error_codes_per_file


def filename_to_mypy_ini_section_string(filename: str) -> str:
    modulename = filename.removesuffix(".py").replace("/", ".").replace(".__init__", "")
    return f"[mypy-{modulename}]"


def print_mypy_ini_sections(error_codes_per_file: dict[str, set[str]]) -> None:
    for filename in error_codes_per_file:
        error_codes = sorted(error_codes_per_file[filename])
        print(filename_to_mypy_ini_section_string(filename))
        print(f'disable_error_code = {",".join(error_codes)}')
        print()


def main() -> None:
    os.system("mypy --no-pretty . > mypy_output")

    output_file_lines = read_mypy_output_file()
    error_lines = extract_error_lines(output_file_lines)
    filename_error_code_pairs: list[tuple[str, str]] = get_filename_error_code_pairs(
        error_lines
    )
    error_codes_per_file: dict[str, set[str]] = get_error_codes_per_file(
        filename_error_code_pairs
    )
    print_mypy_ini_sections(error_codes_per_file)

    os.system("rm -f mypy_output")


if __name__ == "__main__":
    main()
