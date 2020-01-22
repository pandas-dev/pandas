#!/usr/bin/env python
"""
Validate that the titles in the rst files follow the proper capitalization convention.

Print the titles that do not follow the convention.

Usage::
./scripts/validate_rst_title_capitalization.py doc/source/development/contributing.rst
./scripts/validate_rst_title_capitalization.py doc/source/

"""
import argparse
import sys
import re
import os
from typing import Tuple, Generator, List


CAPITALIZATION_EXCEPTIONS = {
    "pandas",
    "Python",
    "IPython",
    "PyTables",
    "Excel",
    "JSON",
    "HTML",
    "SAS",
    "SQL",
    "BigQuery",
    "STATA",
    "Interval",
    "PEP8",
    "Period",
    "Series",
    "Index",
    "DataFrame",
    "C",
    "Git",
    "GitHub",
    "NumPy",
    "Apache",
    "Arrow",
    "Parquet",
    "MultiIndex",
    "NumFOCUS",
    "sklearn",
    "Docker",
}

CAP_EXCEPTIONS_DICT = {word.lower(): word for word in CAPITALIZATION_EXCEPTIONS}

bad_title_dict = {}

err_msg = "Heading capitalization formatted incorrectly. Please correctly capitalize"


def correct_title_capitalization(title: str) -> str:
    """
    Algorithm to create the correct capitalization for a given title

    Parameters
    ----------
    title : str
        Heading string to correct

    Returns
    -------
    correct_title : str
        Correctly capitalized heading

    """

    correct_title: str = re.sub(r"^\W*", "", title).capitalize()

    removed_https_title = re.sub(r"<https?:\/\/.*[\r\n]*>", "", correct_title)

    word_list = re.split(r"\W", removed_https_title)

    for word in word_list:
        if word.lower() in CAP_EXCEPTIONS_DICT:
            correct_title = re.sub(
                r"\b" + word + r"\b", CAP_EXCEPTIONS_DICT[word.lower()], correct_title
            )

    return correct_title


def is_following_capitalization_convention(title: str) -> bool:
    """
    Function to return if a given title is capitalized correctly

    Parameters
    ----------
    title : str
        Heading string to validate

    Returns
    -------
    bool
        True if title capitalized correctly, False if not

    """

    correct_title = correct_title_capitalization(title)

    if title != correct_title:
        return False
    else:
        return True


def find_titles(rst_file: str) -> Generator[Tuple[str, int], None, None]:
    """
    Algorithm to identify particular text that should be considered headings in an
    RST file

    See <https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html> for details
    on what constitutes a string as a heading in RST

    Parameters
    ----------
    rst_file : str
        RST file to scan through for headings

    Yields
    -------
    title : str
        A heading found in the rst file

    line_number : int
        The corresponding line number of the heading

    """

    with open(rst_file, "r") as file_obj:
        lines = file_obj.read().split("\n")

    regex = {
        "*": r"^(?:\*{1})*$",
        "=": r"^(?:={1})*$",
        "-": r"^(?:-{1})*$",
        "^": r"^(?:\^{1})*$",
        "~": r"^(?:~{1})*$",
        "#": r"^(?:#{1})*$",
        '"': r'^(?:"{1})*$',
    }

    table = str.maketrans("", "", "*`_")

    for line_no in range(1, len(lines)):
        if len(lines[line_no]) != 0 and len(lines[line_no - 1]) != 0:
            for key in regex:
                match = re.search(regex[key], lines[line_no])
                if match is not None:
                    if line_no >= 2:
                        if lines[line_no] == lines[line_no - 2]:
                            if len(lines[line_no]) == len(lines[line_no - 1]):
                                yield lines[line_no - 1].translate(table), line_no
                            break
                    if len(lines[line_no]) >= len(lines[line_no - 1]):
                        yield lines[line_no - 1].translate(table), line_no


def fill_bad_title_dict(rst_file: str) -> None:
    """
    Method that fills up the bad_title_dict with incorrectly capitalized headings

    Parameters
    ----------
    rst_file : str
        Directory address of a .rst file as a string

    """

    if rst_file in bad_title_dict:
        return

    for title, line_number in find_titles(rst_file):
        if not is_following_capitalization_convention(title):
            if rst_file not in bad_title_dict:
                bad_title_dict[rst_file] = [(title, line_number)]
            else:
                bad_title_dict[rst_file].append((title, line_number))


def find_rst_files(source_paths: List[str]) -> Generator[str, None, None]:
    """
    Given the command line arguments of directory paths, this method
    yields the strings of the .rst file directories that these paths contain

    Parameters
    ----------
    source_paths : str
        List of directories to validate, provided through command line arguments

    Yields
    -------
    directory_address : str
        Directory address of a .rst files found in command line argument directories

    """

    for directory_address in source_paths:
        if not os.path.exists(directory_address):
            raise ValueError(
                "Please enter a valid path, pointing to a valid file/directory."
            )
        elif directory_address.endswith(".rst"):
            yield directory_address
        else:
            for (dirpath, _, filenames) in os.walk(directory_address):
                for file in filenames:
                    if file.endswith(".rst"):
                        yield os.path.join(dirpath, file)


def main(source_paths: List[str], output_format: str) -> bool:
    """
    The main method to print all headings with incorrect capitalization

    Parameters
    ----------
    source_paths : str
        List of directories to validate, provided through command line arguments
    output_format : str
        Output format of the script.

    Returns
    -------
    number_of_errors : int
        True if there are headings that are printed, False if not

    """

    number_of_errors: int = 0

    directory_list = find_rst_files(source_paths)

    for filename in directory_list:
        fill_bad_title_dict(filename)

    if len(bad_title_dict) == 0:
        return number_of_errors

    for key in bad_title_dict:
        for line in bad_title_dict[key]:
            print(
                f"""{key}:{line[1]}:{err_msg} "{line[0]}" to "{
                correct_title_capitalization(line[0])}" """
            )
            number_of_errors += 1

    return number_of_errors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate heading capitalization")

    parser.add_argument(
        "paths", nargs="+", default=".", help="Source paths of file/directory to check."
    )

    parser.add_argument(
        "--format",
        "-f",
        default="{source_path}:{line_number}:{msg}:{heading}",
        help="Output format of incorrectly capitalized titles",
    )

    args = parser.parse_args()

    sys.exit(main(args.paths, args.format))
