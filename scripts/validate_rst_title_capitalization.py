#!/usr/bin/env python

"""
GH #29641

Collect the titles in the rst files and validate if they follow the proper
capitalization convention.

Prints the titles that do not follow the convention.

Usage::
./scripts/validate_rst_title_capitalization.py doc/source/development/contributing.rst
./scripts/validate_rst_title_capitalization.py doc/source/

"""

import argparse
import sys
import re
import os
from os import walk
from typing import Generator, List


# Keynames that would not follow capitalization convention
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

# Lowercase representation of CAPITALIZATION_EXCEPTIONS
CAPITALIZATION_EXCEPTIONS_LOWER = {word.lower() for word in CAPITALIZATION_EXCEPTIONS}

# Dictionary of bad titles that will be printed later along with line numbers
# Key: Document Directory, Value: Pair(Bad Title, Line Number)
bad_title_dict = {}

# Error Message:
err_msg = "Heading capitalization formatted incorrectly. Please correctly capitalize"


def is_following_capitalization_convention(title: str) -> bool:
    """
    Algorithm to determine if a heading follows the capitalization convention

    This method returns true if the title follows the convention
    and false if it does not

    Parameters
    ----------
    title : str
        Heading string to validate

    Returns
    -------
    bool
        True if capitalization is correct, False if not

    """

    # Remove https link if present in heading
    title = re.sub(r"<https?:\/\/.*[\r\n]*>", "", title)

    # Split with delimiters comma, semicolon and space, parentheses, colon, slashes
    word_list = re.split(r"[;,-/():\s]\s*", title)

    # Edge Case: First word is an empty string
    if len(word_list[0]) == 0:
        return False

    # Dealing with the first word of the title
    if word_list[0] not in CAPITALIZATION_EXCEPTIONS:
        # word is not in CAPITALIZATION_EXCEPTIONS but has different capitalization
        if word_list[0].lower() in CAPITALIZATION_EXCEPTIONS_LOWER:
            return False
        # First letter of first word must be uppercase
        if not word_list[0][0].isupper():
            return False
        # Remaining letters of first word must not be uppercase
        for j in range(1, len(word_list[0])):
            if word_list[0][j].isupper():
                return False

    # Remaining letters must not be uppercase letters
    for i in range(1, len(word_list)):
        if word_list[i] not in CAPITALIZATION_EXCEPTIONS:
            # word is not in CAPITALIZATION_EXCEPTIONS but has different capitalization
            if word_list[i].lower() in CAPITALIZATION_EXCEPTIONS_LOWER:
                return False
            # Remaining letters must not be uppercase
            for j in range(len(word_list[i])):
                if word_list[i][j].isupper():
                    return False

    # Returning True if the heading follows the capitalization convention
    return True


def findTitles(rst_file: str) -> Generator[List[str], List[int], None]:
    """
    Algorithm to identify particular text that should be considered headings in an
    RST file

    See <https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html> for details
    on what constitutes a string as a heading in RST

    Parameters
    ----------
    rst_file : str
        RST file to scan through for headings

    Returns
    -------
    title_list : List[str]
        A list of heading strings found in the document tree

    line_number_list : List[int]
        The corresponding line numbers of the headings in title_list

    """

    # title_list is the list of headings that is encountered in the doctree
    title_list: List[str] = []

    # List of line numbers that corresponding headings in title_list can be found at
    line_number_list: List[int] = []

    # Open and read the .rst file and store the string of data into lines
    with open(rst_file, "r") as file_obj:
        lines = file_obj.read().split("\n")

    # Regular expressions that denote a title beforehand
    regex = {
        "*": r"^(?:\*{1})*$",
        "=": r"^(?:={1})*$",
        "-": r"^(?:-{1})*$",
        "^": r"^(?:\^{1})*$",
        "~": r"^(?:~{1})*$",
        "#": r"^(?:#{1})*$",
        '"': r'^(?:"{1})*$',
    }

    # '*`_' markers are removed from original string text.
    table = str.maketrans("", "", "*`_")

    # Loop through lines lines, appending if they are considered headings
    for lineno in range(1, len(lines)):
        if len(lines[lineno]) != 0 and len(lines[lineno - 1]) != 0:
            for key in regex:
                match = re.search(regex[key], lines[lineno])
                if match is not None:
                    if lineno >= 2:
                        if lines[lineno] == lines[lineno - 2]:
                            if len(lines[lineno]) == len(lines[lineno - 1]):
                                title_list.append(lines[lineno - 1].translate(table))
                                line_number_list.append(lineno)
                            break
                    if len(lines[lineno]) >= len(lines[lineno - 1]):
                        title_list.append(lines[lineno - 1].translate(table))
                        line_number_list.append(lineno)

    return title_list, line_number_list


def fill_bad_title_dict(rst_file: str) -> None:
    """
    Method that fills up the bad_title_dict with incorrectly capitalized headings

    Parameters
    ----------
    rst_file : str
        Directory address of a .rst file as a string

    """

    # Ensure this file doesn't already have a bad_title_dict slot
    if rst_file in bad_title_dict:
        return

    # Make a list of headings along with their line numbers
    title_list, line_number_list = findTitles(rst_file)

    # Append the bad_title_dict if the capitalization convention not followed
    for i in range(len(title_list)):
        if not is_following_capitalization_convention(title_list[i]):
            if rst_file not in bad_title_dict:
                bad_title_dict[rst_file] = [(title_list[i], line_number_list[i])]
            else:
                bad_title_dict[rst_file].append((title_list[i], line_number_list[i]))


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

    # Loop through source_paths, recursively looking for .rst files
    for directory_address in source_paths:
        if not os.path.exists(directory_address):
            raise ValueError(
                "Please enter a valid path, pointing to a valid file/directory."
            )
        elif directory_address.endswith(".rst"):
            yield directory_address
        else:
            for (dirpath, dirnames, filenames) in walk(directory_address):
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
    is_failed : bool
        True if there are headings that are printed, False if not

    """

    is_failed: bool = False

    # Make a list of all RST files from command line directory list
    directory_list = find_rst_files(source_paths)

    # Fill the bad_title_dict, which contains all incorrectly capitalized headings
    for filename in directory_list:
        fill_bad_title_dict(filename)

    # Return an exit status of 0 if there are no bad titles in the dictionary
    if len(bad_title_dict) == 0:
        return is_failed

    # Print bad_title_dict Results
    is_failed = True
    for key in bad_title_dict:
        for line in bad_title_dict[key]:
            print(key + ":" + str(line[1]) + ": " + err_msg + ' "' + line[0] + '"')

    # Exit status of 0
    return is_failed


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate heading capitalization")

    parser.add_argument(
        "paths", nargs="+", default=".", help="Source paths of file/directory to check."
    )

    parser.add_argument(
        "--format",
        "-f",
        default="{source_path}:{line_number}:{heading}:{msg}",
        help="Output format of incorrectly capitalized titles",
    )

    args = parser.parse_args()

    sys.exit(main(args.paths, args.format))
