#!/usr/bin/env python
"""
Validate that the titles in the rst files follow the proper capitalization convention.

Print the titles that do not follow the convention.

Usage::
./scripts/validate_rst_title_capitalization.py doc/source/development/contributing.rst
./scripts/validate_rst_title_capitalization.py doc/source/

"""
import argparse
import glob
import os
import re
import sys
from typing import Iterable, List, Tuple

CAPITALIZATION_EXCEPTIONS = {
    "pandas",
    "pd",
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
    "DataFrames",
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
    "PeriodIndex",
    "NA",
    "NaN",
    "NaT",
    "ValueError",
    "Boolean",
    "BooleanArray",
    "KeyError",
    "API",
    "FAQ",
    "IO",
    "Timedelta",
    "TimedeltaIndex",
    "DatetimeIndex",
    "IntervalIndex",
    "CategoricalIndex",
    ("categorical", "Categorical"),
    "GroupBy",
    "SPSS",
    "ORC",
    "R",
    "HDF5",
    "HDFStore",
    "CDay",
    "CBMonthBegin",
    "CBMonthEnd",
    "BMonthBegin",
    "BMonthEnd",
    "BDay",
    "FY5253Quarter",
    "FY5253",
    "YearBegin",
    "YearEnd",
    "BYearBegin",
    "BYearEnd",
    "YearOffset",
    "QuarterBegin",
    "QuarterEnd",
    "BQuarterBegin",
    "BQuarterEnd",
    "QuarterOffset",
    "LastWeekOfMonth",
    "WeekOfMonth",
    "SemiMonthBegin",
    "SemiMonthEnd",
    "SemiMonthOffset",
    "CustomBusinessMonthBegin",
    "CustomBusinessMonthEnd",
    "BusinessMonthBegin",
    "BusinessMonthEnd",
    "MonthBegin",
    "MonthEnd",
    "MonthOffset",
    "CustomBusinessHour",
    "CustomBusinessDay",
    "BusinessHour",
    "BusinessDay",
    "DateOffset",
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December",
    "Float64Index",
    "FloatIndex",
    "TZ",
    "GIL",
    "strftime",
    "XPORT",
    "Unicode",
    "East",
    "Asian",
    "None",
    "URLs",
    "UInt64",
    "SciPy",
    "Matplotlib",
    "PyPy",
    "SparseDataFrame",
    "Google",
    "CategoricalDtype",
    "UTC",
    "Panel",
}

CAP_EXCEPTIONS_DICT = {
    exception.lower() if type(exception) is str else exception[0].lower(): exception
    for exception in CAPITALIZATION_EXCEPTIONS
}


err_msg = "Heading capitalization formatted incorrectly. Please correctly capitalize"

symbols = ("*", "=", "-", "^", "~", "#", '"')


def correct_title_capitalization(title: str) -> (str, list):
    """
    Algorithm to create the correct capitalization for a given title.

    Parameters
    ----------
    title : str
        Heading string to correct.

    Returns
    -------
    str
        Correctly capitalized heading.
    """

    # Skip modification no matter what if title begins by ":" to exclude specific
    # syntax that is needed to build links.
    if title[0] == ":":
        return title

    # Keep record of current words capitalization before modification
    copy_title_list: list = [
        el for el in title.split(" ") if el.lower() in CAP_EXCEPTIONS_DICT
    ]

    # Strip all non-word characters from the beginning of the title to the
    # first word character.
    correct_title: str = re.sub(r"^\W*", "", title).capitalize()

    # Remove a URL from the title. We do this because words in a URL must
    # stay lowercase, even if they are a capitalization exception.
    removed_https_title = re.sub(r"<https?:\/\/.*[\r\n]*>", "", correct_title)

    # Split a title into a list using non-word character delimiters.
    word_list = re.split(r"\W", removed_https_title)

    tuple_exceptions_list: list = []
    for word in word_list:
        if word:
            for exception in CAP_EXCEPTIONS_DICT:

                # If type of exception is a string, and word is equal to that exception,
                # correct with the right syntax
                if type(CAP_EXCEPTIONS_DICT[exception]) is str:
                    if (
                        word.lower() == exception
                        and word != CAP_EXCEPTIONS_DICT[exception]
                    ):
                        correct_title = re.sub(
                            rf"\b{word}\b",
                            CAP_EXCEPTIONS_DICT[exception],
                            correct_title,
                        )

                # Else if type of exception is a tuple and word is equal to that exception,
                # apply particular algorithm to check if word is in one of the element of the
                # tuple, if so correct with that syntax, else correct with the first element
                # syntax
                elif type(CAP_EXCEPTIONS_DICT[exception]) is tuple:
                    if word.lower() == exception:
                        for el in copy_title_list:
                            if word.lower() == el.lower():
                                if el not in CAP_EXCEPTIONS_DICT[exception]:
                                    correct_title = re.sub(
                                        rf"\b{word}\b",
                                        CAP_EXCEPTIONS_DICT[exception][0],
                                        correct_title,
                                    )
                                elif el in CAP_EXCEPTIONS_DICT[exception]:
                                    correct_title = re.sub(
                                        rf"\b{word}\b", el, correct_title
                                    )

                    # If type of exception is tuple, add that tuple to a list
                    if CAP_EXCEPTIONS_DICT[exception] not in tuple_exceptions_list:
                        tuple_exceptions_list.append(CAP_EXCEPTIONS_DICT[exception])

    # Return the correct title and the list corresponding to tuple exceptions
    # that match a word in title
    return correct_title, tuple_exceptions_list


def find_titles(rst_file: str) -> Iterable[Tuple[str, int]]:
    """
    Algorithm to identify particular text that should be considered headings in an
    RST file.

    See <https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html> for details
    on what constitutes a string as a heading in RST.

    Parameters
    ----------
    rst_file : str
        RST file to scan through for headings.

    Yields
    -------
    title : str
        A heading found in the rst file.

    line_number : int
        The corresponding line number of the heading.
    """

    with open(rst_file, "r") as fd:
        previous_line = ""
        for i, line in enumerate(fd):
            line = line[:-1]
            line_chars = set(line)
            if (
                len(line_chars) == 1
                and line_chars.pop() in symbols
                and len(line) == len(previous_line)
            ):
                yield re.sub(r"[`\*_]", "", previous_line), i
            previous_line = line


def find_rst_files(source_paths: List[str]) -> Iterable[str]:
    """
    Given the command line arguments of directory paths, this method
    yields the strings of the .rst file directories that these paths contain.

    Parameters
    ----------
    source_paths : str
        List of directories to validate, provided through command line arguments.

    Yields
    -------
    str
        Directory address of a .rst files found in command line argument directories.
    """

    for directory_address in source_paths:
        if not os.path.exists(directory_address):
            raise ValueError(
                "Please enter a valid path, pointing to a valid file/directory."
            )
        elif directory_address.endswith(".rst"):
            yield directory_address
        else:
            for filename in glob.glob(
                pathname=f"{directory_address}/**/*.rst", recursive=True
            ):
                yield filename


def main(source_paths: List[str], output_format: str) -> int:
    """
    The main method to print all headings with incorrect capitalization.

    Parameters
    ----------
    source_paths : str
        List of directories to validate, provided through command line arguments.
    output_format : str
        Output format of the script.

    Returns
    -------
    int
        Number of incorrect headings found overall.
    """

    number_of_errors: int = 0

    for filename in find_rst_files(source_paths):
        for title, line_number in find_titles(filename):
            if title != correct_title_capitalization(title)[0]:
                print(
                    f"""{filename}:{line_number}:{err_msg} "{title}" to "{
                    correct_title_capitalization(title)[0]}" """
                )
                number_of_errors += 1

                # Check for different wording
                if correct_title_capitalization(title)[1]:
                    for el in correct_title_capitalization(title)[1]:
                        for word in title.split(" "):
                            if word in el:
                                print(
                                    f"\n\nBe careful. In '{title}', '{word}' can have different writings:",
                                    end=" ",
                                )
                                for elem in el:
                                    print(elem, end=" ")
                                print("\n")

            # Check for different wording
            elif correct_title_capitalization(title)[1]:
                for el in correct_title_capitalization(title)[1]:
                    for word in title.split(" "):
                        if word.lower() in el:
                            print(
                                f"\n\nBe careful. In '{title}', '{word}' can have different writings:",
                                end=" ",
                            )
                            for elem in el:
                                print(elem, end=" ")
                            print("\n")

    return number_of_errors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate heading capitalization")

    parser.add_argument(
        "paths", nargs="+", default=".", help="Source paths of file/directory to check."
    )

    parser.add_argument(
        "--format",
        "-f",
        default="{source_path}:{line_number}:{msg}:{heading}:{correct_heading}",
        help="Output format of incorrectly capitalized titles",
    )

    args = parser.parse_args()

    sys.exit(main(args.paths, args.format))
