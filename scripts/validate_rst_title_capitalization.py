"""
Validate that the titles in the rst files follow the proper capitalization convention.

Print the titles that do not follow the convention.

Usage::

As pre-commit hook (recommended):
    pre-commit run title-capitalization --all-files

From the command-line:
    python scripts/validate_rst_title_capitalization.py <rst file>
"""
from __future__ import annotations

import argparse
import re
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable

CAPITALIZATION_EXCLUSIONS = {
    "ADBC",
    "API",
    "Apache",
    "April",
    "Arrow",
    "Asian",
    "August",
    "BDay",
    "BMonthBegin",
    "BMonthEnd",
    "BQuarterBegin",
    "BQuarterEnd",
    "BYearBegin",
    "BYearEnd",
    "BigQuery",
    "Boolean",
    "BooleanArray",
    "BusinessDay",
    "BusinessHour",
    "BusinessMonthBegin",
    "BusinessMonthEnd",
    "C",
    "CBMonthBegin",
    "CBMonthEnd",
    "CDay",
    "Categorical",
    "CategoricalDtype",
    "CategoricalIndex",
    "Copy-on-Write",
    "CustomBusinessDay",
    "CustomBusinessHour",
    "CustomBusinessMonthBegin",
    "CustomBusinessMonthEnd",
    "DataFrame",
    "DataFrameGroupBy",
    "DataFrames",
    "DateOffset",
    "DatetimeIndex",
    "December",
    "Docker",
    "East",
    "Excel",
    "ExtensionArray",
    "FAQ",
    "FY5253",
    "FY5253Quarter",
    "False",
    "February",
    "Float64Index",
    "FloatIndex",
    "GIL",
    "Git",
    "GitHub",
    "Gitpod",
    "Google",
    "GroupBy",
    "HDF5",
    "HDFStore",
    "HTML",
    "I",
    "IO",
    "I/O",
    "IPython",
    "Index",
    "Interval",
    "IntervalArray",
    "IntervalIndex",
    "JSON",
    "January",
    "July",
    "June",
    "KeyError",
    "LZMA",
    "LastWeekOfMonth",
    "Liveserve",
    "M",
    "ME",
    "March",
    "Matplotlib",
    "May",
    "Month",
    "MonthBegin",
    "MonthEnd",
    "MonthOffset",
    "MultiIndex",
    "NA",
    "NaN",
    "NaT",
    "None",
    "November",
    "NumFOCUS",
    "NumPy",
    "Numba",
    "ORC",
    "October",
    "PEP8",
    "Parquet",
    "Period",
    "PeriodIndex",
    "PyArrow",
    "PyPy",
    "PyTables",
    "Python",
    "Q",
    "QE",
    "QuarterBegin",
    "QuarterEnd",
    "QuarterOffset",
    "R",
    "SAS",
    "SPSS",
    "SQL",
    "STATA",
    "SciPy",
    "SemiMonthBegin",
    "SemiMonthEnd",
    "SemiMonthOffset",
    "September",
    "Series",
    "SeriesGroupBy",
    "SparseDataFrame",
    "Styler",
    "TZ",
    "Timedelta",
    "TimedeltaIndex",
    "Timestamp",
    "UInt64",
    "URLs",
    "UTC",
    "Unicode",
    "VSCode",
    "ValueError",
    "WeekOfMonth",
    "XPORT",
    "XX",
    "Y",
    "YE",
    "YearBegin",
    "YearEnd",
    "YearOffset",
    "msgpack",
    "os",
    "pandas",
    "pd",
    "sklearn",
    "str",
    "strftime",
    "tonumpy"
}

CAP_EXCEPTIONS_DICT = {word.lower(): word for word in CAPITALIZATION_EXCLUSIONS}

err_msg = "Heading capitalization formatted incorrectly. Please correctly capitalize"

symbols = ("*", "=", "-", "^", "~", "#", '"')


def correct_title_capitalization(title: str) -> str:
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

    # Strip all non-word characters from the beginning of the title to the
    # first word character.
    correct_title: str = re.sub(r"^\W*", "", title).capitalize()

    # Remove a URL from the title. We do this because words in a URL must
    # stay lowercase, even if they are a capitalization exception.
    removed_https_title = re.sub(r"<https?:\/\/.*[\r\n]*>", "", correct_title)

    # Split a title into a list using non-word character delimiters.
    word_list = re.split(r"\W", removed_https_title)

    for word in word_list:
        if word.lower() in CAP_EXCEPTIONS_DICT:
            correct_title = re.sub(
                rf"\b{word}\b", CAP_EXCEPTIONS_DICT[word.lower()], correct_title
            )

    return correct_title


def find_titles(rst_file: str) -> Iterable[tuple[str, int]]:
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

    with open(rst_file, encoding="utf-8") as fd:
        previous_line = ""
        for i, line in enumerate(fd):
            line_no_last_elem = line[:-1]
            line_chars = set(line_no_last_elem)
            if (
                    len(line_chars) == 1
                    and line_chars.pop() in symbols
                    and len(line_no_last_elem) == len(previous_line)
            ):
                yield re.sub(r"[`\*_]", "", previous_line), i
            previous_line = line_no_last_elem


def main(source_paths: list[str]) -> int:
    """
    The main method to print all headings with incorrect capitalization.

    Parameters
    ----------
    source_paths : str
        List of directories to validate, provided through command line arguments.

    Returns
    -------
    int
        Number of incorrect headings found overall.
    """

    number_of_errors: int = 0

    for filename in source_paths:
        for title, line_number in find_titles(filename):
            if title != correct_title_capitalization(title):
                print(
                    f"""{filename}:{line_number}:{err_msg} "{title}" to "{
                    correct_title_capitalization(title)}" """
                )
                number_of_errors += 1

    return number_of_errors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate heading capitalization")

    parser.add_argument(
        "paths", nargs="*", help="Source paths of file/directory to check."
    )

    args = parser.parse_args()

    sys.exit(main(args.paths))
