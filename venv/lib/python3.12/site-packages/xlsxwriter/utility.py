###############################################################################
#
# Worksheet - A class for writing Excel Worksheets.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#
import datetime
import re
from typing import Dict, Optional, Tuple, Union
from warnings import warn

from xlsxwriter.color import Color

COL_NAMES: Dict[int, str] = {}

CHAR_WIDTHS = {
    " ": 3,
    "!": 5,
    '"': 6,
    "#": 7,
    "$": 7,
    "%": 11,
    "&": 10,
    "'": 3,
    "(": 5,
    ")": 5,
    "*": 7,
    "+": 7,
    ",": 4,
    "-": 5,
    ".": 4,
    "/": 6,
    "0": 7,
    "1": 7,
    "2": 7,
    "3": 7,
    "4": 7,
    "5": 7,
    "6": 7,
    "7": 7,
    "8": 7,
    "9": 7,
    ":": 4,
    ";": 4,
    "<": 7,
    "=": 7,
    ">": 7,
    "?": 7,
    "@": 13,
    "A": 9,
    "B": 8,
    "C": 8,
    "D": 9,
    "E": 7,
    "F": 7,
    "G": 9,
    "H": 9,
    "I": 4,
    "J": 5,
    "K": 8,
    "L": 6,
    "M": 12,
    "N": 10,
    "O": 10,
    "P": 8,
    "Q": 10,
    "R": 8,
    "S": 7,
    "T": 7,
    "U": 9,
    "V": 9,
    "W": 13,
    "X": 8,
    "Y": 7,
    "Z": 7,
    "[": 5,
    "\\": 6,
    "]": 5,
    "^": 7,
    "_": 7,
    "`": 4,
    "a": 7,
    "b": 8,
    "c": 6,
    "d": 8,
    "e": 8,
    "f": 5,
    "g": 7,
    "h": 8,
    "i": 4,
    "j": 4,
    "k": 7,
    "l": 4,
    "m": 12,
    "n": 8,
    "o": 8,
    "p": 8,
    "q": 8,
    "r": 5,
    "s": 6,
    "t": 5,
    "u": 8,
    "v": 7,
    "w": 11,
    "x": 7,
    "y": 7,
    "z": 6,
    "{": 5,
    "|": 7,
    "}": 5,
    "~": 7,
}

# The following is a list of Emojis used to decide if worksheet names require
# quoting since there is (currently) no native support for matching them in
# Python regular expressions. It is probably unnecessary to exclude them since
# the default quoting is safe in Excel even when unnecessary (the reverse isn't
# true). The Emoji list was generated from:
#
# https://util.unicode.org/UnicodeJsps/list-unicodeset.jsp?a=%5B%3AEmoji%3DYes%3A%5D&abb=on&esc=on&g=&i=
#
# pylint: disable-next=line-too-long
EMOJIS = "\u00a9\u00ae\u203c\u2049\u2122\u2139\u2194-\u2199\u21a9\u21aa\u231a\u231b\u2328\u23cf\u23e9-\u23f3\u23f8-\u23fa\u24c2\u25aa\u25ab\u25b6\u25c0\u25fb-\u25fe\u2600-\u2604\u260e\u2611\u2614\u2615\u2618\u261d\u2620\u2622\u2623\u2626\u262a\u262e\u262f\u2638-\u263a\u2640\u2642\u2648-\u2653\u265f\u2660\u2663\u2665\u2666\u2668\u267b\u267e\u267f\u2692-\u2697\u2699\u269b\u269c\u26a0\u26a1\u26a7\u26aa\u26ab\u26b0\u26b1\u26bd\u26be\u26c4\u26c5\u26c8\u26ce\u26cf\u26d1\u26d3\u26d4\u26e9\u26ea\u26f0-\u26f5\u26f7-\u26fa\u26fd\u2702\u2705\u2708-\u270d\u270f\u2712\u2714\u2716\u271d\u2721\u2728\u2733\u2734\u2744\u2747\u274c\u274e\u2753-\u2755\u2757\u2763\u2764\u2795-\u2797\u27a1\u27b0\u27bf\u2934\u2935\u2b05-\u2b07\u2b1b\u2b1c\u2b50\u2b55\u3030\u303d\u3297\u3299\U0001f004\U0001f0cf\U0001f170\U0001f171\U0001f17e\U0001f17f\U0001f18e\U0001f191-\U0001f19a\U0001f1e6-\U0001f1ff\U0001f201\U0001f202\U0001f21a\U0001f22f\U0001f232-\U0001f23a\U0001f250\U0001f251\U0001f300-\U0001f321\U0001f324-\U0001f393\U0001f396\U0001f397\U0001f399-\U0001f39b\U0001f39e-\U0001f3f0\U0001f3f3-\U0001f3f5\U0001f3f7-\U0001f4fd\U0001f4ff-\U0001f53d\U0001f549-\U0001f54e\U0001f550-\U0001f567\U0001f56f\U0001f570\U0001f573-\U0001f57a\U0001f587\U0001f58a-\U0001f58d\U0001f590\U0001f595\U0001f596\U0001f5a4\U0001f5a5\U0001f5a8\U0001f5b1\U0001f5b2\U0001f5bc\U0001f5c2-\U0001f5c4\U0001f5d1-\U0001f5d3\U0001f5dc-\U0001f5de\U0001f5e1\U0001f5e3\U0001f5e8\U0001f5ef\U0001f5f3\U0001f5fa-\U0001f64f\U0001f680-\U0001f6c5\U0001f6cb-\U0001f6d2\U0001f6d5-\U0001f6d7\U0001f6dc-\U0001f6e5\U0001f6e9\U0001f6eb\U0001f6ec\U0001f6f0\U0001f6f3-\U0001f6fc\U0001f7e0-\U0001f7eb\U0001f7f0\U0001f90c-\U0001f93a\U0001f93c-\U0001f945\U0001f947-\U0001f9ff\U0001fa70-\U0001fa7c\U0001fa80-\U0001fa88\U0001fa90-\U0001fabd\U0001fabf-\U0001fac5\U0001face-\U0001fadb\U0001fae0-\U0001fae8\U0001faf0-\U0001faf8"  # noqa

# Compile performance critical regular expressions.
RE_LEADING_WHITESPACE = re.compile(r"^\s")
RE_TRAILING_WHITESPACE = re.compile(r"\s$")
RE_RANGE_PARTS = re.compile(r"(\$?)([A-Z]{1,3})(\$?)(\d+)")
RE_QUOTE_RULE1 = re.compile(rf"[^\w\.{EMOJIS}]")
RE_QUOTE_RULE2 = re.compile(rf"^[\d\.{EMOJIS}]")
RE_QUOTE_RULE3 = re.compile(r"^([A-Z]{1,3}\d+)$")
RE_QUOTE_RULE4_ROW = re.compile(r"^R(\d+)")
RE_QUOTE_RULE4_COLUMN = re.compile(r"^R?C(\d+)")


def xl_rowcol_to_cell(
    row: int,
    col: int,
    row_abs: bool = False,
    col_abs: bool = False,
) -> str:
    """
    Convert a zero indexed row and column cell reference to a A1 style string.

    Args:
       row:     The cell row.    Int.
       col:     The cell column. Int.
       row_abs: Optional flag to make the row absolute.    Bool.
       col_abs: Optional flag to make the column absolute. Bool.

    Returns:
        A1 style string.

    """
    if row < 0:
        warn(f"Row number '{row}' must be >= 0")
        return ""

    if col < 0:
        warn(f"Col number '{col}' must be >= 0")
        return ""

    row += 1  # Change to 1-index.
    row_abs_str = "$" if row_abs else ""

    col_str = xl_col_to_name(col, col_abs)

    return col_str + row_abs_str + str(row)


def xl_rowcol_to_cell_fast(row: int, col: int) -> str:
    """
    Optimized version of the xl_rowcol_to_cell function. Only used internally.

    Args:
       row: The cell row.    Int.
       col: The cell column. Int.

    Returns:
        A1 style string.

    """
    if col in COL_NAMES:
        col_str = COL_NAMES[col]
    else:
        col_str = xl_col_to_name(col)
        COL_NAMES[col] = col_str

    return col_str + str(row + 1)


def xl_col_to_name(col: int, col_abs: bool = False) -> str:
    """
    Convert a zero indexed column cell reference to a string.

    Args:
       col:     The cell column. Int.
       col_abs: Optional flag to make the column absolute. Bool.

    Returns:
        Column style string.

    """
    col_num = col
    if col_num < 0:
        warn(f"Col number '{col_num}' must be >= 0")
        return ""

    col_num += 1  # Change to 1-index.
    col_str = ""
    col_abs_str = "$" if col_abs else ""

    while col_num:
        # Set remainder from 1 .. 26
        remainder = col_num % 26

        if remainder == 0:
            remainder = 26

        # Convert the remainder to a character.
        col_letter = chr(ord("A") + remainder - 1)

        # Accumulate the column letters, right to left.
        col_str = col_letter + col_str

        # Get the next order of magnitude.
        col_num = int((col_num - 1) / 26)

    return col_abs_str + col_str


def xl_cell_to_rowcol(cell_str: str) -> Tuple[int, int]:
    """
    Convert a cell reference in A1 notation to a zero indexed row and column.

    Args:
       cell_str:  A1 style string.

    Returns:
        row, col: Zero indexed cell row and column indices.

    """
    if not cell_str:
        return 0, 0

    match = RE_RANGE_PARTS.match(cell_str)
    if match is None:
        warn(f"Invalid cell reference '{cell_str}'")
        return 0, 0

    col_str = match.group(2)
    row_str = match.group(4)

    # Convert base26 column string to number.
    expn = 0
    col = 0
    for char in reversed(col_str):
        col += (ord(char) - ord("A") + 1) * (26**expn)
        expn += 1

    # Convert 1-index to zero-index
    row = int(row_str) - 1
    col -= 1

    return row, col


def xl_cell_to_rowcol_abs(cell_str: str) -> Tuple[int, int, bool, bool]:
    """
    Convert an absolute cell reference in A1 notation to a zero indexed
    row and column, with True/False values for absolute rows or columns.

    Args:
       cell_str: A1 style string.

    Returns:
        row, col, row_abs, col_abs:  Zero indexed cell row and column indices.

    """
    if not cell_str:
        return 0, 0, False, False

    match = RE_RANGE_PARTS.match(cell_str)
    if match is None:
        warn(f"Invalid cell reference '{cell_str}'")
        return 0, 0, False, False

    col_abs = bool(match.group(1))
    col_str = match.group(2)
    row_abs = bool(match.group(3))
    row_str = match.group(4)

    # Convert base26 column string to number.
    expn = 0
    col = 0
    for char in reversed(col_str):
        col += (ord(char) - ord("A") + 1) * (26**expn)
        expn += 1

    # Convert 1-index to zero-index
    row = int(row_str) - 1
    col -= 1

    return row, col, row_abs, col_abs


def xl_range(first_row: int, first_col: int, last_row: int, last_col: int) -> str:
    """
    Convert zero indexed row and col cell references to a A1:B1 range string.

    Args:
       first_row: The first cell row.    Int.
       first_col: The first cell column. Int.
       last_row:  The last cell row.     Int.
       last_col:  The last cell column.  Int.

    Returns:
        A1:B1 style range string.

    """
    range1 = xl_rowcol_to_cell(first_row, first_col)
    range2 = xl_rowcol_to_cell(last_row, last_col)

    if range1 == "" or range2 == "":
        warn("Row and column numbers must be >= 0")
        return ""

    if range1 == range2:
        return range1

    return range1 + ":" + range2


def xl_range_abs(first_row: int, first_col: int, last_row: int, last_col: int) -> str:
    """
    Convert zero indexed row and col cell references to a $A$1:$B$1 absolute
    range string.

    Args:
       first_row: The first cell row.    Int.
       first_col: The first cell column. Int.
       last_row:  The last cell row.     Int.
       last_col:  The last cell column.  Int.

    Returns:
        $A$1:$B$1 style range string.

    """
    range1 = xl_rowcol_to_cell(first_row, first_col, True, True)
    range2 = xl_rowcol_to_cell(last_row, last_col, True, True)

    if range1 == "" or range2 == "":
        warn("Row and column numbers must be >= 0")
        return ""

    if range1 == range2:
        return range1

    return range1 + ":" + range2


def xl_range_formula(
    sheetname: str, first_row: int, first_col: int, last_row: int, last_col: int
) -> str:
    """
    Convert worksheet name and zero indexed row and col cell references to
    a Sheet1!A1:B1 range formula string.

    Args:
       sheetname: The worksheet name.    String.
       first_row: The first cell row.    Int.
       first_col: The first cell column. Int.
       last_row:  The last cell row.     Int.
       last_col:  The last cell column.  Int.

    Returns:
        A1:B1 style range string.

    """
    cell_range = xl_range_abs(first_row, first_col, last_row, last_col)
    sheetname = quote_sheetname(sheetname)

    return sheetname + "!" + cell_range


def quote_sheetname(sheetname: str) -> str:
    """
    Sheetnames used in references should be quoted if they contain any spaces,
    special characters or if they look like a A1 or RC cell reference. The rules
    are shown inline below.

    Args:
       sheetname: The worksheet name. String.

    Returns:
        A quoted worksheet string.

    """
    uppercase_sheetname = sheetname.upper()
    requires_quoting = False
    col_max = 163_84
    row_max = 1048576

    # Don't quote sheetname if it is already quoted by the user.
    if not sheetname.startswith("'"):

        match_rule3 = RE_QUOTE_RULE3.match(uppercase_sheetname)
        match_rule4_row = RE_QUOTE_RULE4_ROW.match(uppercase_sheetname)
        match_rule4_column = RE_QUOTE_RULE4_COLUMN.match(uppercase_sheetname)

        # --------------------------------------------------------------------
        # Rule 1. Sheet names that contain anything other than \w and "."
        # characters must be quoted.
        # --------------------------------------------------------------------
        if RE_QUOTE_RULE1.search(sheetname):
            requires_quoting = True

        # --------------------------------------------------------------------
        # Rule 2. Sheet names that start with a digit or "." must be quoted.
        # --------------------------------------------------------------------
        elif RE_QUOTE_RULE2.search(sheetname):
            requires_quoting = True

        # --------------------------------------------------------------------
        # Rule 3. Sheet names must not be a valid A1 style cell reference.
        # Valid means that the row and column range values must also be within
        # Excel row and column limits.
        # --------------------------------------------------------------------
        elif match_rule3:
            cell = match_rule3.group(1)
            (row, col) = xl_cell_to_rowcol(cell)

            if 0 <= row < row_max and 0 <= col < col_max:
                requires_quoting = True

        # --------------------------------------------------------------------
        # Rule 4. Sheet names must not *start* with a valid RC style cell
        # reference. Other characters after the valid RC reference are ignored
        # by Excel. Valid means that the row and column range values must also
        # be within Excel row and column limits.
        #
        # Note: references without trailing characters like R12345 or C12345
        # are caught by Rule 3. Negative references like R-12345 are caught by
        # Rule 1 due to the dash.
        # --------------------------------------------------------------------

        # Rule 4a. Check for sheet names that start with R1 style references.
        elif match_rule4_row:
            row = int(match_rule4_row.group(1))

            if 0 < row <= row_max:
                requires_quoting = True

        # Rule 4b. Check for sheet names that start with C1 or RC1 style
        elif match_rule4_column:
            col = int(match_rule4_column.group(1))

            if 0 < col <= col_max:
                requires_quoting = True

        # Rule 4c. Check for some single R/C references.
        elif uppercase_sheetname in ("R", "C", "RC"):
            requires_quoting = True

    if requires_quoting:
        # Double quote any single quotes.
        sheetname = sheetname.replace("'", "''")

        # Single quote the sheet name.
        sheetname = f"'{sheetname}'"

    return sheetname


def cell_autofit_width(string: str) -> int:
    """
    Calculate the width required to auto-fit a string in a cell.

    Args:
       string: The string to calculate the cell width for. String.

    Returns:
        The string autofit width in pixels. Returns 0 if the string is empty.

    """
    if not string or len(string) == 0:
        return 0

    # Excel adds an additional 7 pixels of padding to the cell boundary.
    return xl_pixel_width(string) + 7


def xl_pixel_width(string: str) -> int:
    """
    Get the pixel width of a string based on individual character widths taken
    from Excel. UTF8 characters, and other unhandled characters, are given a
    default width of 8.

    Args:
       string: The string to calculate the width for. String.

    Returns:
        The string width in pixels. Note, Excel adds an additional 7 pixels of
        padding in the cell.

    """
    length = 0
    for char in string:
        length += CHAR_WIDTHS.get(char, 8)

    return length


def _get_sparkline_style(style_id: int) -> Dict[str, Dict[str, str]]:
    """
    Get the numbered sparkline styles.

    """
    styles = [
        {  # 0
            "low": Color.theme(4, 0),
            "high": Color.theme(4, 0),
            "last": Color.theme(4, 3),
            "first": Color.theme(4, 3),
            "series": Color.theme(4, 5),
            "markers": Color.theme(4, 5),
            "negative": Color.theme(5, 0),
        },
        {  # 1
            "low": Color.theme(4, 0),
            "high": Color.theme(4, 0),
            "last": Color.theme(4, 3),
            "first": Color.theme(4, 3),
            "series": Color.theme(4, 5),
            "markers": Color.theme(4, 5),
            "negative": Color.theme(5, 0),
        },
        {  # 2
            "low": Color.theme(5, 0),
            "high": Color.theme(5, 0),
            "last": Color.theme(5, 3),
            "first": Color.theme(5, 3),
            "series": Color.theme(5, 5),
            "markers": Color.theme(5, 5),
            "negative": Color.theme(6, 0),
        },
        {  # 3
            "low": Color.theme(6, 0),
            "high": Color.theme(6, 0),
            "last": Color.theme(6, 3),
            "first": Color.theme(6, 3),
            "series": Color.theme(6, 5),
            "markers": Color.theme(6, 5),
            "negative": Color.theme(7, 0),
        },
        {  # 4
            "low": Color.theme(7, 0),
            "high": Color.theme(7, 0),
            "last": Color.theme(7, 3),
            "first": Color.theme(7, 3),
            "series": Color.theme(7, 5),
            "markers": Color.theme(7, 5),
            "negative": Color.theme(8, 0),
        },
        {  # 5
            "low": Color.theme(8, 0),
            "high": Color.theme(8, 0),
            "last": Color.theme(8, 3),
            "first": Color.theme(8, 3),
            "series": Color.theme(8, 5),
            "markers": Color.theme(8, 5),
            "negative": Color.theme(9, 0),
        },
        {  # 6
            "low": Color.theme(9, 0),
            "high": Color.theme(9, 0),
            "last": Color.theme(9, 3),
            "first": Color.theme(9, 3),
            "series": Color.theme(9, 5),
            "markers": Color.theme(9, 5),
            "negative": Color.theme(4, 0),
        },
        {  # 7
            "low": Color.theme(5, 4),
            "high": Color.theme(5, 4),
            "last": Color.theme(5, 4),
            "first": Color.theme(5, 4),
            "series": Color.theme(4, 4),
            "markers": Color.theme(5, 4),
            "negative": Color.theme(5, 0),
        },
        {  # 8
            "low": Color.theme(6, 4),
            "high": Color.theme(6, 4),
            "last": Color.theme(6, 4),
            "first": Color.theme(6, 4),
            "series": Color.theme(5, 4),
            "markers": Color.theme(6, 4),
            "negative": Color.theme(6, 0),
        },
        {  # 9
            "low": Color.theme(7, 4),
            "high": Color.theme(7, 4),
            "last": Color.theme(7, 4),
            "first": Color.theme(7, 4),
            "series": Color.theme(6, 4),
            "markers": Color.theme(7, 4),
            "negative": Color.theme(7, 0),
        },
        {  # 10
            "low": Color.theme(8, 4),
            "high": Color.theme(8, 4),
            "last": Color.theme(8, 4),
            "first": Color.theme(8, 4),
            "series": Color.theme(7, 4),
            "markers": Color.theme(8, 4),
            "negative": Color.theme(8, 0),
        },
        {  # 11
            "low": Color.theme(9, 4),
            "high": Color.theme(9, 4),
            "last": Color.theme(9, 4),
            "first": Color.theme(9, 4),
            "series": Color.theme(8, 4),
            "markers": Color.theme(9, 4),
            "negative": Color.theme(9, 0),
        },
        {  # 12
            "low": Color.theme(4, 4),
            "high": Color.theme(4, 4),
            "last": Color.theme(4, 4),
            "first": Color.theme(4, 4),
            "series": Color.theme(9, 4),
            "markers": Color.theme(4, 4),
            "negative": Color.theme(4, 0),
        },
        {  # 13
            "low": Color.theme(4, 4),
            "high": Color.theme(4, 4),
            "last": Color.theme(4, 4),
            "first": Color.theme(4, 4),
            "series": Color.theme(4, 0),
            "markers": Color.theme(4, 4),
            "negative": Color.theme(5, 0),
        },
        {  # 14
            "low": Color.theme(5, 4),
            "high": Color.theme(5, 4),
            "last": Color.theme(5, 4),
            "first": Color.theme(5, 4),
            "series": Color.theme(5, 0),
            "markers": Color.theme(5, 4),
            "negative": Color.theme(6, 0),
        },
        {  # 15
            "low": Color.theme(6, 4),
            "high": Color.theme(6, 4),
            "last": Color.theme(6, 4),
            "first": Color.theme(6, 4),
            "series": Color.theme(6, 0),
            "markers": Color.theme(6, 4),
            "negative": Color.theme(7, 0),
        },
        {  # 16
            "low": Color.theme(7, 4),
            "high": Color.theme(7, 4),
            "last": Color.theme(7, 4),
            "first": Color.theme(7, 4),
            "series": Color.theme(7, 0),
            "markers": Color.theme(7, 4),
            "negative": Color.theme(8, 0),
        },
        {  # 17
            "low": Color.theme(8, 4),
            "high": Color.theme(8, 4),
            "last": Color.theme(8, 4),
            "first": Color.theme(8, 4),
            "series": Color.theme(8, 0),
            "markers": Color.theme(8, 4),
            "negative": Color.theme(9, 0),
        },
        {  # 18
            "low": Color.theme(9, 4),
            "high": Color.theme(9, 4),
            "last": Color.theme(9, 4),
            "first": Color.theme(9, 4),
            "series": Color.theme(9, 0),
            "markers": Color.theme(9, 4),
            "negative": Color.theme(4, 0),
        },
        {  # 19
            "low": Color.theme(4, 5),
            "high": Color.theme(4, 5),
            "last": Color.theme(4, 4),
            "first": Color.theme(4, 4),
            "series": Color.theme(4, 3),
            "markers": Color.theme(4, 1),
            "negative": Color.theme(0, 5),
        },
        {  # 20
            "low": Color.theme(5, 5),
            "high": Color.theme(5, 5),
            "last": Color.theme(5, 4),
            "first": Color.theme(5, 4),
            "series": Color.theme(5, 3),
            "markers": Color.theme(5, 1),
            "negative": Color.theme(0, 5),
        },
        {  # 21
            "low": Color.theme(6, 5),
            "high": Color.theme(6, 5),
            "last": Color.theme(6, 4),
            "first": Color.theme(6, 4),
            "series": Color.theme(6, 3),
            "markers": Color.theme(6, 1),
            "negative": Color.theme(0, 5),
        },
        {  # 22
            "low": Color.theme(7, 5),
            "high": Color.theme(7, 5),
            "last": Color.theme(7, 4),
            "first": Color.theme(7, 4),
            "series": Color.theme(7, 3),
            "markers": Color.theme(7, 1),
            "negative": Color.theme(0, 5),
        },
        {  # 23
            "low": Color.theme(8, 5),
            "high": Color.theme(8, 5),
            "last": Color.theme(8, 4),
            "first": Color.theme(8, 4),
            "series": Color.theme(8, 3),
            "markers": Color.theme(8, 1),
            "negative": Color.theme(0, 5),
        },
        {  # 24
            "low": Color.theme(9, 5),
            "high": Color.theme(9, 5),
            "last": Color.theme(9, 4),
            "first": Color.theme(9, 4),
            "series": Color.theme(9, 3),
            "markers": Color.theme(9, 1),
            "negative": Color.theme(0, 5),
        },
        {  # 25
            "low": Color.theme(1, 3),
            "high": Color.theme(1, 3),
            "last": Color.theme(1, 3),
            "first": Color.theme(1, 3),
            "series": Color.theme(1, 1),
            "markers": Color.theme(1, 3),
            "negative": Color.theme(1, 3),
        },
        {  # 26
            "low": Color.theme(0, 3),
            "high": Color.theme(0, 3),
            "last": Color.theme(0, 3),
            "first": Color.theme(0, 3),
            "series": Color.theme(1, 2),
            "markers": Color.theme(0, 3),
            "negative": Color.theme(0, 3),
        },
        {  # 27
            "low": Color("#D00000"),
            "high": Color("#D00000"),
            "last": Color("#D00000"),
            "first": Color("#D00000"),
            "series": Color("#323232"),
            "markers": Color("#D00000"),
            "negative": Color("#D00000"),
        },
        {  # 28
            "low": Color("#0070C0"),
            "high": Color("#0070C0"),
            "last": Color("#0070C0"),
            "first": Color("#0070C0"),
            "series": Color("#000000"),
            "markers": Color("#0070C0"),
            "negative": Color("#0070C0"),
        },
        {  # 29
            "low": Color("#D00000"),
            "high": Color("#D00000"),
            "last": Color("#D00000"),
            "first": Color("#D00000"),
            "series": Color("#376092"),
            "markers": Color("#D00000"),
            "negative": Color("#D00000"),
        },
        {  # 30
            "low": Color("#000000"),
            "high": Color("#000000"),
            "last": Color("#000000"),
            "first": Color("#000000"),
            "series": Color("#0070C0"),
            "markers": Color("#000000"),
            "negative": Color("#000000"),
        },
        {  # 31
            "low": Color("#FF5055"),
            "high": Color("#56BE79"),
            "last": Color("#359CEB"),
            "first": Color("#5687C2"),
            "series": Color("#5F5F5F"),
            "markers": Color("#D70077"),
            "negative": Color("#FFB620"),
        },
        {  # 32
            "low": Color("#FF5055"),
            "high": Color("#56BE79"),
            "last": Color("#359CEB"),
            "first": Color("#777777"),
            "series": Color("#5687C2"),
            "markers": Color("#D70077"),
            "negative": Color("#FFB620"),
        },
        {  # 33
            "low": Color("#FF5367"),
            "high": Color("#60D276"),
            "last": Color("#FFEB9C"),
            "first": Color("#FFDC47"),
            "series": Color("#C6EFCE"),
            "markers": Color("#8CADD6"),
            "negative": Color("#FFC7CE"),
        },
        {  # 34
            "low": Color("#FF0000"),
            "high": Color("#00B050"),
            "last": Color("#FFC000"),
            "first": Color("#FFC000"),
            "series": Color("#00B050"),
            "markers": Color("#0070C0"),
            "negative": Color("#FF0000"),
        },
        {  # 35
            "low": Color.theme(7, 0),
            "high": Color.theme(6, 0),
            "last": Color.theme(5, 0),
            "first": Color.theme(4, 0),
            "series": Color.theme(3, 0),
            "markers": Color.theme(8, 0),
            "negative": Color.theme(9, 0),
        },
        {  # 36
            "low": Color.theme(7, 0),
            "high": Color.theme(6, 0),
            "last": Color.theme(5, 0),
            "first": Color.theme(4, 0),
            "series": Color.theme(1, 0),
            "markers": Color.theme(8, 0),
            "negative": Color.theme(9, 0),
        },
    ]

    return styles[style_id]


def _supported_datetime(
    dt: Union[datetime.datetime, datetime.time, datetime.date],
) -> bool:
    # Determine is an argument is a supported datetime object.
    return isinstance(
        dt, (datetime.datetime, datetime.date, datetime.time, datetime.timedelta)
    )


def _remove_datetime_timezone(
    dt_obj: datetime.datetime, remove_timezone: bool
) -> datetime.datetime:
    # Excel doesn't support timezones in datetimes/times so we remove the
    # tzinfo from the object if the user has specified that option in the
    # constructor.
    if remove_timezone:
        dt_obj = dt_obj.replace(tzinfo=None)
    else:
        if dt_obj.tzinfo:
            raise TypeError(
                "Excel doesn't support timezones in datetimes. "
                "Set the tzinfo in the datetime/time object to None or "
                "use the 'remove_timezone' Workbook() option"
            )

    return dt_obj


def _datetime_to_excel_datetime(
    dt_obj: Union[datetime.time, datetime.datetime, datetime.timedelta, datetime.date],
    date_1904: bool,
    remove_timezone: bool,
) -> float:
    # Convert a datetime object to an Excel serial date and time. The integer
    # part of the number stores the number of days since the epoch and the
    # fractional part stores the percentage of the day.
    date_type = dt_obj
    is_timedelta = False

    if date_1904:
        # Excel for Mac date epoch.
        epoch = datetime.datetime(1904, 1, 1)
    else:
        # Default Excel epoch.
        epoch = datetime.datetime(1899, 12, 31)

    # We handle datetime .datetime, .date and .time objects but convert
    # them to datetime.datetime objects and process them in the same way.
    if isinstance(dt_obj, datetime.datetime):
        dt_obj = _remove_datetime_timezone(dt_obj, remove_timezone)
        delta = dt_obj - epoch
    elif isinstance(dt_obj, datetime.date):
        dt_obj = datetime.datetime.fromordinal(dt_obj.toordinal())
        delta = dt_obj - epoch
    elif isinstance(dt_obj, datetime.time):
        dt_obj = datetime.datetime.combine(epoch, dt_obj)
        dt_obj = _remove_datetime_timezone(dt_obj, remove_timezone)
        delta = dt_obj - epoch
    elif isinstance(dt_obj, datetime.timedelta):
        is_timedelta = True
        delta = dt_obj
    else:
        raise TypeError("Unknown or unsupported datetime type")

    # Convert a Python datetime.datetime value to an Excel date number.
    excel_time = delta.days + (
        float(delta.seconds) + float(delta.microseconds) / 1e6
    ) / (60 * 60 * 24)

    # The following is a workaround for the fact that in Excel a time only
    # value is represented as 1899-12-31+time whereas in datetime.datetime()
    # it is 1900-1-1+time so we need to subtract the 1 day difference.
    if (
        isinstance(date_type, datetime.datetime)
        and not isinstance(dt_obj, datetime.timedelta)
        and dt_obj.isocalendar()
        == (
            1900,
            1,
            1,
        )
    ):
        excel_time -= 1

    # Account for Excel erroneously treating 1900 as a leap year.
    if not date_1904 and not is_timedelta and excel_time > 59:
        excel_time += 1

    return excel_time


def _preserve_whitespace(string: str) -> Optional[re.Match]:
    # Check if a string has leading or trailing whitespace that requires a
    # "preserve" attribute.
    return RE_LEADING_WHITESPACE.search(string) or RE_TRAILING_WHITESPACE.search(string)
