###############################################################################
#
# Worksheet - A class for writing Excel Worksheets.
#
# SPDX-License-Identifier: BSD-2-Clause
# Copyright 2013-2024, John McNamara, jmcnamara@cpan.org
#
import datetime
import hashlib
import os
import re
from struct import unpack
from warnings import warn
from .exceptions import UndefinedImageSize
from .exceptions import UnsupportedImageFormat

COL_NAMES = {}

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

# Compile performance critical regular expressions.
re_leading = re.compile(r"^\s")
re_trailing = re.compile(r"\s$")
re_range_parts = re.compile(r"(\$?)([A-Z]{1,3})(\$?)(\d+)")


def xl_rowcol_to_cell(row, col, row_abs=False, col_abs=False):
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
        warn("Row number %d must be >= 0" % row)
        return None

    if col < 0:
        warn("Col number %d must be >= 0" % col)
        return None

    row += 1  # Change to 1-index.
    row_abs = "$" if row_abs else ""

    col_str = xl_col_to_name(col, col_abs)

    return col_str + row_abs + str(row)


def xl_rowcol_to_cell_fast(row, col):
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


def xl_col_to_name(col, col_abs=False):
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
        warn("Col number %d must be >= 0" % col_num)
        return None

    col_num += 1  # Change to 1-index.
    col_str = ""
    col_abs = "$" if col_abs else ""

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

    return col_abs + col_str


def xl_cell_to_rowcol(cell_str):
    """
    Convert a cell reference in A1 notation to a zero indexed row and column.

    Args:
       cell_str:  A1 style string.

    Returns:
        row, col: Zero indexed cell row and column indices.

    """
    if not cell_str:
        return 0, 0

    match = re_range_parts.match(cell_str)
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


def xl_cell_to_rowcol_abs(cell_str):
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

    match = re_range_parts.match(cell_str)

    col_abs = match.group(1)
    col_str = match.group(2)
    row_abs = match.group(3)
    row_str = match.group(4)

    if col_abs:
        col_abs = True
    else:
        col_abs = False

    if row_abs:
        row_abs = True
    else:
        row_abs = False

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


def xl_range(first_row, first_col, last_row, last_col):
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

    if range1 is None or range2 is None:
        warn("Row and column numbers must be >= 0")
        return None

    if range1 == range2:
        return range1
    else:
        return range1 + ":" + range2


def xl_range_abs(first_row, first_col, last_row, last_col):
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

    if range1 is None or range2 is None:
        warn("Row and column numbers must be >= 0")
        return None

    if range1 == range2:
        return range1
    else:
        return range1 + ":" + range2


def xl_range_formula(sheetname, first_row, first_col, last_row, last_col):
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


def quote_sheetname(sheetname):
    """
    Convert a worksheet name to a quoted  name if it contains spaces or
    special characters.

    Args:
       sheetname: The worksheet name. String.

    Returns:
        A quoted worksheet string.

    """

    if not sheetname.isalnum() and not sheetname.startswith("'"):
        # Double quote any single quotes.
        sheetname = sheetname.replace("'", "''")

        # Single quote the sheet name.
        sheetname = "'%s'" % sheetname

    return sheetname


def xl_pixel_width(string):
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


def xl_color(color):
    # Used in conjunction with the XlsxWriter *color() methods to convert
    # a color name into an RGB formatted string. These colors are for
    # backward compatibility with older versions of Excel.
    named_colors = {
        "black": "#000000",
        "blue": "#0000FF",
        "brown": "#800000",
        "cyan": "#00FFFF",
        "gray": "#808080",
        "green": "#008000",
        "lime": "#00FF00",
        "magenta": "#FF00FF",
        "navy": "#000080",
        "orange": "#FF6600",
        "pink": "#FF00FF",
        "purple": "#800080",
        "red": "#FF0000",
        "silver": "#C0C0C0",
        "white": "#FFFFFF",
        "yellow": "#FFFF00",
    }

    if color in named_colors:
        color = named_colors[color]

    if not re.match("#[0-9a-fA-F]{6}", color):
        warn("Color '%s' isn't a valid Excel color" % color)

    # Convert the RGB color to the Excel ARGB format.
    return "FF" + color.lstrip("#").upper()


def get_rgb_color(color):
    # Convert the user specified color to an RGB color.
    rgb_color = xl_color(color)

    # Remove leading FF from RGB color for charts.
    rgb_color = re.sub(r"^FF", "", rgb_color)

    return rgb_color


def get_sparkline_style(style_id):
    styles = [
        {
            "series": {"theme": "4", "tint": "-0.499984740745262"},
            "negative": {"theme": "5"},
            "markers": {"theme": "4", "tint": "-0.499984740745262"},
            "first": {"theme": "4", "tint": "0.39997558519241921"},
            "last": {"theme": "4", "tint": "0.39997558519241921"},
            "high": {"theme": "4"},
            "low": {"theme": "4"},
        },  # 0
        {
            "series": {"theme": "4", "tint": "-0.499984740745262"},
            "negative": {"theme": "5"},
            "markers": {"theme": "4", "tint": "-0.499984740745262"},
            "first": {"theme": "4", "tint": "0.39997558519241921"},
            "last": {"theme": "4", "tint": "0.39997558519241921"},
            "high": {"theme": "4"},
            "low": {"theme": "4"},
        },  # 1
        {
            "series": {"theme": "5", "tint": "-0.499984740745262"},
            "negative": {"theme": "6"},
            "markers": {"theme": "5", "tint": "-0.499984740745262"},
            "first": {"theme": "5", "tint": "0.39997558519241921"},
            "last": {"theme": "5", "tint": "0.39997558519241921"},
            "high": {"theme": "5"},
            "low": {"theme": "5"},
        },  # 2
        {
            "series": {"theme": "6", "tint": "-0.499984740745262"},
            "negative": {"theme": "7"},
            "markers": {"theme": "6", "tint": "-0.499984740745262"},
            "first": {"theme": "6", "tint": "0.39997558519241921"},
            "last": {"theme": "6", "tint": "0.39997558519241921"},
            "high": {"theme": "6"},
            "low": {"theme": "6"},
        },  # 3
        {
            "series": {"theme": "7", "tint": "-0.499984740745262"},
            "negative": {"theme": "8"},
            "markers": {"theme": "7", "tint": "-0.499984740745262"},
            "first": {"theme": "7", "tint": "0.39997558519241921"},
            "last": {"theme": "7", "tint": "0.39997558519241921"},
            "high": {"theme": "7"},
            "low": {"theme": "7"},
        },  # 4
        {
            "series": {"theme": "8", "tint": "-0.499984740745262"},
            "negative": {"theme": "9"},
            "markers": {"theme": "8", "tint": "-0.499984740745262"},
            "first": {"theme": "8", "tint": "0.39997558519241921"},
            "last": {"theme": "8", "tint": "0.39997558519241921"},
            "high": {"theme": "8"},
            "low": {"theme": "8"},
        },  # 5
        {
            "series": {"theme": "9", "tint": "-0.499984740745262"},
            "negative": {"theme": "4"},
            "markers": {"theme": "9", "tint": "-0.499984740745262"},
            "first": {"theme": "9", "tint": "0.39997558519241921"},
            "last": {"theme": "9", "tint": "0.39997558519241921"},
            "high": {"theme": "9"},
            "low": {"theme": "9"},
        },  # 6
        {
            "series": {"theme": "4", "tint": "-0.249977111117893"},
            "negative": {"theme": "5"},
            "markers": {"theme": "5", "tint": "-0.249977111117893"},
            "first": {"theme": "5", "tint": "-0.249977111117893"},
            "last": {"theme": "5", "tint": "-0.249977111117893"},
            "high": {"theme": "5", "tint": "-0.249977111117893"},
            "low": {"theme": "5", "tint": "-0.249977111117893"},
        },  # 7
        {
            "series": {"theme": "5", "tint": "-0.249977111117893"},
            "negative": {"theme": "6"},
            "markers": {"theme": "6", "tint": "-0.249977111117893"},
            "first": {"theme": "6", "tint": "-0.249977111117893"},
            "last": {"theme": "6", "tint": "-0.249977111117893"},
            "high": {"theme": "6", "tint": "-0.249977111117893"},
            "low": {"theme": "6", "tint": "-0.249977111117893"},
        },  # 8
        {
            "series": {"theme": "6", "tint": "-0.249977111117893"},
            "negative": {"theme": "7"},
            "markers": {"theme": "7", "tint": "-0.249977111117893"},
            "first": {"theme": "7", "tint": "-0.249977111117893"},
            "last": {"theme": "7", "tint": "-0.249977111117893"},
            "high": {"theme": "7", "tint": "-0.249977111117893"},
            "low": {"theme": "7", "tint": "-0.249977111117893"},
        },  # 9
        {
            "series": {"theme": "7", "tint": "-0.249977111117893"},
            "negative": {"theme": "8"},
            "markers": {"theme": "8", "tint": "-0.249977111117893"},
            "first": {"theme": "8", "tint": "-0.249977111117893"},
            "last": {"theme": "8", "tint": "-0.249977111117893"},
            "high": {"theme": "8", "tint": "-0.249977111117893"},
            "low": {"theme": "8", "tint": "-0.249977111117893"},
        },  # 10
        {
            "series": {"theme": "8", "tint": "-0.249977111117893"},
            "negative": {"theme": "9"},
            "markers": {"theme": "9", "tint": "-0.249977111117893"},
            "first": {"theme": "9", "tint": "-0.249977111117893"},
            "last": {"theme": "9", "tint": "-0.249977111117893"},
            "high": {"theme": "9", "tint": "-0.249977111117893"},
            "low": {"theme": "9", "tint": "-0.249977111117893"},
        },  # 11
        {
            "series": {"theme": "9", "tint": "-0.249977111117893"},
            "negative": {"theme": "4"},
            "markers": {"theme": "4", "tint": "-0.249977111117893"},
            "first": {"theme": "4", "tint": "-0.249977111117893"},
            "last": {"theme": "4", "tint": "-0.249977111117893"},
            "high": {"theme": "4", "tint": "-0.249977111117893"},
            "low": {"theme": "4", "tint": "-0.249977111117893"},
        },  # 12
        {
            "series": {"theme": "4"},
            "negative": {"theme": "5"},
            "markers": {"theme": "4", "tint": "-0.249977111117893"},
            "first": {"theme": "4", "tint": "-0.249977111117893"},
            "last": {"theme": "4", "tint": "-0.249977111117893"},
            "high": {"theme": "4", "tint": "-0.249977111117893"},
            "low": {"theme": "4", "tint": "-0.249977111117893"},
        },  # 13
        {
            "series": {"theme": "5"},
            "negative": {"theme": "6"},
            "markers": {"theme": "5", "tint": "-0.249977111117893"},
            "first": {"theme": "5", "tint": "-0.249977111117893"},
            "last": {"theme": "5", "tint": "-0.249977111117893"},
            "high": {"theme": "5", "tint": "-0.249977111117893"},
            "low": {"theme": "5", "tint": "-0.249977111117893"},
        },  # 14
        {
            "series": {"theme": "6"},
            "negative": {"theme": "7"},
            "markers": {"theme": "6", "tint": "-0.249977111117893"},
            "first": {"theme": "6", "tint": "-0.249977111117893"},
            "last": {"theme": "6", "tint": "-0.249977111117893"},
            "high": {"theme": "6", "tint": "-0.249977111117893"},
            "low": {"theme": "6", "tint": "-0.249977111117893"},
        },  # 15
        {
            "series": {"theme": "7"},
            "negative": {"theme": "8"},
            "markers": {"theme": "7", "tint": "-0.249977111117893"},
            "first": {"theme": "7", "tint": "-0.249977111117893"},
            "last": {"theme": "7", "tint": "-0.249977111117893"},
            "high": {"theme": "7", "tint": "-0.249977111117893"},
            "low": {"theme": "7", "tint": "-0.249977111117893"},
        },  # 16
        {
            "series": {"theme": "8"},
            "negative": {"theme": "9"},
            "markers": {"theme": "8", "tint": "-0.249977111117893"},
            "first": {"theme": "8", "tint": "-0.249977111117893"},
            "last": {"theme": "8", "tint": "-0.249977111117893"},
            "high": {"theme": "8", "tint": "-0.249977111117893"},
            "low": {"theme": "8", "tint": "-0.249977111117893"},
        },  # 17
        {
            "series": {"theme": "9"},
            "negative": {"theme": "4"},
            "markers": {"theme": "9", "tint": "-0.249977111117893"},
            "first": {"theme": "9", "tint": "-0.249977111117893"},
            "last": {"theme": "9", "tint": "-0.249977111117893"},
            "high": {"theme": "9", "tint": "-0.249977111117893"},
            "low": {"theme": "9", "tint": "-0.249977111117893"},
        },  # 18
        {
            "series": {"theme": "4", "tint": "0.39997558519241921"},
            "negative": {"theme": "0", "tint": "-0.499984740745262"},
            "markers": {"theme": "4", "tint": "0.79998168889431442"},
            "first": {"theme": "4", "tint": "-0.249977111117893"},
            "last": {"theme": "4", "tint": "-0.249977111117893"},
            "high": {"theme": "4", "tint": "-0.499984740745262"},
            "low": {"theme": "4", "tint": "-0.499984740745262"},
        },  # 19
        {
            "series": {"theme": "5", "tint": "0.39997558519241921"},
            "negative": {"theme": "0", "tint": "-0.499984740745262"},
            "markers": {"theme": "5", "tint": "0.79998168889431442"},
            "first": {"theme": "5", "tint": "-0.249977111117893"},
            "last": {"theme": "5", "tint": "-0.249977111117893"},
            "high": {"theme": "5", "tint": "-0.499984740745262"},
            "low": {"theme": "5", "tint": "-0.499984740745262"},
        },  # 20
        {
            "series": {"theme": "6", "tint": "0.39997558519241921"},
            "negative": {"theme": "0", "tint": "-0.499984740745262"},
            "markers": {"theme": "6", "tint": "0.79998168889431442"},
            "first": {"theme": "6", "tint": "-0.249977111117893"},
            "last": {"theme": "6", "tint": "-0.249977111117893"},
            "high": {"theme": "6", "tint": "-0.499984740745262"},
            "low": {"theme": "6", "tint": "-0.499984740745262"},
        },  # 21
        {
            "series": {"theme": "7", "tint": "0.39997558519241921"},
            "negative": {"theme": "0", "tint": "-0.499984740745262"},
            "markers": {"theme": "7", "tint": "0.79998168889431442"},
            "first": {"theme": "7", "tint": "-0.249977111117893"},
            "last": {"theme": "7", "tint": "-0.249977111117893"},
            "high": {"theme": "7", "tint": "-0.499984740745262"},
            "low": {"theme": "7", "tint": "-0.499984740745262"},
        },  # 22
        {
            "series": {"theme": "8", "tint": "0.39997558519241921"},
            "negative": {"theme": "0", "tint": "-0.499984740745262"},
            "markers": {"theme": "8", "tint": "0.79998168889431442"},
            "first": {"theme": "8", "tint": "-0.249977111117893"},
            "last": {"theme": "8", "tint": "-0.249977111117893"},
            "high": {"theme": "8", "tint": "-0.499984740745262"},
            "low": {"theme": "8", "tint": "-0.499984740745262"},
        },  # 23
        {
            "series": {"theme": "9", "tint": "0.39997558519241921"},
            "negative": {"theme": "0", "tint": "-0.499984740745262"},
            "markers": {"theme": "9", "tint": "0.79998168889431442"},
            "first": {"theme": "9", "tint": "-0.249977111117893"},
            "last": {"theme": "9", "tint": "-0.249977111117893"},
            "high": {"theme": "9", "tint": "-0.499984740745262"},
            "low": {"theme": "9", "tint": "-0.499984740745262"},
        },  # 24
        {
            "series": {"theme": "1", "tint": "0.499984740745262"},
            "negative": {"theme": "1", "tint": "0.249977111117893"},
            "markers": {"theme": "1", "tint": "0.249977111117893"},
            "first": {"theme": "1", "tint": "0.249977111117893"},
            "last": {"theme": "1", "tint": "0.249977111117893"},
            "high": {"theme": "1", "tint": "0.249977111117893"},
            "low": {"theme": "1", "tint": "0.249977111117893"},
        },  # 25
        {
            "series": {"theme": "1", "tint": "0.34998626667073579"},
            "negative": {"theme": "0", "tint": "-0.249977111117893"},
            "markers": {"theme": "0", "tint": "-0.249977111117893"},
            "first": {"theme": "0", "tint": "-0.249977111117893"},
            "last": {"theme": "0", "tint": "-0.249977111117893"},
            "high": {"theme": "0", "tint": "-0.249977111117893"},
            "low": {"theme": "0", "tint": "-0.249977111117893"},
        },  # 26
        {
            "series": {"rgb": "FF323232"},
            "negative": {"rgb": "FFD00000"},
            "markers": {"rgb": "FFD00000"},
            "first": {"rgb": "FFD00000"},
            "last": {"rgb": "FFD00000"},
            "high": {"rgb": "FFD00000"},
            "low": {"rgb": "FFD00000"},
        },  # 27
        {
            "series": {"rgb": "FF000000"},
            "negative": {"rgb": "FF0070C0"},
            "markers": {"rgb": "FF0070C0"},
            "first": {"rgb": "FF0070C0"},
            "last": {"rgb": "FF0070C0"},
            "high": {"rgb": "FF0070C0"},
            "low": {"rgb": "FF0070C0"},
        },  # 28
        {
            "series": {"rgb": "FF376092"},
            "negative": {"rgb": "FFD00000"},
            "markers": {"rgb": "FFD00000"},
            "first": {"rgb": "FFD00000"},
            "last": {"rgb": "FFD00000"},
            "high": {"rgb": "FFD00000"},
            "low": {"rgb": "FFD00000"},
        },  # 29
        {
            "series": {"rgb": "FF0070C0"},
            "negative": {"rgb": "FF000000"},
            "markers": {"rgb": "FF000000"},
            "first": {"rgb": "FF000000"},
            "last": {"rgb": "FF000000"},
            "high": {"rgb": "FF000000"},
            "low": {"rgb": "FF000000"},
        },  # 30
        {
            "series": {"rgb": "FF5F5F5F"},
            "negative": {"rgb": "FFFFB620"},
            "markers": {"rgb": "FFD70077"},
            "first": {"rgb": "FF5687C2"},
            "last": {"rgb": "FF359CEB"},
            "high": {"rgb": "FF56BE79"},
            "low": {"rgb": "FFFF5055"},
        },  # 31
        {
            "series": {"rgb": "FF5687C2"},
            "negative": {"rgb": "FFFFB620"},
            "markers": {"rgb": "FFD70077"},
            "first": {"rgb": "FF777777"},
            "last": {"rgb": "FF359CEB"},
            "high": {"rgb": "FF56BE79"},
            "low": {"rgb": "FFFF5055"},
        },  # 32
        {
            "series": {"rgb": "FFC6EFCE"},
            "negative": {"rgb": "FFFFC7CE"},
            "markers": {"rgb": "FF8CADD6"},
            "first": {"rgb": "FFFFDC47"},
            "last": {"rgb": "FFFFEB9C"},
            "high": {"rgb": "FF60D276"},
            "low": {"rgb": "FFFF5367"},
        },  # 33
        {
            "series": {"rgb": "FF00B050"},
            "negative": {"rgb": "FFFF0000"},
            "markers": {"rgb": "FF0070C0"},
            "first": {"rgb": "FFFFC000"},
            "last": {"rgb": "FFFFC000"},
            "high": {"rgb": "FF00B050"},
            "low": {"rgb": "FFFF0000"},
        },  # 34
        {
            "series": {"theme": "3"},
            "negative": {"theme": "9"},
            "markers": {"theme": "8"},
            "first": {"theme": "4"},
            "last": {"theme": "5"},
            "high": {"theme": "6"},
            "low": {"theme": "7"},
        },  # 35
        {
            "series": {"theme": "1"},
            "negative": {"theme": "9"},
            "markers": {"theme": "8"},
            "first": {"theme": "4"},
            "last": {"theme": "5"},
            "high": {"theme": "6"},
            "low": {"theme": "7"},
        },  # 36
    ]

    return styles[style_id]


def supported_datetime(dt):
    # Determine is an argument is a supported datetime object.
    return isinstance(
        dt, (datetime.datetime, datetime.date, datetime.time, datetime.timedelta)
    )


def remove_datetime_timezone(dt_obj, remove_timezone):
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


def datetime_to_excel_datetime(dt_obj, date_1904, remove_timezone):
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
        dt_obj = remove_datetime_timezone(dt_obj, remove_timezone)
        delta = dt_obj - epoch
    elif isinstance(dt_obj, datetime.date):
        dt_obj = datetime.datetime.fromordinal(dt_obj.toordinal())
        delta = dt_obj - epoch
    elif isinstance(dt_obj, datetime.time):
        dt_obj = datetime.datetime.combine(epoch, dt_obj)
        dt_obj = remove_datetime_timezone(dt_obj, remove_timezone)
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
    if isinstance(date_type, datetime.datetime) and dt_obj.isocalendar() == (
        1900,
        1,
        1,
    ):
        excel_time -= 1

    # Account for Excel erroneously treating 1900 as a leap year.
    if not date_1904 and not is_timedelta and excel_time > 59:
        excel_time += 1

    return excel_time


def preserve_whitespace(string):
    # Check if a string has leading or trailing whitespace that requires a
    # "preserve" attribute.
    if re_leading.search(string) or re_trailing.search(string):
        return True
    else:
        return False


def get_image_properties(filename, image_data):
    # Extract dimension information from the image file.
    height = 0
    width = 0
    x_dpi = 96
    y_dpi = 96

    if not image_data:
        # Open the image file and read in the data.
        fh = open(filename, "rb")
        data = fh.read()
    else:
        # Read the image data from the user supplied byte stream.
        data = image_data.getvalue()

    digest = hashlib.sha256(data).hexdigest()

    # Get the image filename without the path.
    image_name = os.path.basename(filename)

    # Look for some common image file markers.
    marker1 = unpack("3s", data[1:4])[0]
    marker2 = unpack(">H", data[:2])[0]
    marker3 = unpack("2s", data[:2])[0]
    marker4 = unpack("<L", data[:4])[0]
    marker5 = (unpack("4s", data[40:44]))[0]
    marker6 = unpack("4s", data[:4])[0]

    png_marker = b"PNG"
    bmp_marker = b"BM"
    emf_marker = b" EMF"
    gif_marker = b"GIF8"

    if marker1 == png_marker:
        (image_type, width, height, x_dpi, y_dpi) = _process_png(data)

    elif marker2 == 0xFFD8:
        (image_type, width, height, x_dpi, y_dpi) = _process_jpg(data)

    elif marker3 == bmp_marker:
        (image_type, width, height) = _process_bmp(data)

    elif marker4 == 0x9AC6CDD7:
        (image_type, width, height, x_dpi, y_dpi) = _process_wmf(data)

    elif marker4 == 1 and marker5 == emf_marker:
        (image_type, width, height, x_dpi, y_dpi) = _process_emf(data)

    elif marker6 == gif_marker:
        (image_type, width, height, x_dpi, y_dpi) = _process_gif(data)

    else:
        raise UnsupportedImageFormat(
            "%s: Unknown or unsupported image file format." % filename
        )

    # Check that we found the required data.
    if not height or not width:
        raise UndefinedImageSize("%s: no size data found in image file." % filename)

    if not image_data:
        fh.close()

    # Set a default dpi for images with 0 dpi.
    if x_dpi == 0:
        x_dpi = 96
    if y_dpi == 0:
        y_dpi = 96

    return image_type, width, height, image_name, x_dpi, y_dpi, digest


def _process_png(data):
    # Extract width and height information from a PNG file.
    offset = 8
    data_length = len(data)
    end_marker = False
    width = 0
    height = 0
    x_dpi = 96
    y_dpi = 96

    # Search through the image data to read the height and width in the
    # IHDR element. Also read the DPI in the pHYs element.
    while not end_marker and offset < data_length:
        length = unpack(">I", data[offset + 0 : offset + 4])[0]
        marker = unpack("4s", data[offset + 4 : offset + 8])[0]

        # Read the image dimensions.
        if marker == b"IHDR":
            width = unpack(">I", data[offset + 8 : offset + 12])[0]
            height = unpack(">I", data[offset + 12 : offset + 16])[0]

        # Read the image DPI.
        if marker == b"pHYs":
            x_density = unpack(">I", data[offset + 8 : offset + 12])[0]
            y_density = unpack(">I", data[offset + 12 : offset + 16])[0]
            units = unpack("b", data[offset + 16 : offset + 17])[0]

            if units == 1:
                x_dpi = x_density * 0.0254
                y_dpi = y_density * 0.0254

        if marker == b"IEND":
            end_marker = True
            continue

        offset = offset + length + 12

    return "png", width, height, x_dpi, y_dpi


def _process_jpg(data):
    # Extract width and height information from a JPEG file.
    offset = 2
    data_length = len(data)
    end_marker = False
    width = 0
    height = 0
    x_dpi = 96
    y_dpi = 96

    # Search through the image data to read the JPEG markers.
    while not end_marker and offset < data_length:
        marker = unpack(">H", data[offset + 0 : offset + 2])[0]
        length = unpack(">H", data[offset + 2 : offset + 4])[0]

        # Read the height and width in the 0xFFCn elements (except C4, C8
        # and CC which aren't SOF markers).
        if (
            (marker & 0xFFF0) == 0xFFC0
            and marker != 0xFFC4
            and marker != 0xFFC8
            and marker != 0xFFCC
        ):
            height = unpack(">H", data[offset + 5 : offset + 7])[0]
            width = unpack(">H", data[offset + 7 : offset + 9])[0]

        # Read the DPI in the 0xFFE0 element.
        if marker == 0xFFE0:
            units = unpack("b", data[offset + 11 : offset + 12])[0]
            x_density = unpack(">H", data[offset + 12 : offset + 14])[0]
            y_density = unpack(">H", data[offset + 14 : offset + 16])[0]

            if units == 1:
                x_dpi = x_density
                y_dpi = y_density

            if units == 2:
                x_dpi = x_density * 2.54
                y_dpi = y_density * 2.54

            # Workaround for incorrect dpi.
            if x_dpi == 1:
                x_dpi = 96
            if y_dpi == 1:
                y_dpi = 96

        if marker == 0xFFDA:
            end_marker = True
            continue

        offset = offset + length + 2

    return "jpeg", width, height, x_dpi, y_dpi


def _process_gif(data):
    # Extract width and height information from a GIF file.
    x_dpi = 96
    y_dpi = 96

    width = unpack("<h", data[6:8])[0]
    height = unpack("<h", data[8:10])[0]

    return "gif", width, height, x_dpi, y_dpi


def _process_bmp(data):
    # Extract width and height information from a BMP file.
    width = unpack("<L", data[18:22])[0]
    height = unpack("<L", data[22:26])[0]
    return "bmp", width, height


def _process_wmf(data):
    # Extract width and height information from a WMF file.
    x_dpi = 96
    y_dpi = 96

    # Read the bounding box, measured in logical units.
    x1 = unpack("<h", data[6:8])[0]
    y1 = unpack("<h", data[8:10])[0]
    x2 = unpack("<h", data[10:12])[0]
    y2 = unpack("<h", data[12:14])[0]

    # Read the number of logical units per inch. Used to scale the image.
    inch = unpack("<H", data[14:16])[0]

    # Convert to rendered height and width.
    width = float((x2 - x1) * x_dpi) / inch
    height = float((y2 - y1) * y_dpi) / inch

    return "wmf", width, height, x_dpi, y_dpi


def _process_emf(data):
    # Extract width and height information from a EMF file.

    # Read the bounding box, measured in logical units.
    bound_x1 = unpack("<l", data[8:12])[0]
    bound_y1 = unpack("<l", data[12:16])[0]
    bound_x2 = unpack("<l", data[16:20])[0]
    bound_y2 = unpack("<l", data[20:24])[0]

    # Convert the bounds to width and height.
    width = bound_x2 - bound_x1
    height = bound_y2 - bound_y1

    # Read the rectangular frame in units of 0.01mm.
    frame_x1 = unpack("<l", data[24:28])[0]
    frame_y1 = unpack("<l", data[28:32])[0]
    frame_x2 = unpack("<l", data[32:36])[0]
    frame_y2 = unpack("<l", data[36:40])[0]

    # Convert the frame bounds to mm width and height.
    width_mm = 0.01 * (frame_x2 - frame_x1)
    height_mm = 0.01 * (frame_y2 - frame_y1)

    # Get the dpi based on the logical size.
    x_dpi = width * 25.4 / width_mm
    y_dpi = height * 25.4 / height_mm

    # This is to match Excel's calculation. It is probably to account for
    # the fact that the bounding box is inclusive-inclusive. Or a bug.
    width += 1
    height += 1

    return "emf", width, height, x_dpi, y_dpi
