from collections.abc import Callable
from typing import Final, Literal

from .biffh import BaseObject
from .book import Book
from .timemachine import *

DEBUG: Final[int]
excel_default_palette_b5: Final[tuple[tuple[int, int, int], ...]]
excel_default_palette_b2: Final[tuple[tuple[int, int, int], ...]]
excel_default_palette_b8: Final[tuple[tuple[int, int, int], ...]]
default_palette: Final[dict[int, tuple[tuple[int, int, int], ...]]]
built_in_style_names: Final[list[str]]

def initialise_colour_map(book: Book) -> None: ...
def nearest_colour_index(
    colour_map: dict[int, tuple[int, int, int] | None], rgb: tuple[int, int, int] | None, debug: int = 0
) -> int: ...

class EqNeAttrs:
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...

class Font(BaseObject, EqNeAttrs):
    bold: Literal[0, 1]
    character_set: int
    colour_index: int
    escapement: Literal[0, 1, 2]
    family: Literal[0, 1, 2, 3, 4, 5]
    font_index: int
    height: int
    italic: Literal[0, 1]
    name: str
    struck_out: Literal[0, 1]
    underline_type: Literal[0, 1, 33, 34]
    underlined: Literal[0, 1]
    weight: int
    outline: Literal[0, 1]
    shadow: Literal[0, 1]

def handle_efont(book: Book, data: bytes) -> None: ...
def handle_font(book: Book, data: bytes) -> None: ...

class Format(BaseObject, EqNeAttrs):
    format_key: int
    type: int
    format_str: str
    def __init__(self, format_key: int, ty: int, format_str: str) -> None: ...

std_format_strings: dict[int, str]
fmt_code_ranges: list[tuple[int, int, int]]
std_format_code_types: dict[int, int]
date_char_dict: dict[str, Literal[5]]
skip_char_dict: dict[str, Literal[1]]
num_char_dict: dict[str, Literal[5]]
non_date_formats: dict[str, Literal[1]]
fmt_bracketed_sub: Callable[[str, str], str]

def is_date_format_string(book: Book, fmt: str) -> bool: ...
def handle_format(self: Book, data: bytes, rectype: int = 1054) -> None: ...
def handle_palette(book: Book, data: bytes) -> None: ...
def palette_epilogue(book: Book) -> None: ...
def handle_style(book: Book, data: bytes) -> None: ...
def check_colour_indexes_in_obj(book: Book, obj: object, orig_index: int) -> None: ...
def fill_in_standard_formats(book: Book) -> None: ...
def handle_xf(self: Book, data: bytes) -> None: ...
def xf_epilogue(self: Book) -> None: ...
def initialise_book(book: Book) -> None: ...

class XFBorder(BaseObject, EqNeAttrs):
    top_colour_index: int
    bottom_colour_index: int
    left_colour_index: int
    right_colour_index: int
    diag_colour_index: int
    top_line_style: int
    bottom_line_style: int
    left_line_style: int
    right_line_style: int
    diag_line_style: int
    diag_down: Literal[0, 1]
    diag_up: Literal[0, 1]

class XFBackground(BaseObject, EqNeAttrs):
    fill_pattern: int
    background_colour_index: int
    pattern_colour_index: int

class XFAlignment(BaseObject, EqNeAttrs):
    hor_align: int
    vert_align: int
    rotation: int
    text_wrapped: Literal[0, 1]
    indent_level: int
    shrink_to_fit: Literal[0, 1]
    text_direction: Literal[0, 1, 2]

class XFProtection(BaseObject, EqNeAttrs):
    cell_locked: Literal[0, 1]
    formula_hidden: Literal[0, 1]

class XF(BaseObject):
    is_style: Literal[0, 1]
    parent_style_index: int
    xf_index: int
    font_index: int
    format_key: int
    protection: XFProtection | None
    background: XFBackground | None
    alignment: XFAlignment | None
    border: XFBorder | None
