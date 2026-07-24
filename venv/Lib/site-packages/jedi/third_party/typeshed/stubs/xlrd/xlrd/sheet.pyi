from _typeshed import SupportsWrite
from array import array
from collections.abc import Callable, Generator, Sequence
from typing import Any, Final, Literal, overload

from .biffh import *
from .book import Book
from .formatting import XF
from .timemachine import *

OBJ_MSO_DEBUG: Final[int]

class MSODrawing(BaseObject): ...
class MSObj(BaseObject): ...
class MSTxo(BaseObject): ...

class Note(BaseObject):
    author: str
    col_hidden: int
    colx: int
    rich_text_runlist: list[tuple[str, int]] | None
    row_hidden: int
    rowx: int
    show: int
    text: str

class Hyperlink(BaseObject):
    frowx: int | None
    lrowx: int | None
    fcolx: int | None
    lcolx: int | None
    type: str | None
    url_or_path: bytes | str | None
    desc: str | None
    target: str | None
    textmark: str | None
    quicktip: str | None

def unpack_RK(rk_str: bytes) -> float: ...

cellty_from_fmtty: Final[dict[int, int]]
ctype_text: Final[dict[int, str]]

class Cell(BaseObject):
    __slots__ = ["ctype", "value", "xf_index"]
    ctype: Literal[0, 1, 2, 3, 4, 5, 6]
    value: str | float
    xf_index: int | None
    def __init__(self, ctype: Literal[0, 1, 2, 3, 4, 5, 6], value: str, xf_index: int | None = None) -> None: ...

empty_cell: Final[Cell]

class Colinfo(BaseObject):
    width: int
    xf_index: int
    hidden: int
    bit1_flag: int
    outline_level: int
    collapsed: int

class Rowinfo(BaseObject):
    __slots__ = (
        "height",
        "has_default_height",
        "outline_level",
        "outline_group_starts_ends",
        "hidden",
        "height_mismatch",
        "has_default_xf_index",
        "xf_index",
        "additional_space_above",
        "additional_space_below",
    )
    height: int | None
    has_default_height: int | None
    outline_level: int | None
    outline_group_starts_ends: int | None
    hidden: int | None
    height_mismatch: int | None
    has_default_xf_index: int | None
    xf_index: int | None
    additional_space_above: int | None
    additional_space_below: int | None
    def __init__(self) -> None: ...
    def __getstate__(self) -> tuple[int | None, ...]: ...
    def __setstate__(self, state: tuple[int | None, ...]) -> None: ...

class Sheet(BaseObject):
    name: str
    book: Book | None
    nrows: int
    ncols: int
    colinfo_map: dict[int, Colinfo]
    rowinfo_map: dict[int, Rowinfo]
    col_label_ranges: list[tuple[int, int, int, int]]
    row_label_ranges: list[tuple[int, int, int, int]]
    merged_cells: list[tuple[int, int, int, int]]
    rich_text_runlist_map: dict[tuple[int, int], list[tuple[int, int]]]
    defcolwidth: float | None
    standardwidth: float | None
    default_row_height: int | None
    default_row_height_mismatch: int | None
    default_row_hidden: int | None
    default_additional_space_above: int | None
    default_additional_space_below: int | None
    visibility: Literal[0, 1, 2]
    gcw: tuple[int, ...]
    hyperlink_list: list[Hyperlink]
    hyperlink_map: dict[tuple[int, int], Hyperlink]
    cell_note_map: dict[tuple[int, int], Note]
    vert_split_pos: int
    horz_split_pos: int
    horz_split_first_visible: int
    vert_split_first_visible: int
    split_active_pane: int
    has_pane_record: int
    horizontal_page_breaks: list[tuple[int, int, int]]
    vertical_page_breaks: list[tuple[int, int, int]]
    biff_version: int
    logfile: SupportsWrite[str]
    bt: array[int]
    bf: array[int]
    number: int
    verbosity: int
    formatting_info: bool
    ragged_rows: bool
    put_cell: Callable[[int, int, int | None, str, int | None], None]
    first_visible_rowx: int
    first_visible_colx: int
    gridline_colour_index: int
    gridline_colour_rgb: tuple[int, int, int] | None
    cooked_page_break_preview_mag_factor: int
    cooked_normal_view_mag_factor: int
    cached_page_break_preview_mag_factor: int
    cached_normal_view_mag_factor: int
    scl_mag_factor: int | None
    utter_max_rows: int
    utter_max_cols: int
    def __init__(self, book: Book, position: int, name: str, number: int) -> None: ...
    def cell(self, rowx: int, colx: int) -> Cell: ...
    def cell_value(self, rowx: int, colx: int) -> str: ...
    def cell_type(self, rowx: int, colx: int) -> int: ...
    def cell_xf_index(self, rowx: int, colx: int) -> int: ...
    def row_len(self, rowx: int) -> int: ...
    def row(self, rowx: int) -> list[Cell]: ...
    @overload
    def __getitem__(self, item: int) -> list[Cell]: ...
    @overload
    def __getitem__(self, item: tuple[int, int]) -> Cell: ...
    def get_rows(self) -> Generator[list[Cell]]: ...
    __iter__ = get_rows
    def row_types(self, rowx: int, start_colx: int = 0, end_colx: int | None = None) -> Sequence[int]: ...
    def row_values(self, rowx: int, start_colx: int = 0, end_colx: int | None = None) -> Sequence[str]: ...
    def row_slice(self, rowx: int, start_colx: int = 0, end_colx: int | None = None) -> list[Cell]: ...
    def col_slice(self, colx: int, start_rowx: int = 0, end_rowx: int | None = None) -> list[Cell]: ...
    def col_values(self, colx: int, start_rowx: int = 0, end_rowx: int | None = None) -> list[str]: ...
    def col_types(self, colx: int, start_rowx: int = 0, end_rowx: int | None = None) -> list[int]: ...
    col = col_slice
    def tidy_dimensions(self) -> None: ...
    def put_cell_ragged(self, rowx: int, colx: int, ctype: int | None, value: str, xf_index: int | None) -> None: ...
    def put_cell_unragged(self, rowx: int, colx: int, ctype: int | None, value: str, xf_index: int | None) -> None: ...
    def read(self, bk: Book) -> Literal[1]: ...
    def string_record_contents(self, data: bytes) -> str | None: ...
    def update_cooked_mag_factors(self) -> None: ...
    def fixed_BIFF2_xfindex(self, cell_attr: bytes, rowx: int, colx: int, true_xfx: int | None = None) -> int: ...
    def insert_new_BIFF20_xf(self, cell_attr: bytes, style: int = 0) -> int: ...
    def fake_XF_from_BIFF20_cell_attr(self, cell_attr: bytes, style: int = 0) -> XF: ...
    def req_fmt_info(self) -> None: ...
    def computed_column_width(self, colx: int) -> float: ...
    def handle_hlink(self, data: bytes) -> None: ...
    def handle_quicktip(self, data: bytes) -> None: ...
    def handle_msodrawingetc(self, recid: Any, data_len: int, data: bytes) -> None: ...
    def handle_obj(self, data: bytes) -> MSObj | None: ...
    def handle_note(self, data: bytes, txos: dict[int, MSTxo]) -> None: ...
    def handle_txo(self, data: bytes) -> MSTxo | None: ...
    def handle_feat11(self, data: bytes) -> None: ...
