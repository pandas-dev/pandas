import sys
from collections.abc import Callable
from typing import Any, Final, TextIO

from .timemachine import *

DEBUG: Final[int]

class XLRDError(Exception): ...

class BaseObject:
    _repr_these: list[str]
    def dump(self, f: TextIO | None = None, header: str | None = None, footer: str | None = None, indent: int = 0) -> None: ...

FUN: Final[int]
FDT: Final[int]
FNU: Final[int]
FGE: Final[int]
FTX: Final[int]
DATEFORMAT: Final[int]
NUMBERFORMAT: Final[int]
XL_CELL_EMPTY: Final[int]
XL_CELL_TEXT: Final[int]
XL_CELL_NUMBER: Final[int]
XL_CELL_DATE: Final[int]
XL_CELL_BOOLEAN: Final[int]
XL_CELL_ERROR: Final[int]
XL_CELL_BLANK: Final[int]
biff_text_from_num: Final[dict[int, str]]
error_text_from_code: Final[dict[int, str]]
BIFF_FIRST_UNICODE: Final[int]
XL_WORKBOOK_GLOBALS: Final[int]
WBKBLOBAL: Final[int]
XL_WORKBOOK_GLOBALS_4W: Final[int]
XL_WORKSHEET: Final[int]
WRKSHEET: Final[int]
XL_BOUNDSHEET_WORKSHEET: Final[int]
XL_BOUNDSHEET_CHART: Final[int]
XL_BOUNDSHEET_VB_MODULE: Final[int]
XL_ARRAY: Final[int]
XL_ARRAY2: Final[int]
XL_BLANK: Final[int]
XL_BLANK_B2: Final[int]
XL_BOF: Final[int]
XL_BOOLERR: Final[int]
XL_BOOLERR_B2: Final[int]
XL_BOUNDSHEET: Final[int]
XL_BUILTINFMTCOUNT: Final[int]
XL_CF: Final[int]
XL_CODEPAGE: Final[int]
XL_COLINFO: Final[int]
XL_COLUMNDEFAULT: Final[int]
XL_COLWIDTH: Final[int]
XL_CONDFMT: Final[int]
XL_CONTINUE: Final[int]
XL_COUNTRY: Final[int]
XL_DATEMODE: Final[int]
XL_DEFAULTROWHEIGHT: Final[int]
XL_DEFCOLWIDTH: Final[int]
XL_DIMENSION: Final[int]
XL_DIMENSION2: Final[int]
XL_EFONT: Final[int]
XL_EOF: Final[int]
XL_EXTERNNAME: Final[int]
XL_EXTERNSHEET: Final[int]
XL_EXTSST: Final[int]
XL_FEAT11: Final[int]
XL_FILEPASS: Final[int]
XL_FONT: Final[int]
XL_FONT_B3B4: Final[int]
XL_FORMAT: Final[int]
XL_FORMAT2: Final[int]
XL_FORMULA: Final[int]
XL_FORMULA3: Final[int]
XL_FORMULA4: Final[int]
XL_GCW: Final[int]
XL_HLINK: Final[int]
XL_QUICKTIP: Final[int]
XL_HORIZONTALPAGEBREAKS: Final[int]
XL_INDEX: Final[int]
XL_INTEGER: Final[int]
XL_IXFE: Final[int]
XL_LABEL: Final[int]
XL_LABEL_B2: Final[int]
XL_LABELRANGES: Final[int]
XL_LABELSST: Final[int]
XL_LEFTMARGIN: Final[int]
XL_TOPMARGIN: Final[int]
XL_RIGHTMARGIN: Final[int]
XL_BOTTOMMARGIN: Final[int]
XL_HEADER: Final[int]
XL_FOOTER: Final[int]
XL_HCENTER: Final[int]
XL_VCENTER: Final[int]
XL_MERGEDCELLS: Final[int]
XL_MSO_DRAWING: Final[int]
XL_MSO_DRAWING_GROUP: Final[int]
XL_MSO_DRAWING_SELECTION: Final[int]
XL_MULRK: Final[int]
XL_MULBLANK: Final[int]
XL_NAME: Final[int]
XL_NOTE: Final[int]
XL_NUMBER: Final[int]
XL_NUMBER_B2: Final[int]
XL_OBJ: Final[int]
XL_PAGESETUP: Final[int]
XL_PALETTE: Final[int]
XL_PANE: Final[int]
XL_PRINTGRIDLINES: Final[int]
XL_PRINTHEADERS: Final[int]
XL_RK: Final[int]
XL_ROW: Final[int]
XL_ROW_B2: Final[int]
XL_RSTRING: Final[int]
XL_SCL: Final[int]
XL_SHEETHDR: Final[int]
XL_SHEETPR: Final[int]
XL_SHEETSOFFSET: Final[int]
XL_SHRFMLA: Final[int]
XL_SST: Final[int]
XL_STANDARDWIDTH: Final[int]
XL_STRING: Final[int]
XL_STRING_B2: Final[int]
XL_STYLE: Final[int]
XL_SUPBOOK: Final[int]
XL_TABLEOP: Final[int]
XL_TABLEOP2: Final[int]
XL_TABLEOP_B2: Final[int]
XL_TXO: Final[int]
XL_UNCALCED: Final[int]
XL_UNKNOWN: Final[int]
XL_VERTICALPAGEBREAKS: Final[int]
XL_WINDOW2: Final[int]
XL_WINDOW2_B2: Final[int]
XL_WRITEACCESS: Final[int]
XL_WSBOOL: Final[int]
XL_XF: Final[int]
XL_XF2: Final[int]
XL_XF3: Final[int]
XL_XF4: Final[int]
boflen: Final[dict[int, int]]
bofcodes: Final[tuple[int, int, int, int]]
XL_FORMULA_OPCODES: Final[tuple[int, int, int]]

def is_cell_opcode(c: int) -> bool: ...
def upkbits(
    tgt_obj: object, src: int, manifest: list[tuple[int, int, str]], local_setattr: Callable[[Any, str, Any], None] = ...
) -> None: ...
def upkbitsL(
    tgt_obj: object,
    src: int,
    manifest: list[tuple[int, int, str]],
    local_setattr: Callable[[Any, str, Any], None] = ...,
    local_int: Callable[[Any], int] = ...,
) -> None: ...
def unpack_string(data: bytes, pos: int, encoding: str, lenlen: int = 1) -> str: ...
def unpack_string_update_pos(
    data: bytes, pos: int, encoding: str, lenlen: int = 1, known_len: int | None = None
) -> tuple[str, int]: ...
def unpack_unicode(data: bytes, pos: int, lenlen: int = 2) -> str: ...
def unpack_unicode_update_pos(data: bytes, pos: int, lenlen: int = 2, known_len: int | None = None) -> tuple[str, int]: ...
def unpack_cell_range_address_list_update_pos(
    output_list: list[tuple[int, int, int, int]], data: bytes, pos: int, biff_version: int, addr_size: int = 6
) -> int: ...

biff_rec_name_dict: Final[dict[int, str]]

def hex_char_dump(
    strg: bytes, ofs: int, dlen: int, base: int = 0, fout: TextIO = sys.stdout, unnumbered: bool = False
) -> None: ...
def biff_dump(
    mem: bytes, stream_offset: int, stream_len: int, base: int = 0, fout: TextIO = sys.stdout, unnumbered: bool = False
) -> None: ...
def biff_count_records(mem: bytes, stream_offset: int, stream_len: int, fout: TextIO = sys.stdout) -> None: ...

encoding_from_codepage: Final[dict[int, str]]
