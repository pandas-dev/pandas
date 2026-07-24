from typing import Final

from .book import Book, Name
from .timemachine import *

__all__ = [
    "oBOOL",
    "oERR",
    "oNUM",
    "oREF",
    "oREL",
    "oSTRG",
    "oUNK",
    "decompile_formula",
    "dump_formula",
    "evaluate_name_formula",
    "okind_dict",
    "rangename3d",
    "rangename3drel",
    "cellname",
    "cellnameabs",
    "colname",
    "FMLA_TYPE_CELL",
    "FMLA_TYPE_SHARED",
    "FMLA_TYPE_ARRAY",
    "FMLA_TYPE_COND_FMT",
    "FMLA_TYPE_DATA_VAL",
    "FMLA_TYPE_NAME",
    "Operand",
    "Ref3D",
]

FMLA_TYPE_CELL: Final[int]
FMLA_TYPE_SHARED: Final[int]
FMLA_TYPE_ARRAY: Final[int]
FMLA_TYPE_COND_FMT: Final[int]
FMLA_TYPE_DATA_VAL: Final[int]
FMLA_TYPE_NAME: Final[int]
oBOOL: Final[int]
oERR: Final[int]
oNUM: Final[int]
oREF: Final[int]
oREL: Final[int]
oSTRG: Final[int]
oUNK: Final[int]
okind_dict: Final[dict[int, str]]

class FormulaError(Exception): ...

class Operand:
    value: float | str | None
    kind: int
    text: str
    rank: int
    def __init__(self, akind: int | None = None, avalue: float | str | None = None, arank: int = 0, atext: str = "?") -> None: ...

class Ref3D(tuple[int, int, int, int, int, int, int, int, int, int, int, int]):
    coords: tuple[int, int, int, int, int, int]
    relflags: tuple[int, int, int, int, int, int]
    shtxlo: int
    shtxhi: int
    rowxlo: int
    rowxhi: int
    colxlo: int
    colxhi: int
    def __init__(  # pyright: ignore[reportInconsistentConstructor]
        self, atuple: tuple[int, int, int, int, int, int, int, int, int, int, int, int]
    ) -> None: ...

def evaluate_name_formula(bk: Book, nobj: Name, namex: str, blah: int = 0, level: int = 0) -> None: ...
def decompile_formula(
    bk: Book,
    fmla: bytes,
    fmlalen: int,
    fmlatype: int | None = None,
    browx: int | None = None,
    bcolx: int | None = None,
    blah: int = 0,
    level: int = 0,
    r1c1: int = 0,
) -> str | None: ...
def dump_formula(bk: Book, data: bytes, fmlalen: int, bv: int, reldelta: int, blah: int = 0, isname: int = 0) -> None: ...
def cellname(rowx: int, colx: int) -> str: ...
def cellnameabs(rowx: int, colx: int, r1c1: int = 0) -> str: ...
def colname(colx: int) -> str: ...
def rangename3d(book: Book, ref3d: Ref3D) -> str: ...
def rangename3drel(book: Book, ref3d: Ref3D, browx: int | None = None, bcolx: int | None = None, r1c1: int = 0) -> str: ...
