from _typeshed import Incomplete, ReadableBuffer, StrOrBytesPath
from collections.abc import Sequence
from typing import Final, Literal, NamedTuple
from typing_extensions import Self
from weakref import WeakKeyDictionary

from reportlab.pdfbase import pdfdoc, pdfmetrics

__version__: Final[str]

class TTFError(pdfdoc.PDFError): ...

def SUBSETN(n, table: ReadableBuffer | None = ...) -> bytes: ...
def makeToUnicodeCMap(fontname: str, subset) -> str: ...
def splice(stream, offset, value): ...

GF_ARG_1_AND_2_ARE_WORDS: Final = 1
GF_ARGS_ARE_XY_VALUES: Final = 2
GF_ROUND_XY_TO_GRID: Final = 4
GF_WE_HAVE_A_SCALE: Final = 8
GF_RESERVED: Final = 16
GF_MORE_COMPONENTS: Final = 32
GF_WE_HAVE_AN_X_AND_Y_SCALE: Final = 64
GF_WE_HAVE_A_TWO_BY_TWO: Final = 128
GF_WE_HAVE_INSTRUCTIONS: Final = 256
GF_USE_MY_METRICS: Final = 512
GF_OVERLAP_COMPOUND: Final = 1024
GF_SCALED_COMPONENT_OFFSET: Final = 2048
GF_UNSCALED_COMPONENT_OFFSET: Final = 4096

def TTFOpenFile(fn: StrOrBytesPath) -> tuple[StrOrBytesPath,]: ...

class TTFontParser:
    ttfVersions: tuple[int, ...]
    ttcVersions: tuple[int, ...]
    fileKind: str
    validate: bool | Literal[0, 1]
    subfontNameX: bytes
    def __init__(self, file, validate: bool | Literal[0, 1] = 0, subfontIndex: int = 0) -> None: ...
    ttcVersion: int
    numSubfonts: int
    subfontOffsets: list[int]
    def readTTCHeader(self) -> None: ...
    def getSubfont(self, subfontIndex: int) -> None: ...
    numTables: int
    searchRange: int
    entrySelector: int
    rangeShift: int
    table: dict[Incomplete, Incomplete]
    tables: list[Incomplete]
    def readTableDirectory(self) -> None: ...
    version: int
    def readHeader(self) -> bool: ...
    filename: Incomplete
    def readFile(self, f) -> None: ...
    def checksumTables(self) -> None: ...
    def checksumFile(self) -> None: ...
    def get_table_pos(self, tag) -> tuple[Incomplete, Incomplete]: ...
    def seek(self, pos: int) -> None: ...
    def skip(self, delta: int) -> None: ...
    def seek_table(self, tag, offset_in_table: int = 0) -> int: ...
    def read_tag(self) -> str: ...
    def get_chunk(self, pos: int, length: int) -> bytes: ...
    def read_uint8(self) -> int: ...
    def read_ushort(self) -> int: ...
    def read_ulong(self) -> int: ...
    def read_short(self) -> int: ...
    def get_ushort(self, pos: int) -> int: ...
    def get_ulong(self, pos: int) -> int: ...
    def get_table(self, tag): ...

class TTFontMaker:
    tables: dict[Incomplete, Incomplete]
    def __init__(self) -> None: ...
    def add(self, tag, data) -> None: ...
    def makeStream(self) -> bytes: ...

class CMapFmt2SubHeader(NamedTuple):
    firstCode: int
    entryCount: int
    idDelta: int
    idRangeOffset: int

class TTFNameBytes(bytes):
    ustr: Incomplete
    def __new__(cls, b, enc: str = "utf8") -> Self: ...

class TTFontFile(TTFontParser):
    def __init__(
        self, file, charInfo: bool | Literal[0, 1] = 1, validate: bool | Literal[0, 1] = 0, subfontIndex: int | str | bytes = 0
    ) -> None: ...
    name: Incomplete
    familyName: Incomplete
    styleName: Incomplete
    fullName: Incomplete
    uniqueFontID: Incomplete
    fontRevision: Incomplete
    unitsPerEm: Incomplete
    bbox: Incomplete
    ascent: Incomplete
    descent: Incomplete
    capHeight: Incomplete
    stemV: Incomplete
    italicAngle: Incomplete
    underlinePosition: Incomplete
    underlineThickness: Incomplete
    flags: Incomplete
    numGlyphs: Incomplete
    charToGlyph: Incomplete
    defaultWidth: Incomplete
    charWidths: Incomplete
    hmetrics: Incomplete
    glyphPos: Incomplete
    def extractInfo(self, charInfo: bool | Literal[0, 1] = 1) -> None: ...
    def makeSubset(self, subset: Sequence[Incomplete]) -> bytes: ...

FF_FIXED: Final = 1
FF_SERIF: Final = 2
FF_SYMBOLIC: Final = 4
FF_SCRIPT: Final = 8
FF_NONSYMBOLIC: Final = 32
FF_ITALIC: Final = 64
FF_ALLCAP: Final = 65536
FF_SMALLCAP: Final = 131072
FF_FORCEBOLD: Final = 262144

class TTFontFace(TTFontFile, pdfmetrics.TypeFace):
    def __init__(self, filename, validate: bool | Literal[0, 1] = 0, subfontIndex: int | str | bytes = 0) -> None: ...
    def getCharWidth(self, code): ...
    def addSubsetObjects(self, doc, fontname, subset): ...

class TTEncoding:
    name: str
    def __init__(self) -> None: ...

class TTFont:
    class State:
        namePrefix: str
        nextCode: int
        internalName: Incomplete
        frozen: bool | Literal[0, 1]
        subsets: Incomplete
        def __init__(self, asciiReadable: bool | Literal[0, 1] | None = None, ttf=None) -> None: ...

    fontName: str
    face: TTFontFace
    encoding: TTEncoding
    state: WeakKeyDictionary[Incomplete, State]
    def __init__(
        self,
        name: str,
        filename,
        validate: bool | Literal[0, 1] = 0,
        subfontIndex: int | str | bytes = 0,
        asciiReadable: bool | Literal[0, 1] | None = None,
        shapable: bool = True,
    ) -> None: ...
    def stringWidth(self, text, size, encoding: str = "utf8") -> float: ...
    def splitString(self, text, doc, encoding: str = "utf-8") -> list[tuple[int, bytes]]: ...
    def getSubsetInternalName(self, subset, doc) -> str: ...
    def addObjects(self, doc) -> None: ...
    @property
    def hbFace(self) -> Incomplete | None: ...
    def hbFont(self, fontSize: float = 10): ...
    @property
    def shapable(self) -> bool: ...
    @shapable.setter
    def shapable(self, v) -> None: ...
    def pdfScale(self, v): ...
    def unregister(self) -> None: ...

class ShapedFragWord(list[Incomplete]): ...

class ShapeData(NamedTuple):
    cluster: int
    x_advance: float
    y_advance: float
    x_offset: float
    y_offset: float
    width: float

class ShapedStr(str):
    def __new__(cls, s, shapeData: ShapeData | None = None) -> Self: ...
    def __add__(self, other) -> ShapedStr: ...
    def __radd__(self, other) -> ShapedStr: ...

def shapeStr(s: str, fontName: str, fontSize: float, force: bool = False): ...
def freshTTFont(ttfn, ttfpath, **kwds) -> TTFont: ...
def makeShapedFragWord(w, K: list[Incomplete] = [], V: list[Incomplete] = []) -> type[ShapedFragWord]: ...
def shapeFragWord(w, features=None, force: bool = False): ...
