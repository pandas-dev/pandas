from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum, Flag, IntEnum, IntFlag
from typing import Final, Literal
from typing_extensions import Self, TypeAlias

from .drawing import DeviceCMYK, DeviceGray, DeviceRGB
from .syntax import Name

_Color: TypeAlias = str | int | Sequence[int] | DeviceCMYK | DeviceGray | DeviceRGB

class SignatureFlag(IntEnum):
    SIGNATURES_EXIST = 1
    APPEND_ONLY = 2

class CoerciveEnum(Enum):  # type: ignore[misc]  # Enum with no members
    @classmethod
    def coerce(cls, value: Self | str, case_sensitive: bool = False) -> Self: ...

class CoerciveIntEnum(IntEnum):  # type: ignore[misc]  # Enum with no members
    @classmethod
    def coerce(cls, value: Self | str | int) -> Self: ...

class CoerciveIntFlag(IntFlag):  # type: ignore[misc]  # Enum with no members
    @classmethod
    def coerce(cls, value: Self | str | int) -> Self: ...

class WrapMode(CoerciveEnum):
    WORD = "WORD"
    CHAR = "CHAR"

class CharVPos(CoerciveEnum):
    SUP = "SUP"
    SUB = "SUB"
    NOM = "NOM"
    DENOM = "DENOM"
    LINE = "LINE"

class Align(CoerciveEnum):
    C = "CENTER"
    X = "X_CENTER"
    L = "LEFT"
    R = "RIGHT"
    J = "JUSTIFY"

    @classmethod
    def coerce(cls, value: Self | str) -> Self: ...  # type: ignore[override]

_Align: TypeAlias = Align | Literal["CENTER", "X_CENTER", "LEFT", "RIGHT", "JUSTIFY"]  # noqa: Y047

class VAlign(CoerciveEnum):
    M = "MIDDLE"
    T = "TOP"
    B = "BOTTOM"

    @classmethod
    def coerce(cls, value: Self | str) -> Self: ...  # type: ignore[override]

class TextEmphasis(CoerciveIntFlag):
    NONE = 0
    B = 1
    I = 2
    U = 4
    S = 8

    @property
    def style(self) -> str: ...
    def add(self, value: TextEmphasis) -> TextEmphasis: ...
    def remove(self, value: TextEmphasis) -> TextEmphasis: ...

class MethodReturnValue(CoerciveIntFlag):
    PAGE_BREAK = 1
    LINES = 2
    HEIGHT = 4

class CellBordersLayout(CoerciveIntFlag):
    NONE = 0
    LEFT = 1
    RIGHT = 2
    TOP = 4
    BOTTOM = 8
    ALL = 15
    INHERIT = 16

@dataclass
class TableBorderStyle:
    thickness: float | None = None
    color: int | tuple[int, int, int] | None = None
    dash: float | None = None
    gap: float = 0.0
    phase: float = 0.0

    @staticmethod
    def from_bool(should_draw: TableBorderStyle | bool | None) -> TableBorderStyle: ...
    @property
    def dash_dict(self) -> dict[str, float | None]: ...
    def changes_stroke(self, pdf) -> bool: ...
    def should_render(self) -> bool: ...
    def get_change_stroke_commands(self, scale: float) -> list[str]: ...
    @staticmethod
    def get_line_command(x1: float, y1: float, x2: float, y2: float) -> list[str]: ...
    def get_draw_commands(self, pdf, x1: float, y1: float, x2: float, y2: float) -> list[str]: ...

@dataclass
class TableCellStyle:
    left: bool | TableBorderStyle = False
    bottom: bool | TableBorderStyle = False
    right: bool | TableBorderStyle = False
    top: bool | TableBorderStyle = False

    @staticmethod
    def get_change_fill_color_command(color: _Color | None) -> list[str]: ...
    def get_draw_commands(
        self, pdf, x1: float, y1: float, x2: float, y2: float, fill_color: _Color | None = None
    ) -> list[str]: ...
    def override_cell_border(self, cell_border: CellBordersLayout) -> Self: ...
    def draw_cell_border(self, pdf, x1: float, y1: float, x2: float, y2: float, fill_color: _Color | None = None) -> None: ...

class TableBordersLayout(ABC):
    ALL: Final[TableBordersLayoutAll]
    NONE: Final[TableBordersLayoutNone]
    INTERNAL: Final[TableBordersLayoutInternal]
    MINIMAL: Final[TableBordersLayoutMinimal]
    HORIZONTAL_LINES: Final[TableBordersLayoutHorizontalLines]
    NO_HORIZONTAL_LINES: Final[TableBordersLayoutNoHorizontalLines]
    SINGLE_TOP_LINE: Final[TableBordersLayoutSingleTopLine]
    @abstractmethod
    def cell_style_getter(
        self, row_idx: int, col_idx: int, col_pos: int, num_heading_rows: int, num_rows: int, num_col_idx: int, num_col_pos: int
    ) -> TableCellStyle: ...
    @classmethod
    def coerce(cls, value: Self | str) -> Self: ...

class TableBordersLayoutAll(TableBordersLayout):
    def cell_style_getter(
        self, row_idx: int, col_idx: int, col_pos: int, num_heading_rows: int, num_rows: int, num_col_idx: int, num_col_pos: int
    ) -> TableCellStyle: ...

class TableBordersLayoutNone(TableBordersLayout):
    def cell_style_getter(
        self, row_idx: int, col_idx: int, col_pos: int, num_heading_rows: int, num_rows: int, num_col_idx: int, num_col_pos: int
    ) -> TableCellStyle: ...

class TableBordersLayoutInternal(TableBordersLayout):
    def cell_style_getter(
        self, row_idx: int, col_idx: int, col_pos: int, num_heading_rows: int, num_rows: int, num_col_idx: int, num_col_pos: int
    ) -> TableCellStyle: ...

class TableBordersLayoutMinimal(TableBordersLayout):
    def cell_style_getter(
        self, row_idx: int, col_idx: int, col_pos: int, num_heading_rows: int, num_rows: int, num_col_idx: int, num_col_pos: int
    ) -> TableCellStyle: ...

class TableBordersLayoutHorizontalLines(TableBordersLayout):
    def cell_style_getter(
        self, row_idx: int, col_idx: int, col_pos: int, num_heading_rows: int, num_rows: int, num_col_idx: int, num_col_pos: int
    ) -> TableCellStyle: ...

class TableBordersLayoutNoHorizontalLines(TableBordersLayout):
    def cell_style_getter(
        self, row_idx: int, col_idx: int, col_pos: int, num_heading_rows: int, num_rows: int, num_col_idx: int, num_col_pos: int
    ) -> TableCellStyle: ...

class TableBordersLayoutSingleTopLine(TableBordersLayout):
    def cell_style_getter(
        self, row_idx: int, col_idx: int, col_pos: int, num_heading_rows: int, num_rows: int, num_col_idx: int, num_col_pos: int
    ) -> TableCellStyle: ...

class TableCellFillMode(CoerciveEnum):
    NONE = "NONE"
    ALL = "ALL"
    ROWS = "ROWS"
    COLUMNS = "COLUMNS"
    EVEN_ROWS = "EVEN_ROWS"
    EVEN_COLUMNS = "EVEN_COLUMNS"

    def should_fill_cell(self, i: int, j: int) -> bool: ...
    @classmethod
    def coerce(cls, value: Self | str) -> Self: ...  # type: ignore[override]

class TableSpan(CoerciveEnum):
    ROW = "ROW"
    COL = "COL"

class TableHeadingsDisplay(CoerciveIntEnum):
    NONE = 0
    ON_TOP_OF_EVERY_PAGE = 1

class RenderStyle(CoerciveEnum):
    D = "DRAW"
    F = "FILL"
    DF = "DRAW_FILL"
    @property
    def operator(self) -> str: ...
    @property
    def is_draw(self) -> bool: ...
    @property
    def is_fill(self) -> bool: ...
    @classmethod
    def coerce(cls, value: Self | str) -> Self: ...  # type: ignore[override]

class TextMode(CoerciveIntEnum):
    FILL = 0
    STROKE = 1
    FILL_STROKE = 2
    INVISIBLE = 3
    FILL_CLIP = 4
    STROKE_CLIP = 5
    FILL_STROKE_CLIP = 6
    CLIP = 7

class XPos(CoerciveEnum):
    LEFT = "LEFT"
    RIGHT = "RIGHT"
    START = "START"
    END = "END"
    WCONT = "WCONT"
    CENTER = "CENTER"
    LMARGIN = "LMARGIN"
    RMARGIN = "RMARGIN"

class YPos(CoerciveEnum):
    TOP = "TOP"
    LAST = "LAST"
    NEXT = "NEXT"
    TMARGIN = "TMARGIN"
    BMARGIN = "BMARGIN"

class Angle(CoerciveIntEnum):
    NORTH = 90
    EAST = 0
    SOUTH = 270
    WEST = 180
    NORTHEAST = 45
    SOUTHEAST = 315
    SOUTHWEST = 225
    NORTHWEST = 135

class PageLayout(CoerciveEnum):
    SINGLE_PAGE = Name("SinglePage")
    ONE_COLUMN = Name("OneColumn")
    TWO_COLUMN_LEFT = Name("TwoColumnLeft")
    TWO_COLUMN_RIGHT = Name("TwoColumnRight")
    TWO_PAGE_LEFT = Name("TwoPageLeft")
    TWO_PAGE_RIGHT = Name("TwoPageRight")

class PageMode(CoerciveEnum):
    USE_NONE = Name("UseNone")
    USE_OUTLINES = Name("UseOutlines")
    USE_THUMBS = Name("UseThumbs")
    FULL_SCREEN = Name("FullScreen")
    USE_OC = Name("UseOC")
    USE_ATTACHMENTS = Name("UseAttachments")

class TextMarkupType(CoerciveEnum):
    HIGHLIGHT = Name("Highlight")
    UNDERLINE = Name("Underline")
    SQUIGGLY = Name("Squiggly")
    STRIKE_OUT = Name("StrikeOut")

class BlendMode(CoerciveEnum):
    NORMAL = Name("Normal")
    MULTIPLY = Name("Multiply")
    SCREEN = Name("Screen")
    OVERLAY = Name("Overlay")
    DARKEN = Name("Darken")
    LIGHTEN = Name("Lighten")
    COLOR_DODGE = Name("ColorDodge")
    COLOR_BURN = Name("ColorBurn")
    HARD_LIGHT = Name("HardLight")
    SOFT_LIGHT = Name("SoftLight")
    DIFFERENCE = Name("Difference")
    EXCLUSION = Name("Exclusion")
    HUE = Name("Hue")
    SATURATION = Name("Saturation")
    COLOR = Name("Color")
    LUMINOSITY = Name("Luminosity")

class AnnotationFlag(CoerciveIntEnum):
    INVISIBLE = 1
    HIDDEN = 2
    PRINT = 4
    NO_ZOOM = 8
    NO_ROTATE = 16
    NO_VIEW = 32
    READ_ONLY = 64
    LOCKED = 128
    TOGGLE_NO_VIEW = 256
    LOCKED_CONTENTS = 512

class AnnotationName(CoerciveEnum):
    NOTE = Name("Note")
    COMMENT = Name("Comment")
    HELP = Name("Help")
    PARAGRAPH = Name("Paragraph")
    NEW_PARAGRAPH = Name("NewParagraph")
    INSERT = Name("Insert")

class FileAttachmentAnnotationName(CoerciveEnum):
    PUSH_PIN = Name("PushPin")
    GRAPH_PUSH_PIN = Name("GraphPushPin")
    PAPERCLIP_TAG = Name("PaperclipTag")

class IntersectionRule(CoerciveEnum):
    NONZERO = "nonzero"
    EVENODD = "evenodd"

class PathPaintRule(CoerciveEnum):
    STROKE = "S"
    FILL_NONZERO = "f"
    FILL_EVENODD = "f*"
    STROKE_FILL_NONZERO = "B"
    STROKE_FILL_EVENODD = "B*"
    DONT_PAINT = "n"
    AUTO = "auto"

class ClippingPathIntersectionRule(CoerciveEnum):
    NONZERO = "W"
    EVENODD = "W*"

class StrokeCapStyle(CoerciveIntEnum):
    BUTT = 0
    ROUND = 1
    SQUARE = 2

class StrokeJoinStyle(CoerciveIntEnum):
    MITER = 0
    ROUND = 1
    BEVEL = 2

class PDFStyleKeys(Enum):
    FILL_ALPHA = Name("ca")
    BLEND_MODE = Name("BM")
    STROKE_ALPHA = Name("CA")
    STROKE_ADJUSTMENT = Name("SA")
    STROKE_WIDTH = Name("LW")
    STROKE_CAP_STYLE = Name("LC")
    STROKE_JOIN_STYLE = Name("LJ")
    STROKE_MITER_LIMIT = Name("ML")
    STROKE_DASH_PATTERN = Name("D")

class Corner(CoerciveEnum):
    TOP_RIGHT = "TOP_RIGHT"
    TOP_LEFT = "TOP_LEFT"
    BOTTOM_RIGHT = "BOTTOM_RIGHT"
    BOTTOM_LEFT = "BOTTOM_LEFT"

class FontDescriptorFlags(Flag):
    FIXED_PITCH = 1
    SYMBOLIC = 4
    ITALIC = 64
    FORCE_BOLD = 262144

class AccessPermission(IntFlag):
    PRINT_LOW_RES = 4
    MODIFY = 8
    COPY = 16
    ANNOTATION = 32
    FILL_FORMS = 256
    COPY_FOR_ACCESSIBILITY = 512
    ASSEMBLE = 1024
    PRINT_HIGH_RES = 2048
    @classmethod
    def all(cls) -> int: ...
    @classmethod
    def none(cls) -> Literal[0]: ...

class EncryptionMethod(Enum):
    NO_ENCRYPTION = 0
    RC4 = 1
    AES_128 = 2
    AES_256 = 3

class TextDirection(CoerciveEnum):
    LTR = "LTR"
    RTL = "RTL"
    TTB = "TTB"
    BTT = "BTT"

class OutputIntentSubType(CoerciveEnum):
    PDFX = "GTS_PDFX"
    PDFA = "GTS_PDFA1"
    ISOPDF = "ISO_PDFE1"

class PageLabelStyle(CoerciveEnum):
    NUMBER = "D"
    UPPER_ROMAN = "R"
    LOWER_ROMAN = "r"
    UPPER_LETTER = "A"
    LOWER_LETTER = "a"
    NONE = None

class Duplex(CoerciveEnum):
    SIMPLEX = "Simplex"
    DUPLEX_FLIP_SHORT_EDGE = "DuplexFlipShortEdge"
    DUPLEX_FLIP_LONG_EDGE = "DuplexFlipLongEdge"

class PageBoundaries(CoerciveEnum):
    ART_BOX = "ArtBox"
    BLEED_BOX = "BleedBox"
    CROP_BOX = "CropBox"
    MEDIA_BOX = "MediaBox"
    TRIM_BOX = "TrimBox"

class PageOrientation(CoerciveEnum):
    PORTRAIT = "P"
    LANDSCAPE = "L"

    @classmethod
    def coerce(cls, value: Self | str) -> Self: ...  # type: ignore[override]

class PDFResourceType(Enum):
    EXT_G_STATE = "ExtGState"
    COLOR_SPACE = "ColorSpace"
    PATTERN = "Pattern"
    SHADDING = "Shading"
    X_OBJECT = "XObject"
    FONT = "Font"
    PROC_SET = "ProcSet"
    PROPERTIES = "Properties"
