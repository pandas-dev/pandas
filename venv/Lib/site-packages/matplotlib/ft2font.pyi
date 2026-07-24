from enum import Enum, Flag
from os import PathLike
import sys
from typing import BinaryIO, Literal, NewType, NotRequired, TypeAlias, TypedDict, cast, final, overload
from typing_extensions import Buffer  # < Py 3.12

import numpy as np
from numpy.typing import NDArray

__freetype_build_type__: str
__freetype_version__: str
__libraqm_version__: str

# We can't change the type hints for standard library chr/ord, so character codes are a
# simple type alias.
CharacterCodeType: TypeAlias = int
# But glyph indices are internal, so use a distinct type hint.
GlyphIndexType = NewType('GlyphIndexType', int)

class FaceFlags(Flag):
    SCALABLE = cast(int, ...)
    FIXED_SIZES = cast(int, ...)
    FIXED_WIDTH = cast(int, ...)
    SFNT = cast(int, ...)
    HORIZONTAL = cast(int, ...)
    VERTICAL = cast(int, ...)
    KERNING = cast(int, ...)
    FAST_GLYPHS = cast(int, ...)
    MULTIPLE_MASTERS = cast(int, ...)
    GLYPH_NAMES = cast(int, ...)
    EXTERNAL_STREAM = cast(int, ...)
    HINTER = cast(int, ...)
    CID_KEYED = cast(int, ...)
    TRICKY = cast(int, ...)
    COLOR = cast(int, ...)
    VARIATION = cast(int, ...)
    SVG = cast(int, ...)
    SBIX = cast(int, ...)
    SBIX_OVERLAY = cast(int, ...)

class Kerning(Enum):
    DEFAULT = cast(int, ...)
    UNFITTED = cast(int, ...)
    UNSCALED = cast(int, ...)

class LoadFlags(Flag):
    DEFAULT = cast(int, ...)
    NO_SCALE = cast(int, ...)
    NO_HINTING = cast(int, ...)
    RENDER = cast(int, ...)
    NO_BITMAP = cast(int, ...)
    VERTICAL_LAYOUT = cast(int, ...)
    FORCE_AUTOHINT = cast(int, ...)
    CROP_BITMAP = cast(int, ...)
    PEDANTIC = cast(int, ...)
    IGNORE_GLOBAL_ADVANCE_WIDTH = cast(int, ...)
    NO_RECURSE = cast(int, ...)
    IGNORE_TRANSFORM = cast(int, ...)
    MONOCHROME = cast(int, ...)
    LINEAR_DESIGN = cast(int, ...)
    NO_AUTOHINT = cast(int, ...)
    COLOR = cast(int, ...)
    COMPUTE_METRICS = cast(int, ...)
    BITMAP_METRICS_ONLY = cast(int, ...)
    NO_SVG = cast(int, ...)
    # The following should be unique, but the above can be OR'd together.
    TARGET_NORMAL = cast(int, ...)
    TARGET_LIGHT = cast(int, ...)
    TARGET_MONO = cast(int, ...)
    TARGET_LCD = cast(int, ...)
    TARGET_LCD_V = cast(int, ...)

class RenderMode(Enum):
    NORMAL = cast(int, ...)
    LIGHT = cast(int, ...)
    MONO = cast(int, ...)
    LCD = cast(int, ...)
    LCD_V = cast(int, ...)
    SDF = cast(int, ...)

class StyleFlags(Flag):
    NORMAL = cast(int, ...)
    ITALIC = cast(int, ...)
    BOLD = cast(int, ...)

class _SfntHeadDict(TypedDict):
    version: tuple[int, int]
    fontRevision: tuple[int, int]
    checkSumAdjustment: int
    magicNumber: int
    flags: int
    unitsPerEm: int
    created: tuple[int, int]
    modified: tuple[int, int]
    xMin: int
    yMin: int
    xMax: int
    yMax: int
    macStyle: int
    lowestRecPPEM: int
    fontDirectionHint: int
    indexToLocFormat: int
    glyphDataFormat: int

class _SfntMaxpDict(TypedDict):
    version: tuple[int, int]
    numGlyphs: int
    maxPoints: int
    maxContours: int
    maxComponentPoints: int
    maxComponentContours: int
    maxZones: int
    maxTwilightPoints: int
    maxStorage: int
    maxFunctionDefs: int
    maxInstructionDefs: int
    maxStackElements: int
    maxSizeOfInstructions: int
    maxComponentElements: int
    maxComponentDepth: int

class _SfntOs2Dict(TypedDict):
    version: int
    xAvgCharWidth: int
    usWeightClass: int
    usWidthClass: int
    fsType: int
    ySubscriptXSize: int
    ySubscriptYSize: int
    ySubscriptXOffset: int
    ySubscriptYOffset: int
    ySuperscriptXSize: int
    ySuperscriptYSize: int
    ySuperscriptXOffset: int
    ySuperscriptYOffset: int
    yStrikeoutSize: int
    yStrikeoutPosition: int
    sFamilyClass: int
    panose: bytes
    ulUnicodeRange: tuple[int, int, int, int]
    achVendID: bytes
    fsSelection: int
    usFirstCharIndex: int
    usLastCharIndex: int
    sTypoAscender: int
    sTypoDescender: int
    sTypoLineGap: int
    usWinAscent: int
    usWinDescent: int
    # version >= 1
    ulCodePageRange: NotRequired[tuple[int, int]]
    # version >= 2
    sxHeight: NotRequired[int]
    sCapHeight: NotRequired[int]
    usDefaultChar: NotRequired[int]
    usBreakChar: NotRequired[int]
    usMaxContext: NotRequired[int]
    # version >= 5
    usLowerOpticalPointSize: NotRequired[int]
    usUpperOpticalPointSize: NotRequired[int]

class _SfntHheaDict(TypedDict):
    version: tuple[int, int]
    ascent: int
    descent: int
    lineGap: int
    advanceWidthMax: int
    minLeftBearing: int
    minRightBearing: int
    xMaxExtent: int
    caretSlopeRise: int
    caretSlopeRun: int
    caretOffset: int
    metricDataFormat: int
    numOfLongHorMetrics: int

class _SfntVheaDict(TypedDict):
    version: tuple[int, int]
    vertTypoAscender: int
    vertTypoDescender: int
    vertTypoLineGap: int
    advanceHeightMax: int
    minTopSideBearing: int
    minBottomSideBearing: int
    yMaxExtent: int
    caretSlopeRise: int
    caretSlopeRun: int
    caretOffset: int
    metricDataFormat: int
    numOfLongVerMetrics: int

class _SfntPostDict(TypedDict):
    format: tuple[int, int]
    italicAngle: tuple[int, int]
    underlinePosition: int
    underlineThickness: int
    isFixedPitch: int
    minMemType42: int
    maxMemType42: int
    minMemType1: int
    maxMemType1: int

class _SfntPcltDict(TypedDict):
    version: tuple[int, int]
    fontNumber: int
    pitch: int
    xHeight: int
    style: int
    typeFamily: int
    capHeight: int
    symbolSet: int
    typeFace: bytes
    characterComplement: bytes
    strokeWeight: int
    widthType: int
    serifStyle: int

@final
class LayoutItem:
    @property
    def ft_object(self) -> FT2Font: ...
    @property
    def char(self) -> str: ...
    @property
    def glyph_index(self) -> GlyphIndexType: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    @property
    def prev_kern(self) -> float: ...
    def __str__(self) -> str: ...

@final
class FT2Font(Buffer):
    def __init__(
        self,
        filename: str | bytes | PathLike | BinaryIO,
        *,
        face_index: int = ...,
        _fallback_list: list[FT2Font] | None = ...,
        _kerning_factor: int | None = ...,
        _warn_if_used: bool = ...,
    ) -> None: ...
    if sys.version_info[:2] >= (3, 12):
        def __buffer__(self, /, flags: int) -> memoryview: ...
    def _layout(
        self,
        text: str,
        flags: LoadFlags,
        features: tuple[str, ...] | None = ...,
        language: str | tuple[tuple[str, int, int], ...] | None = ...,
    ) -> list[LayoutItem]: ...
    def clear(self) -> None: ...
    def draw_glyph_to_bitmap(
        self, image: NDArray[np.uint8], x: int, y: int, glyph: Glyph, antialiased: bool = ...
    ) -> None: ...
    def draw_glyphs_to_bitmap(self, antialiased: bool = ...) -> None: ...
    def get_bitmap_offset(self) -> tuple[int, int]: ...
    def get_char_index(self, codepoint: CharacterCodeType) -> GlyphIndexType: ...
    def get_charmap(self) -> dict[CharacterCodeType, GlyphIndexType]: ...
    def get_descent(self) -> int: ...
    def get_glyph_name(self, index: GlyphIndexType) -> str: ...
    def get_image(self) -> NDArray[np.uint8]: ...
    def get_kerning(self, left: GlyphIndexType, right: GlyphIndexType, mode: Kerning) -> int: ...
    def get_name_index(self, name: str) -> GlyphIndexType: ...
    def get_num_glyphs(self) -> int: ...
    def get_path(self) -> tuple[NDArray[np.float64], NDArray[np.int8]]: ...
    def get_ps_font_info(
        self,
    ) -> tuple[str, str, str, str, str, int, int, int, int]: ...
    def get_sfnt(self) -> dict[tuple[int, int, int, int], bytes]: ...
    @overload
    def get_sfnt_table(self, name: Literal["head"]) -> _SfntHeadDict | None: ...
    @overload
    def get_sfnt_table(self, name: Literal["maxp"]) -> _SfntMaxpDict | None: ...
    @overload
    def get_sfnt_table(self, name: Literal["OS/2"]) -> _SfntOs2Dict | None: ...
    @overload
    def get_sfnt_table(self, name: Literal["hhea"]) -> _SfntHheaDict | None: ...
    @overload
    def get_sfnt_table(self, name: Literal["vhea"]) -> _SfntVheaDict | None: ...
    @overload
    def get_sfnt_table(self, name: Literal["post"]) -> _SfntPostDict | None: ...
    @overload
    def get_sfnt_table(self, name: Literal["pclt"]) -> _SfntPcltDict | None: ...
    def get_width_height(self) -> tuple[int, int]: ...
    def load_char(self, charcode: CharacterCodeType, flags: LoadFlags = ...) -> Glyph: ...
    def load_glyph(self, glyphindex: GlyphIndexType, flags: LoadFlags = ...) -> Glyph: ...
    def select_charmap(self, i: int) -> None: ...
    def set_charmap(self, i: int) -> None: ...
    def set_size(self, ptsize: float, dpi: float) -> None: ...
    def set_text(
        self,
        string: str,
        angle: float = ...,
        flags: LoadFlags = ...,
        *,
        features: tuple[str] | None = ...,
        language: str | list[tuple[str, int, int]] | None = ...,
    ) -> NDArray[np.float64]: ...
    @property
    def ascender(self) -> int: ...
    @property
    def bbox(self) -> tuple[int, int, int, int]: ...
    @property
    def descender(self) -> int: ...
    @property
    def face_flags(self) -> FaceFlags: ...
    @property
    def face_index(self) -> int: ...
    @property
    def family_name(self) -> str: ...
    @property
    def fname(self) -> str | bytes: ...
    @property
    def height(self) -> int: ...
    @property
    def max_advance_height(self) -> int: ...
    @property
    def max_advance_width(self) -> int: ...
    @property
    def num_charmaps(self) -> int: ...
    @property
    def num_faces(self) -> int: ...
    @property
    def num_fixed_sizes(self) -> int: ...
    @property
    def num_glyphs(self) -> int: ...
    @property
    def num_named_instances(self) -> int: ...
    @property
    def postscript_name(self) -> str: ...
    @property
    def scalable(self) -> bool: ...
    @property
    def style_flags(self) -> StyleFlags: ...
    @property
    def style_name(self) -> str: ...
    @property
    def underline_position(self) -> int: ...
    @property
    def underline_thickness(self) -> int: ...
    @property
    def units_per_EM(self) -> int: ...

@final
class FT2Image(Buffer):
    def __init__(self, width: int, height: int) -> None: ...
    def draw_rect_filled(self, x0: int, y0: int, x1: int, y1: int) -> None: ...
    if sys.version_info[:2] >= (3, 12):
        def __buffer__(self, /, flags: int) -> memoryview: ...

@final
class Glyph:
    @property
    def width(self) -> int: ...
    @property
    def height(self) -> int: ...
    @property
    def horiBearingX(self) -> int: ...
    @property
    def horiBearingY(self) -> int: ...
    @property
    def horiAdvance(self) -> int: ...
    @property
    def linearHoriAdvance(self) -> int: ...
    @property
    def vertBearingX(self) -> int: ...
    @property
    def vertBearingY(self) -> int: ...
    @property
    def vertAdvance(self) -> int: ...
    @property
    def bbox(self) -> tuple[int, int, int, int]: ...
