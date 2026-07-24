import dataclasses
from _typeshed import Incomplete, Unused
from collections import defaultdict
from collections.abc import Generator
from dataclasses import dataclass
from logging import Logger
from typing import Final, overload
from typing_extensions import Self, deprecated

from ._fonttools_shims import _TTFont
from .drawing import DeviceGray, DeviceRGB, Number
from .enums import Align, TextEmphasis
from .syntax import PDFObject

LOGGER: Logger

# Only defined if harfbuzz is installed.
class HarfBuzzFont(Incomplete):  # derives from uharfbuzz.Font
    def __deepcopy__(self, _memo: object) -> Self: ...

@dataclass
class FontFace:
    __slots__ = ("family", "emphasis", "size_pt", "color", "fill_color")
    family: str | None
    emphasis: TextEmphasis | None
    size_pt: int | None
    color: DeviceGray | DeviceRGB | None
    fill_color: DeviceGray | DeviceRGB | None

    def __init__(
        self,
        family: str | None = None,
        emphasis=None,
        size_pt: int | None = None,
        color: int | tuple[Number, Number, Number] | DeviceGray | DeviceRGB | None = None,
        fill_color: int | tuple[Number, Number, Number] | DeviceGray | DeviceRGB | None = None,
    ) -> None: ...

    replace = dataclasses.replace

    @overload
    @staticmethod
    def combine(default_style: None, override_style: None) -> None: ...  # type: ignore[misc]
    @overload
    @staticmethod
    def combine(default_style: FontFace | None, override_style: FontFace | None) -> FontFace: ...

class TextStyle(FontFace):
    t_margin: int
    l_margin: int | Align
    b_margin: int
    def __init__(
        self,
        font_family: str | None = None,
        font_style: str | None = None,
        font_size_pt: int | None = None,
        color: int | tuple[int, int, int] | None = None,
        fill_color: int | tuple[int, int, int] | None = None,
        underline: bool = False,
        t_margin: int | None = None,
        l_margin: int | Align | str | None = None,
        b_margin: int | None = None,
    ): ...
    def replace(  # type: ignore[override]
        self,
        /,
        font_family: str | None = None,
        emphasis: TextEmphasis | None = None,
        font_size_pt: int | None = None,
        color: int | tuple[int, int, int] | None = None,
        fill_color: int | tuple[int, int, int] | None = None,
        t_margin: int | None = None,
        l_margin: int | None = None,
        b_margin: int | None = None,
    ) -> TextStyle: ...

@deprecated("fpdf.TitleStyle is deprecated since 2.7.10. It has been replaced by fpdf.TextStyle.")
class TitleStyle(TextStyle): ...

__pdoc__: Final[dict[str, bool]]

class CoreFont:
    __slots__ = ("i", "type", "name", "sp", "ss", "up", "ut", "cw", "fontkey", "emphasis")
    i: int
    type: str
    name: str
    up: int
    ut: int
    sp: int
    ss: int
    cw: int
    fontkey: str
    emphasis: TextEmphasis
    def __init__(self, fpdf, fontkey: str, style: int) -> None: ...
    def get_text_width(self, text: str, font_size_pt: int, _: Unused) -> float: ...
    def encode_text(self, text: str) -> str: ...

class TTFFont:
    __slots__ = (
        "i",
        "type",
        "name",
        "desc",
        "glyph_ids",
        "hbfont",
        "sp",
        "ss",
        "up",
        "ut",
        "cw",
        "ttffile",
        "fontkey",
        "emphasis",
        "scale",
        "subset",
        "cmap",
        "ttfont",
        "missing_glyphs",
    )
    i: int
    type: str
    ttffile: Incomplete
    fontkey: str
    ttfont: _TTFont
    scale: float
    desc: PDFFontDescriptor
    cw: defaultdict[str, int]
    cmap: Incomplete
    glyph_ids: dict[Incomplete, Incomplete]
    missing_glyphs: list[Incomplete]
    name: str
    up: int
    ut: int
    sp: int
    ss: int
    emphasis: TextEmphasis
    subset: SubsetMap
    hbfont: HarfBuzzFont | None  # Not always defined.
    def __init__(self, fpdf, font_file_path, fontkey: str, style: int) -> None: ...
    def __deepcopy__(self, memo) -> Self: ...
    def close(self) -> None: ...
    def get_text_width(self, text: str, font_size_pt: int, text_shaping_params): ...
    def shaped_text_width(self, text: str, font_size_pt: int, text_shaping_params): ...
    def perform_harfbuzz_shaping(self, text: str, font_size_pt: int, text_shaping_params): ...
    def encode_text(self, text: str) -> str: ...
    def shape_text(self, text: str, font_size_pt: int, text_shaping_params): ...

class PDFFontDescriptor(PDFObject):
    type: Incomplete
    ascent: Incomplete
    descent: Incomplete
    cap_height: Incomplete
    flags: Incomplete
    font_b_box: Incomplete
    italic_angle: Incomplete
    stem_v: Incomplete
    missing_width: Incomplete
    font_name: Incomplete
    def __init__(self, ascent, descent, cap_height, flags, font_b_box, italic_angle, stem_v, missing_width) -> None: ...

@dataclass(order=True)
class Glyph:
    __slots__ = ("glyph_id", "unicode", "glyph_name", "glyph_width")
    glyph_id: int
    unicode: tuple[Incomplete, ...]
    glyph_name: str
    glyph_width: int
    def __hash__(self) -> int: ...

class SubsetMap:
    font: TTFFont
    def __init__(self, font: TTFFont) -> None: ...
    def __len__(self) -> int: ...
    def items(self) -> Generator[Incomplete]: ...
    def pick(self, unicode: int): ...
    def pick_glyph(self, glyph): ...
    def get_glyph(self, glyph=None, unicode=None, glyph_name=None, glyph_width=None) -> Glyph: ...
    def get_all_glyph_names(self): ...

CORE_FONTS: dict[str, str]
COURIER_FONT: dict[str, int]
CORE_FONTS_CHARWIDTHS: dict[str, dict[str, int]]
