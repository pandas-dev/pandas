import dataclasses
from pathlib import Path
import io
import os
from enum import Enum
from collections.abc import Generator

from typing import NamedTuple
from typing import Self

from .ft2font import CharacterCodeType, GlyphIndexType


class _dvistate(Enum):
    pre = ...
    outer = ...
    inpage = ...
    post_post = ...
    finale = ...

class Page(NamedTuple):
    text: list[Text]
    boxes: list[Box]
    height: int
    width: int
    descent: int

class Box(NamedTuple):
    x: int
    y: int
    height: int
    width: int

class Text(NamedTuple):
    x: int
    y: int
    font: DviFont
    glyph: CharacterCodeType
    width: int
    @property
    def font_path(self) -> Path: ...
    @property
    def font_size(self) -> float: ...
    @property
    def font_effects(self) -> dict[str, float]: ...
    @property
    def index(self) -> GlyphIndexType: ...  # type: ignore[override]
    @property
    def glyph_name_or_index(self) -> GlyphIndexType | str: ...

class Dvi:
    file: io.BufferedReader
    dpi: float | None
    fonts: dict[int, DviFont]
    state: _dvistate
    def __init__(self, filename: str | os.PathLike, dpi: float | None) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(self, etype, evalue, etrace) -> None: ...
    def __iter__(self) -> Generator[Page, None, None]: ...
    def close(self) -> None: ...

class DviFont:
    texname: bytes
    def __init__(
        self, scale: float, metrics: Tfm | TtfMetrics, texname: bytes, vf: Vf | None
    ) -> None: ...
    @classmethod
    def from_luatex(cls, scale: float, texname: bytes) -> DviFont: ...
    @classmethod
    def from_xetex(
        cls, scale: float, texname: bytes, subfont: int, effects: dict[str, float]
    ) -> DviFont: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def size(self) -> float: ...
    @property
    def widths(self) -> list[int]: ...
    @property
    def fname(self) -> str: ...
    @property
    def face_index(self) -> int: ...
    def resolve_path(self) -> Path: ...
    @property
    def subfont(self) -> int: ...
    @property
    def effects(self) -> dict[str, float]: ...

class Vf(Dvi):
    def __init__(self, filename: str | os.PathLike) -> None: ...
    def __getitem__(self, code: int) -> Page: ...

@dataclasses.dataclass(frozen=True, kw_only=True)
class TexMetrics:
    tex_width: int
    tex_height: int
    tex_depth: int
    # work around mypy not respecting kw_only=True in stub files
    __match_args__ = ()

class Tfm:
    checksum: int
    design_size: int
    def __init__(self, filename: str | os.PathLike) -> None: ...
    def get_metrics(self, idx: int) -> TexMetrics | None: ...
    @property
    def width(self) -> dict[int, int]: ...
    @property
    def height(self) -> dict[int, int]: ...
    @property
    def depth(self) -> dict[int, int]: ...

class TtfMetrics:
    def __init__(self, filename: str | os.PathLike) -> None: ...
    def get_metrics(self, idx: int) -> TexMetrics: ...

class PsFont(NamedTuple):
    texname: bytes
    psname: bytes
    effects: dict[str, float]
    encoding: None | bytes
    filename: str

class PsfontsMap:
    def __new__(cls, filename: str | os.PathLike) -> Self: ...
    def __getitem__(self, texname: bytes) -> PsFont: ...

def find_tex_file(filename: str | os.PathLike) -> str: ...
