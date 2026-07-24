from _typeshed import Incomplete
from dataclasses import dataclass
from typing import Any, Literal
from typing_extensions import TypeAlias

from fpdf.enums import Align

from .image_parsing import _ImageFilter

_AlignLiteral: TypeAlias = Literal[
    "",
    "CENTER",
    "X_CENTER",
    "LEFT",
    "RIGHT",
    "JUSTIFY",
    "center",
    "x_center",
    "left",
    "right",
    "justify",
    "C",
    "X",
    "L",
    "R",
    "J",
    "c",
    "x",
    "l",
    "r",
    "j",
]
_TextAlign: TypeAlias = Align | _AlignLiteral  # noqa: Y047

class ImageInfo(dict[str, Any]):
    @property
    def width(self) -> int: ...
    @property
    def height(self) -> int: ...
    @property
    def rendered_width(self) -> int: ...
    @property
    def rendered_height(self) -> int: ...
    def scale_inside_box(self, x: float, y: float, w: float, h: float) -> tuple[float, float, float, float]: ...

class RasterImageInfo(ImageInfo):
    def size_in_document_units(self, w: float, h: float, scale=1) -> tuple[float, float]: ...

class VectorImageInfo(ImageInfo): ...

@dataclass
class ImageCache:
    images: dict[str, dict[Incomplete, Incomplete]] = ...
    icc_profiles: dict[bytes, int] = ...
    image_filter: _ImageFilter = "AUTO"

    def reset_usages(self) -> None: ...
