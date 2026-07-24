from collections.abc import Iterable
from dataclasses import dataclass
from io import BytesIO
from logging import Logger
from types import TracebackType
from typing import Any, Final, Literal
from typing_extensions import TypeAlias

from PIL import Image

from .image_datastructures import ImageCache, ImageInfo, VectorImageInfo
from .svg import SVGObject

_ImageFilter: TypeAlias = Literal["AUTO", "FlateDecode", "DCTDecode", "JPXDecode", "LZWDecode"]

RESAMPLE: Image.Resampling

@dataclass
class ImageSettings:
    compression_level: int = -1

LOGGER: Logger
SUPPORTED_IMAGE_FILTERS: tuple[_ImageFilter, ...]
SETTINGS: ImageSettings

TIFFBitRevTable: list[int]

LZW_CLEAR_TABLE_MARKER: Final = 256
LZW_EOD_MARKER: Final = 257
LZW_INITIAL_BITS_PER_CODE: Final = 9
LZW_MAX_BITS_PER_CODE: Final = 12

def preload_image(
    image_cache: ImageCache, name: str | BytesIO | Image.Image, dims: tuple[float, float] | None = None
) -> tuple[str, BytesIO | Image.Image | None, ImageInfo]: ...
def load_image(filename): ...
def is_iccp_valid(iccp, filename) -> bool: ...
def get_svg_info(filename: str, img: BytesIO, image_cache: ImageCache) -> tuple[str, SVGObject, VectorImageInfo]: ...

# Returned dict could be typed as a TypedDict.
def get_img_info(
    filename, img: BytesIO | Image.Image | None = None, image_filter: _ImageFilter = "AUTO", dims=None
) -> dict[str, Any]: ...

class temp_attr:
    obj: Any
    field: str
    value: Any
    exists: bool  # defined after __enter__ is called
    def __init__(self, obj: Any, field: str, value: Any) -> None: ...
    def __enter__(self) -> None: ...
    def __exit__(
        self, exctype: type[BaseException] | None, excinst: BaseException | None, exctb: TracebackType | None
    ) -> None: ...

def ccitt_payload_location_from_pil(img: Image.Image) -> tuple[int, int]: ...
def transcode_monochrome(img: Image.Image): ...
def pack_codes_into_bytes(codes: Iterable[int]) -> bytes: ...
def clear_table() -> tuple[dict[bytes, int], int, int, int]: ...
