from _typeshed import SupportsRead
from pathlib import Path

from PIL import Image

from ..._types import Ink
from ..styledpil import StyledPilImage

class QRColorMask:
    back_color: Ink
    has_transparency: bool
    paint_color: Ink
    # image is not actually used by any of the initialize implementations in this project.
    def initialize(self, styledPilImage: StyledPilImage, image: Image.Image) -> None: ...
    def apply_mask(self, image: Image.Image, use_cache: bool = False) -> None: ...
    def get_fg_pixel(self, image: Image.Image, x: int, y: int) -> Ink: ...
    def get_bg_pixel(self, image: Image.Image, x: int, y: int) -> Ink: ...
    def interp_num(self, n1: int, n2: int, norm: float) -> int: ...
    def interp_color(self, col1: Ink, col2: Ink, norm: float) -> Ink: ...
    def extrap_num(self, n1: int, n2: int, interped_num: int) -> float | None: ...
    def extrap_color(self, col1: Ink, col2: Ink, interped_color: Ink) -> float | None: ...

class SolidFillColorMask(QRColorMask):
    front_color: Ink
    def __init__(self, back_color: Ink = (255, 255, 255), front_color: Ink = (0, 0, 0)) -> None: ...
    def apply_mask(self, image: Image.Image) -> None: ...  # type: ignore[override]

class RadialGradiantColorMask(QRColorMask):
    center_color: Ink
    edge_color: Ink
    def __init__(
        self, back_color: Ink = (255, 255, 255), center_color: Ink = (0, 0, 0), edge_color: Ink = (0, 0, 255)
    ) -> None: ...

class SquareGradiantColorMask(QRColorMask):
    center_color: Ink
    edge_color: Ink
    def __init__(
        self, back_color: Ink = (255, 255, 255), center_color: Ink = (0, 0, 0), edge_color: Ink = (0, 0, 255)
    ) -> None: ...

class HorizontalGradiantColorMask(QRColorMask):
    left_color: Ink
    right_color: Ink
    def __init__(
        self, back_color: Ink = (255, 255, 255), left_color: Ink = (0, 0, 0), right_color: Ink = (0, 0, 255)
    ) -> None: ...

class VerticalGradiantColorMask(QRColorMask):
    top_color: Ink
    bottom_color: Ink
    def __init__(
        self, back_color: Ink = (255, 255, 255), top_color: Ink = (0, 0, 0), bottom_color: Ink = (0, 0, 255)
    ) -> None: ...

class ImageColorMask(QRColorMask):
    color_img: Ink
    def __init__(
        self,
        back_color: Ink = (255, 255, 255),
        color_mask_path: str | bytes | Path | SupportsRead[bytes] | None = None,
        color_mask_image: Image.Image | None = None,
    ) -> None: ...
    paint_color: Ink
