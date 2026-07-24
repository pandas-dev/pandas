import abc
from typing import Literal

from PIL import Image, ImageDraw

from ...._types import Box
from ....main import ActiveWithNeighbors
from ...styledpil import StyledPilImage
from .base import QRModuleDrawer

ANTIALIASING_FACTOR: int

class StyledPilQRModuleDrawer(QRModuleDrawer, metaclass=abc.ABCMeta):
    img: StyledPilImage

class SquareModuleDrawer(StyledPilQRModuleDrawer):
    imgDraw: ImageDraw.ImageDraw
    def drawrect(self, box: Box, is_active: bool) -> None: ...  # type: ignore[override]

class GappedSquareModuleDrawer(StyledPilQRModuleDrawer):
    size_ratio: float
    def __init__(self, size_ratio: float = 0.8) -> None: ...
    imgDraw: ImageDraw.ImageDraw
    delta: float
    def drawrect(self, box: Box, is_active: bool) -> None: ...  # type: ignore[override]

class CircleModuleDrawer(StyledPilQRModuleDrawer):
    circle: Image.Image
    def drawrect(self, box: Box, is_active: bool) -> None: ...  # type: ignore[override]

class RoundedModuleDrawer(StyledPilQRModuleDrawer):
    needs_neighbors: Literal[True]
    radius_ratio: float
    def __init__(self, radius_ratio: float = 1) -> None: ...
    corner_width: int
    SQUARE: Image.Image
    NW_ROUND: Image.Image
    SW_ROUND: Image.Image
    SE_ROUND: Image.Image
    NE_ROUND: Image.Image
    def setup_corners(self) -> None: ...
    def drawrect(self, box: Box, is_active: ActiveWithNeighbors) -> None: ...  # type: ignore[override]

class VerticalBarsDrawer(StyledPilQRModuleDrawer):
    needs_neighbors: Literal[True]
    horizontal_shrink: float
    def __init__(self, horizontal_shrink: float = 0.8) -> None: ...
    half_height: int
    delta: int
    SQUARE: Image.Image
    ROUND_TOP: Image.Image
    ROUND_BOTTOM: Image.Image
    def setup_edges(self) -> None: ...
    def drawrect(self, box: Box, is_active: ActiveWithNeighbors) -> None: ...  # type: ignore[override]

class HorizontalBarsDrawer(StyledPilQRModuleDrawer):
    needs_neighbors: Literal[True]
    vertical_shrink: float
    def __init__(self, vertical_shrink: float = 0.8) -> None: ...
    half_width: int
    delta: int
    SQUARE: Image.Image
    ROUND_LEFT: Image.Image
    ROUND_RIGHT: Image.Image
    def setup_edges(self) -> None: ...
    def drawrect(self, box: Box, is_active: ActiveWithNeighbors) -> None: ...  # type: ignore[override]
