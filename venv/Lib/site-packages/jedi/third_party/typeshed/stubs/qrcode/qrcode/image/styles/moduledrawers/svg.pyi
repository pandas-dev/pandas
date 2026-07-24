import abc
from decimal import Decimal
from typing import Any, NamedTuple
from xml.etree.ElementTree import Element, QName

from ...._types import Box
from ...svg import SvgFragmentImage, SvgPathImage
from .base import QRModuleDrawer

ANTIALIASING_FACTOR: int

class Coords(NamedTuple):
    x0: Decimal
    y0: Decimal
    x1: Decimal
    y1: Decimal
    xh: Decimal
    yh: Decimal

class BaseSvgQRModuleDrawer(QRModuleDrawer, metaclass=abc.ABCMeta):
    img: SvgFragmentImage
    size_ratio: Decimal
    # kwargs are used to allow for subclasses with additional keyword arguments
    def __init__(self, *, size_ratio: Decimal = ..., **kwargs: Any) -> None: ...
    box_delta: float
    box_size: Decimal
    box_half: Decimal
    def coords(self, box: Box) -> Coords: ...

class SvgQRModuleDrawer(BaseSvgQRModuleDrawer, metaclass=abc.ABCMeta):
    tag: str
    tag_qname: QName
    def drawrect(self, box: Box, is_active: bool) -> None: ...  # type: ignore[override]
    @abc.abstractmethod
    def el(self, box: Box) -> Element: ...

class SvgSquareDrawer(SvgQRModuleDrawer):
    unit_size: str
    def el(self, box: Box) -> Element: ...

class SvgCircleDrawer(SvgQRModuleDrawer):
    tag: str
    radius: str
    def el(self, box: Box) -> Element: ...

class SvgPathQRModuleDrawer(BaseSvgQRModuleDrawer, metaclass=abc.ABCMeta):
    img: SvgPathImage
    def drawrect(self, box: Box, is_active: bool) -> None: ...  # type: ignore[override]
    @abc.abstractmethod
    def subpath(self, box: Box) -> str: ...

class SvgPathSquareDrawer(SvgPathQRModuleDrawer):
    def subpath(self, box: Box) -> str: ...

class SvgPathCircleDrawer(SvgPathQRModuleDrawer):
    def subpath(self, box: Box) -> str: ...
