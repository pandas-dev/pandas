from collections.abc import Collection
from typing_extensions import Self

from .base import BaseMultipartGeometry
from .linestring import LineString, _ConvertibleToLineString
from .multipoint import MultiPoint

__all__ = ["MultiLineString"]

class MultiLineString(BaseMultipartGeometry[LineString]):
    __slots__: list[str] = []
    def __new__(self, lines: BaseMultipartGeometry | Collection[_ConvertibleToLineString] | None = None) -> Self: ...
    def svg(self, scale_factor: float = 1.0, stroke_color: str | None = None, opacity: float | None = None) -> str: ...  # type: ignore[override]
    # more precise base overrides
    @property
    def boundary(self) -> MultiPoint: ...
