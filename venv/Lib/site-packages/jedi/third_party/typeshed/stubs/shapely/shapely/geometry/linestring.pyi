from collections.abc import Iterable
from typing import Literal, SupportsFloat, SupportsIndex
from typing_extensions import Self, TypeAlias

from .._typing import ArrayLikeSeq
from ..constructive import BufferJoinStyle
from .base import BaseGeometry
from .multilinestring import MultiLineString
from .multipoint import MultiPoint
from .point import Point
from .polygon import Polygon

__all__ = ["LineString"]

_ConvertibleToLineString: TypeAlias = LineString | ArrayLikeSeq[float] | Iterable[Point | Iterable[SupportsFloat]]

class LineString(BaseGeometry):
    __slots__: list[str] = []
    def __new__(self, coordinates: _ConvertibleToLineString | None = None) -> Self: ...
    def svg(self, scale_factor: float = 1.0, stroke_color: str | None = None, opacity: float | None = None) -> str: ...  # type: ignore[override]
    def offset_curve(
        self,
        distance: float,
        quad_segs: SupportsIndex = 16,
        join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = ...,
        mitre_limit: float = 5.0,
    ) -> LineString | MultiLineString: ...
    def parallel_offset(  # to be deprecated
        self,
        distance: float,
        side: str = "right",
        resolution: SupportsIndex = 16,
        join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = ...,
        mitre_limit: float = 5.0,
    ) -> LineString | MultiLineString: ...
    # more precise base overrides
    @property
    def boundary(self) -> MultiPoint: ...
    @property
    def convex_hull(self) -> LineString: ...
    @property
    def envelope(self) -> Polygon: ...
    @property
    def oriented_envelope(self) -> LineString: ...
    @property
    def minimum_rotated_rectangle(self) -> LineString: ...
