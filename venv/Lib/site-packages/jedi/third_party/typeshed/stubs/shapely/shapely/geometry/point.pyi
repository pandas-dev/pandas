from collections.abc import Iterable
from typing import overload
from typing_extensions import Self, TypeAlias

from .._typing import ArrayLikeSeq
from .base import BaseGeometry
from .collection import GeometryCollection

__all__ = ["Point"]

_PointLike: TypeAlias = Point | Iterable[float] | ArrayLikeSeq[float]

class Point(BaseGeometry):
    __slots__: list[str] = []
    @overload  # no args: empty point
    def __new__(self) -> Self: ...
    @overload  # one arg: (x, y[, z]) tuple or a Point instance
    def __new__(self, coords: _PointLike, /) -> Self: ...
    @overload  # two args: (x, y) tuple
    def __new__(self, x: float, y: float, /) -> Self: ...
    @overload  # three args: (x, y, z) tuple
    def __new__(self, x: float, y: float, z: float, /) -> Self: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    @property
    def z(self) -> float: ...
    @property
    def m(self) -> float: ...
    def svg(self, scale_factor: float = 1.0, fill_color: str | None = None, opacity: float | None = None) -> str: ...  # type: ignore[override]
    # more precise base overrides
    @property
    def boundary(self) -> GeometryCollection: ...  # empty geometry collection
    @property
    def convex_hull(self) -> Point: ...
    @property
    def envelope(self) -> Point: ...
    @property
    def oriented_envelope(self) -> Point: ...
    @property
    def minimum_rotated_rectangle(self) -> Point: ...
