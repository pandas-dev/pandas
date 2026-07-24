from collections.abc import Collection
from typing_extensions import Self

from .base import BaseMultipartGeometry
from .multilinestring import MultiLineString
from .polygon import Polygon, _PolygonHolesLike, _PolygonShellLike

__all__ = ["MultiPolygon"]

class MultiPolygon(BaseMultipartGeometry[Polygon]):
    __slots__: list[str] = []
    def __new__(
        self,
        polygons: (
            BaseMultipartGeometry
            | Collection[Polygon | tuple[_PolygonShellLike] | tuple[_PolygonShellLike, _PolygonHolesLike] | None]
            | None
        ) = None,
    ) -> Self: ...
    def svg(self, scale_factor: float = 1.0, fill_color: str | None = None, opacity: float | None = None) -> str: ...  # type: ignore[override]
    # more precise base overrides
    @property
    def boundary(self) -> MultiLineString: ...
