from collections.abc import Callable, Iterable
from typing import Any, Literal, overload

from ._typing import GeoT, OptGeoArrayLike, SupportsGeoInterface
from .algorithms.polylabel import polylabel as polylabel
from .geometry import GeometryCollection, LineString, MultiLineString, Point, Polygon
from .geometry.base import BaseGeometry, BaseMultipartGeometry, GeometrySequence
from .geometry.linestring import _ConvertibleToLineString
from .lib import Geometry

__all__ = [
    "clip_by_rect",
    "linemerge",
    "nearest_points",
    "operator",
    "orient",
    "polygonize",
    "polygonize_full",
    "shared_paths",
    "snap",
    "split",
    "substring",
    "transform",
    "triangulate",
    "unary_union",
    "validate",
    "voronoi_diagram",
]

class CollectionOperator:
    @overload
    def shapeup(self, ob: GeoT) -> GeoT: ...  # type: ignore[overload-overlap]
    @overload
    def shapeup(self, ob: dict[str, Any] | SupportsGeoInterface) -> BaseGeometry: ...  # type: ignore[overload-overlap]
    @overload
    def shapeup(self, ob: _ConvertibleToLineString) -> LineString: ...
    def polygonize(
        self, lines: OptGeoArrayLike | Iterable[_ConvertibleToLineString | None]
    ) -> GeometrySequence[GeometryCollection[Polygon]]: ...
    def polygonize_full(
        self, lines: OptGeoArrayLike | Iterable[_ConvertibleToLineString | None]
    ) -> tuple[
        GeometryCollection[Polygon], GeometryCollection[LineString], GeometryCollection[LineString], GeometryCollection[Polygon]
    ]: ...
    def linemerge(
        self, lines: MultiLineString | BaseMultipartGeometry | Iterable[_ConvertibleToLineString], directed: bool = False
    ) -> LineString | MultiLineString: ...
    def unary_union(self, geoms: OptGeoArrayLike) -> BaseGeometry: ...

operator: CollectionOperator
polygonize = operator.polygonize
polygonize_full = operator.polygonize_full
linemerge = operator.linemerge
unary_union = operator.unary_union

# This is also an alias to operator method but we want to mark it as deprecated
@overload  # edges false
def triangulate(geom: Geometry, tolerance: float = 0.0, edges: Literal[False] = False) -> list[Polygon]: ...
@overload  # edges true (keyword)
def triangulate(geom: Geometry, tolerance: float = 0.0, *, edges: Literal[True]) -> list[LineString]: ...
@overload  # edges true (positional)
def triangulate(geom: Geometry, tolerance: float, edges: Literal[True]) -> list[LineString]: ...
@overload  # fallback
def triangulate(geom: Geometry, tolerance: float = 0.0, edges: bool = False) -> list[Polygon] | list[LineString]: ...
@overload
def voronoi_diagram(
    geom: Geometry, envelope: Geometry | None = None, tolerance: float = 0.0, edges: Literal[False] = False
) -> GeometryCollection[Polygon]: ...
@overload
def voronoi_diagram(
    geom: Geometry, envelope: Geometry | None, tolerance: float, edges: Literal[True]
) -> GeometryCollection[LineString | MultiLineString]: ...
@overload
def voronoi_diagram(
    geom: Geometry, envelope: Geometry | None = None, tolerance: float = 0.0, *, edges: Literal[True]
) -> GeometryCollection[LineString | MultiLineString]: ...
@overload
def voronoi_diagram(
    geom: Geometry, envelope: Geometry | None = None, tolerance: float = 0.0, edges: bool = False
) -> GeometryCollection[Polygon | LineString | MultiLineString]: ...
@overload
def validate(geom: None) -> None: ...
@overload
def validate(geom: Geometry) -> str: ...
@overload
def validate(geom: Geometry | None) -> str | None: ...
def transform(func: Callable[[float, float, float | None], tuple[float, ...]], geom: GeoT) -> GeoT: ...
def nearest_points(g1: Geometry, g2: Geometry) -> tuple[Point, Point]: ...
def snap(g1: GeoT, g2: Geometry, tolerance: float) -> GeoT: ...
def shared_paths(g1: LineString, g2: LineString) -> GeometryCollection[MultiLineString]: ...

class SplitOp:
    @staticmethod
    def split(geom: Geometry, splitter: Geometry) -> GeometryCollection: ...

split = SplitOp.split

def substring(geom: LineString, start_dist: float, end_dist: float, normalized: bool = False) -> Point | LineString: ...
def clip_by_rect(geom: Geometry, xmin: float, ymin: float, xmax: float, ymax: float) -> BaseGeometry: ...
def orient(geom: GeoT, sign: float = 1.0) -> GeoT: ...
