from collections.abc import Sequence
from typing import Any, Literal, SupportsIndex, overload

from ._enum import ParamEnum
from ._typing import ArrayLike, ArrayLikeSeq, GeoArray, OptGeoArrayLike, OptGeoArrayLikeSeq, OptGeoT
from .geometry import GeometryCollection, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon
from .geometry.base import BaseGeometry, BaseMultipartGeometry
from .lib import Geometry

__all__ = [
    "BufferCapStyle",
    "BufferJoinStyle",
    "boundary",
    "buffer",
    "build_area",
    "centroid",
    "clip_by_rect",
    "concave_hull",
    "constrained_delaunay_triangles",
    "convex_hull",
    "delaunay_triangles",
    "envelope",
    "extract_unique_points",
    "make_valid",
    "maximum_inscribed_circle",
    "minimum_bounding_circle",
    "minimum_clearance_line",
    "minimum_rotated_rectangle",
    "node",
    "normalize",
    "offset_curve",
    "orient_polygons",
    "oriented_envelope",
    "point_on_surface",
    "polygonize",
    "polygonize_full",
    "remove_repeated_points",
    "reverse",
    "segmentize",
    "simplify",
    "snap",
    "voronoi_polygons",
]

class BufferCapStyle(ParamEnum):
    round = 1
    flat = 2
    square = 3

class BufferJoinStyle(ParamEnum):
    round = 1
    mitre = 2
    bevel = 3

@overload
def boundary(geometry: Point | MultiPoint, **kwargs) -> GeometryCollection: ...
@overload
def boundary(geometry: LineString | MultiLineString, **kwargs) -> MultiPoint: ...
@overload
def boundary(geometry: Polygon | MultiPolygon, **kwargs) -> MultiLineString: ...
@overload
def boundary(geometry: GeometryCollection | None, **kwargs) -> None: ...
@overload
def boundary(geometry: Geometry, **kwargs) -> BaseMultipartGeometry | Any: ...
@overload
def boundary(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def buffer(
    geometry: Geometry,
    distance: float,
    quad_segs: int = 8,
    cap_style: BufferJoinStyle | Literal["round", "square", "flat"] = "round",
    join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = "round",
    mitre_limit: float = 5.0,
    single_sided: bool = False,
    **kwargs,
) -> Polygon: ...
@overload
def buffer(
    geometry: None,
    distance: float,
    quad_segs: int = 8,
    cap_style: BufferJoinStyle | Literal["round", "square", "flat"] = "round",
    join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = "round",
    mitre_limit: float = 5.0,
    single_sided: bool = False,
    **kwargs,
) -> None: ...
@overload
def buffer(
    geometry: Geometry | None,
    distance: float,
    quad_segs: int = 8,
    cap_style: BufferJoinStyle | Literal["round", "square", "flat"] = "round",
    join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = "round",
    mitre_limit: float = 5.0,
    single_sided: bool = False,
    **kwargs,
) -> Polygon | None: ...
@overload
def buffer(
    geometry: OptGeoArrayLike,
    distance: ArrayLikeSeq[float],
    quad_segs: int = 8,
    cap_style: BufferJoinStyle | Literal["round", "square", "flat"] = "round",
    join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = "round",
    mitre_limit: float = 5.0,
    single_sided: bool = False,
    **kwargs,
) -> GeoArray: ...
@overload
def buffer(
    geometry: OptGeoArrayLikeSeq,
    distance: ArrayLike[float],
    quad_segs: int = 8,
    cap_style: BufferJoinStyle | Literal["round", "square", "flat"] = "round",
    join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = "round",
    mitre_limit: float = 5.0,
    single_sided: bool = False,
    **kwargs,
) -> GeoArray: ...
@overload
def offset_curve(
    geometry: Geometry,
    distance: float,
    quad_segs: SupportsIndex = 8,
    join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = "round",
    mitre_limit: float = 5.0,
    **kwargs,
) -> LineString | MultiLineString: ...
@overload
def offset_curve(
    geometry: None,
    distance: float,
    quad_segs: SupportsIndex = 8,
    join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = "round",
    mitre_limit: float = 5.0,
    **kwargs,
) -> None: ...
@overload
def offset_curve(
    geometry: Geometry | None,
    distance: float,
    quad_segs: SupportsIndex = 8,
    join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = "round",
    mitre_limit: float = 5.0,
    **kwargs,
) -> LineString | MultiLineString | None: ...
@overload
def offset_curve(
    geometry: OptGeoArrayLike,
    distance: ArrayLikeSeq[float],
    quad_segs: SupportsIndex = 8,
    join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = "round",
    mitre_limit: float = 5.0,
    **kwargs,
) -> GeoArray: ...
@overload
def offset_curve(
    geometry: OptGeoArrayLikeSeq,
    distance: ArrayLike[float],
    quad_segs: SupportsIndex = 8,
    join_style: BufferJoinStyle | Literal["round", "mitre", "bevel"] = "round",
    mitre_limit: float = 5.0,
    **kwargs,
) -> GeoArray: ...
@overload
def centroid(geometry: Geometry, **kwargs) -> Point: ...
@overload
def centroid(geometry: None, **kwargs) -> None: ...
@overload
def centroid(geometry: Geometry | None, **kwargs) -> Point | None: ...
@overload
def centroid(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def clip_by_rect(geometry: Geometry, xmin: float, ymin: float, xmax: float, ymax: float, **kwargs) -> BaseGeometry: ...
@overload
def clip_by_rect(geometry: None, xmin: float, ymin: float, xmax: float, ymax: float, **kwargs) -> None: ...
@overload
def clip_by_rect(
    geometry: Geometry | None, xmin: float, ymin: float, xmax: float, ymax: float, **kwargs
) -> BaseGeometry | None: ...
@overload
def clip_by_rect(geometry: OptGeoArrayLikeSeq, xmin: float, ymin: float, xmax: float, ymax: float, **kwargs) -> GeoArray: ...
@overload
def concave_hull(geometry: Geometry, ratio: float = 0.0, allow_holes: bool = False, **kwargs) -> BaseGeometry: ...
@overload
def concave_hull(geometry: None, ratio: float = 0.0, allow_holes: bool = False, **kwargs) -> None: ...
@overload
def concave_hull(geometry: Geometry | None, ratio: float = 0.0, allow_holes: bool = False, **kwargs) -> BaseGeometry | None: ...
@overload
def concave_hull(geometry: OptGeoArrayLikeSeq, ratio: float = 0.0, allow_holes: bool = False, **kwargs) -> GeoArray: ...
@overload
def convex_hull(geometry: Geometry, **kwargs) -> BaseGeometry: ...
@overload
def convex_hull(geometry: None, **kwargs) -> None: ...
@overload
def convex_hull(geometry: Geometry | None, **kwargs) -> BaseGeometry | None: ...
@overload
def convex_hull(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def delaunay_triangles(
    geometry: Geometry, tolerance: float = 0.0, only_edges: Literal[False] = False, **kwargs
) -> GeometryCollection: ...
@overload
def delaunay_triangles(geometry: Geometry, tolerance: float, only_edges: Literal[True], **kwargs) -> MultiLineString: ...
@overload
def delaunay_triangles(geometry: Geometry, tolerance: float = 0.0, *, only_edges: Literal[True], **kwargs) -> MultiLineString: ...
@overload
def delaunay_triangles(
    geometry: Geometry, tolerance: float = 0.0, only_edges: bool = False, **kwargs
) -> GeometryCollection | MultiLineString: ...
@overload
def delaunay_triangles(geometry: None, tolerance: float = 0.0, only_edges: bool = False, **kwargs) -> None: ...
@overload
def delaunay_triangles(
    geometry: Geometry | None, tolerance: float = 0.0, only_edges: bool = False, **kwargs
) -> GeometryCollection | MultiLineString | None: ...
@overload
def delaunay_triangles(
    geometry: OptGeoArrayLike, tolerance: ArrayLike[float], only_edges: ArrayLikeSeq[bool], **kwargs
) -> GeoArray: ...
@overload
def delaunay_triangles(
    geometry: OptGeoArrayLike, tolerance: ArrayLike[float] = 0.0, *, only_edges: ArrayLikeSeq[bool], **kwargs
) -> GeoArray: ...
@overload
def delaunay_triangles(
    geometry: OptGeoArrayLike, tolerance: ArrayLikeSeq[float], only_edges: ArrayLike[bool] = False, **kwargs
) -> GeoArray: ...
@overload
def delaunay_triangles(
    geometry: OptGeoArrayLikeSeq, tolerance: ArrayLike[float] = 0.0, only_edges: ArrayLike[bool] = False, **kwargs
) -> GeoArray: ...
@overload
def constrained_delaunay_triangles(geometry: Geometry, **kwargs) -> GeometryCollection: ...
@overload
def constrained_delaunay_triangles(geometry: None, **kwargs) -> None: ...
@overload
def constrained_delaunay_triangles(geometry: Geometry | None, **kwargs) -> GeometryCollection | None: ...
@overload
def constrained_delaunay_triangles(geometry: OptGeoArrayLikeSeq | OptGeoArrayLike, **kwargs) -> GeoArray: ...
@overload
def envelope(geometry: Point, **kwargs) -> Point: ...
@overload
def envelope(geometry: Geometry, **kwargs) -> BaseGeometry: ...
@overload
def envelope(geometry: None, **kwargs) -> None: ...
@overload
def envelope(geometry: Geometry | None, **kwargs) -> BaseGeometry | None: ...
@overload
def envelope(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def extract_unique_points(geometry: Geometry, **kwargs) -> MultiPoint: ...
@overload
def extract_unique_points(geometry: None, **kwargs) -> None: ...
@overload
def extract_unique_points(geometry: Geometry | None, **kwargs) -> MultiPoint | None: ...
@overload
def extract_unique_points(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def build_area(geometry: Geometry, **kwargs) -> BaseGeometry: ...
@overload
def build_area(geometry: None, **kwargs) -> None: ...
@overload
def build_area(geometry: Geometry | None, **kwargs) -> BaseGeometry | None: ...
@overload
def build_area(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...

# make_valid with `method="linework"` only accepts `keep_collapsed=True`
@overload
def make_valid(
    geometry: Geometry, *, method: Literal["linework"] = "linework", keep_collapsed: Literal[True] = True, **kwargs
) -> BaseGeometry: ...
@overload
def make_valid(
    geometry: None, *, method: Literal["linework"] = "linework", keep_collapsed: Literal[True] = True, **kwargs
) -> None: ...
@overload
def make_valid(
    geometry: Geometry | None, *, method: Literal["linework"] = "linework", keep_collapsed: Literal[True] = True, **kwargs
) -> BaseGeometry | None: ...
@overload
def make_valid(
    geometry: OptGeoArrayLikeSeq, *, method: Literal["linework"] = "linework", keep_collapsed: Literal[True] = True, **kwargs
) -> GeoArray: ...
@overload
def make_valid(geometry: Geometry, *, method: Literal["structure"], keep_collapsed: bool = True, **kwargs) -> BaseGeometry: ...
@overload
def make_valid(geometry: None, *, method: Literal["structure"], keep_collapsed: bool = True, **kwargs) -> None: ...
@overload
def make_valid(
    geometry: Geometry | None, *, method: Literal["structure"], keep_collapsed: bool = True, **kwargs
) -> BaseGeometry | None: ...
@overload
def make_valid(
    geometry: OptGeoArrayLikeSeq, *, method: Literal["structure"], keep_collapsed: bool = True, **kwargs
) -> GeoArray: ...
@overload
def minimum_clearance_line(geometry: Point, **kwargs) -> Point: ...
@overload
def minimum_clearance_line(geometry: LineString | Polygon | BaseMultipartGeometry, **kwargs) -> Polygon: ...
@overload
def minimum_clearance_line(geometry: Geometry, **kwargs) -> Polygon | Point: ...
@overload
def minimum_clearance_line(geometry: None, **kwargs) -> None: ...
@overload
def minimum_clearance_line(geometry: Geometry | None, **kwargs) -> Polygon | Point | None: ...
@overload
def minimum_clearance_line(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def normalize(geometry: OptGeoT, **kwargs) -> OptGeoT: ...
@overload
def normalize(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def point_on_surface(geometry: Geometry, **kwargs) -> Point: ...
@overload
def point_on_surface(geometry: None, **kwargs) -> None: ...
@overload
def point_on_surface(geometry: Geometry | None, **kwargs) -> Point | None: ...
@overload
def point_on_surface(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def node(geometry: Geometry, **kwargs) -> MultiLineString: ...
@overload
def node(geometry: None, **kwargs) -> None: ...
@overload
def node(geometry: Geometry | None, **kwargs) -> MultiLineString | None: ...
@overload
def node(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def polygonize(geometries: Sequence[Geometry | None], **kwargs) -> GeometryCollection: ...
@overload
def polygonize(geometries: Sequence[Sequence[Geometry | None]], **kwargs) -> GeoArray: ...
@overload
def polygonize(geometries: OptGeoArrayLikeSeq, **kwargs) -> GeometryCollection | GeoArray: ...
@overload
def polygonize_full(
    geometries: Sequence[Geometry | None], **kwargs
) -> tuple[GeometryCollection, GeometryCollection, GeometryCollection, GeometryCollection]: ...
@overload
def polygonize_full(
    geometries: Sequence[Sequence[Geometry | None]], **kwargs
) -> tuple[GeoArray, GeoArray, GeoArray, GeoArray]: ...
@overload
def polygonize_full(
    geometries: OptGeoArrayLikeSeq, **kwargs
) -> (
    tuple[GeometryCollection, GeometryCollection, GeometryCollection, GeometryCollection]
    | tuple[GeoArray, GeoArray, GeoArray, GeoArray]
): ...
@overload
def remove_repeated_points(geometry: OptGeoT, tolerance: float = 0.0, **kwargs) -> OptGeoT: ...
@overload
def remove_repeated_points(geometry: OptGeoArrayLikeSeq, tolerance: float = 0.0, **kwargs) -> GeoArray: ...
@overload
def reverse(geometry: OptGeoT, **kwargs) -> OptGeoT: ...
@overload
def reverse(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def segmentize(geometry: OptGeoT, max_segment_length: float, **kwargs) -> OptGeoT: ...
@overload
def segmentize(geometry: OptGeoArrayLike, max_segment_length: ArrayLikeSeq[float], **kwargs) -> GeoArray: ...
@overload
def segmentize(geometry: OptGeoArrayLikeSeq, max_segment_length: ArrayLike[float], **kwargs) -> GeoArray: ...
@overload
def simplify(geometry: OptGeoT, tolerance: float, preserve_topology: bool = True, **kwargs) -> OptGeoT: ...
@overload
def simplify(geometry: OptGeoArrayLike, tolerance: ArrayLikeSeq[float], preserve_topology: bool = True, **kwargs) -> GeoArray: ...
@overload
def simplify(geometry: OptGeoArrayLikeSeq, tolerance: ArrayLike[float], preserve_topology: bool = True, **kwargs) -> GeoArray: ...
@overload
def snap(geometry: OptGeoT, reference: Geometry, tolerance: float, **kwargs) -> OptGeoT: ...
@overload
def snap(geometry: Geometry | None, reference: None, tolerance: float, **kwargs) -> None: ...
@overload
def snap(geometry: OptGeoArrayLikeSeq, reference: OptGeoArrayLike, tolerance: ArrayLike[float], **kwargs) -> GeoArray: ...
@overload
def snap(geometry: OptGeoArrayLike, reference: OptGeoArrayLikeSeq, tolerance: ArrayLike[float], **kwargs) -> GeoArray: ...
@overload
def snap(geometry: OptGeoArrayLike, reference: OptGeoArrayLike, tolerance: ArrayLikeSeq[float], **kwargs) -> GeoArray: ...
@overload
def voronoi_polygons(
    geometry: Geometry,
    tolerance: float = 0.0,
    extend_to: Geometry | None = None,
    only_edges: Literal[False] = False,
    ordered: bool = False,
    **kwargs,
) -> GeometryCollection[Polygon]: ...
@overload
def voronoi_polygons(
    geometry: Geometry, tolerance: float, extend_to: Geometry | None, only_edges: Literal[True], ordered: bool = False, **kwargs
) -> LineString | MultiLineString: ...
@overload
def voronoi_polygons(
    geometry: Geometry,
    tolerance: float = 0.0,
    extend_to: Geometry | None = None,
    *,
    only_edges: Literal[True],
    ordered: bool = False,
    **kwargs,
) -> LineString | MultiLineString: ...
@overload
def voronoi_polygons(
    geometry: Geometry,
    tolerance: float = 0.0,
    extend_to: Geometry | None = None,
    only_edges: bool = False,
    ordered: bool = False,
    **kwargs,
) -> GeometryCollection[Polygon] | LineString | MultiLineString: ...
@overload
def voronoi_polygons(
    geometry: None,
    tolerance: float = 0.0,
    extend_to: Geometry | None = None,
    only_edges: bool = False,
    ordered: bool = False,
    **kwargs,
) -> None: ...
@overload
def voronoi_polygons(
    geometry: Geometry | None,
    tolerance: float = 0.0,
    extend_to: Geometry | None = None,
    only_edges: bool = False,
    ordered: bool = False,
    **kwargs,
) -> GeometryCollection[Polygon] | LineString | MultiLineString | None: ...
@overload  # `geometry` as sequence-like
def voronoi_polygons(
    geometry: OptGeoArrayLikeSeq,
    tolerance: ArrayLike[float] = 0.0,
    extend_to: OptGeoArrayLike = None,
    only_edges: ArrayLike[bool] = False,
    ordered: ArrayLike[bool] = False,
    **kwargs,
) -> GeoArray: ...
@overload  # `tolerance` as sequence-like
def voronoi_polygons(
    geometry: OptGeoArrayLike,
    tolerance: ArrayLikeSeq[float],
    extend_to: OptGeoArrayLike = None,
    only_edges: ArrayLike[bool] = False,
    ordered: ArrayLike[bool] = False,
    **kwargs,
) -> GeoArray: ...
@overload
def voronoi_polygons(  # `extend_to` as positional sequence-like
    geometry: OptGeoArrayLike,
    tolerance: ArrayLike[float],
    extend_to: OptGeoArrayLikeSeq,
    only_edges: ArrayLike[bool] = False,
    ordered: ArrayLike[bool] = False,
    **kwargs,
) -> GeoArray: ...
@overload  # `extend_to` as keyword sequence-like
def voronoi_polygons(
    geometry: OptGeoArrayLike,
    tolerance: ArrayLike[float] = 0.0,
    *,
    extend_to: OptGeoArrayLikeSeq,
    only_edges: ArrayLike[bool] = False,
    ordered: ArrayLike[bool] = False,
    **kwargs,
) -> GeoArray: ...
@overload
def voronoi_polygons(  # `only_edges` as positional sequence-like
    geometry: OptGeoArrayLike,
    tolerance: ArrayLike[float],
    extend_to: OptGeoArrayLike,
    only_edges: ArrayLikeSeq[bool],
    ordered: ArrayLike[bool] = False,
    **kwargs,
) -> GeoArray: ...
@overload
def voronoi_polygons(  # `only_edges` as keyword sequence-like
    geometry: OptGeoArrayLike,
    tolerance: ArrayLike[float] = 0.0,
    extend_to: OptGeoArrayLike = None,
    *,
    only_edges: ArrayLikeSeq[bool],
    ordered: ArrayLike[bool] = False,
    **kwargs,
) -> GeoArray: ...
@overload  # `ordered` as positional sequence-like
def voronoi_polygons(
    geometry: OptGeoArrayLike,
    tolerance: ArrayLike[float],
    extend_to: OptGeoArrayLike,
    only_edges: ArrayLike[bool],
    ordered: ArrayLikeSeq[bool],
    **kwargs,
) -> GeoArray: ...
@overload  # `ordered` as keyword sequence-like
def voronoi_polygons(
    geometry: OptGeoArrayLike,
    tolerance: ArrayLike[float] = 0.0,
    extend_to: OptGeoArrayLike = None,
    *,
    only_edges: ArrayLike[bool] = False,
    ordered: ArrayLikeSeq[bool],
    **kwargs,
) -> GeoArray: ...
@overload
def oriented_envelope(geometry: Point, **kwargs) -> Point: ...
@overload
def oriented_envelope(geometry: Geometry, **kwargs) -> BaseGeometry: ...
@overload
def oriented_envelope(geometry: None, **kwargs) -> None: ...
@overload
def oriented_envelope(geometry: Geometry | None, **kwargs) -> BaseGeometry | None: ...
@overload
def oriented_envelope(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...

minimum_rotated_rectangle = oriented_envelope

@overload
def minimum_bounding_circle(geometry: Point, **kwargs) -> Point: ...
@overload
def minimum_bounding_circle(geometry: LineString | Polygon | BaseMultipartGeometry, **kwargs) -> Polygon: ...
@overload
def minimum_bounding_circle(geometry: Geometry, **kwargs) -> Polygon | Point: ...
@overload
def minimum_bounding_circle(geometry: None, **kwargs) -> None: ...
@overload
def minimum_bounding_circle(geometry: Geometry | None, **kwargs) -> Polygon | Point | None: ...
@overload
def minimum_bounding_circle(geometry: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def maximum_inscribed_circle(geometry: Polygon | MultiPolygon, tolerance: float | None = None, **kwargs) -> LineString: ...
@overload
def maximum_inscribed_circle(geometry: None, tolerance: float | None = None, **kwargs) -> None: ...
@overload
def maximum_inscribed_circle(
    geometry: Polygon | MultiPolygon | None, tolerance: float | None = None, **kwargs
) -> LineString | None: ...
@overload
def maximum_inscribed_circle(geometry: OptGeoArrayLikeSeq, tolerance: ArrayLike[float] | None = None, **kwargs) -> GeoArray: ...
@overload
def maximum_inscribed_circle(geometry: OptGeoArrayLike, tolerance: ArrayLikeSeq[float], **kwargs) -> GeoArray: ...
@overload
def orient_polygons(geometry: OptGeoT, *, exterior_cw: bool = False, **kwargs) -> OptGeoT: ...
@overload
def orient_polygons(geometry: OptGeoArrayLikeSeq, *, exterior_cw: bool = False, **kwargs) -> GeoArray: ...
