from typing import overload

import numpy as np
from numpy.typing import NDArray

from ._typing import ArrayLike, ArrayLikeSeq, GeoArray, OptGeoArrayLike, OptGeoArrayLikeSeq
from .geometry import GeometryCollection, LineString, MultiLineString, Point
from .lib import Geometry

__all__ = ["line_interpolate_point", "line_locate_point", "line_merge", "shared_paths", "shortest_line"]

@overload
def line_interpolate_point(line: None, distance: float, normalized: bool = False, **kwargs) -> None: ...
@overload
def line_interpolate_point(
    line: None, distance: ArrayLikeSeq[float], normalized: bool = False, **kwargs
) -> NDArray[np.object_]: ...  # Array of None
@overload
def line_interpolate_point(
    line: LineString | MultiLineString | GeometryCollection, distance: float, normalized: bool = False, **kwargs
) -> Point: ...
@overload
def line_interpolate_point(
    line: LineString | MultiLineString | GeometryCollection, distance: ArrayLikeSeq[float], normalized: bool = False, **kwargs
) -> GeoArray: ...
@overload
def line_interpolate_point(
    line: OptGeoArrayLikeSeq, distance: ArrayLike[float], normalized: bool = False, **kwargs
) -> GeoArray: ...
@overload
def line_locate_point(
    line: LineString | MultiLineString | GeometryCollection | None, other: Point | None, normalized: bool = False, **kwargs
) -> np.float64: ...
@overload
def line_locate_point(
    line: LineString | MultiLineString | GeometryCollection | None, other: OptGeoArrayLikeSeq, normalized: bool = False, **kwargs
) -> NDArray[np.float64]: ...
@overload
def line_locate_point(
    line: OptGeoArrayLikeSeq, other: OptGeoArrayLike, normalized: bool = False, **kwargs
) -> NDArray[np.float64]: ...
@overload
def line_merge(line: None, directed: bool = False, **kwargs) -> None: ...
@overload
def line_merge(line: Geometry, directed: bool = False, **kwargs) -> LineString | MultiLineString | GeometryCollection: ...
@overload
def line_merge(line: OptGeoArrayLikeSeq, directed: bool = False, **kwargs) -> GeoArray: ...
@overload
def shared_paths(a: LineString | MultiLineString | None, b: None, **kwargs) -> None: ...
@overload
def shared_paths(a: None, b: LineString | MultiLineString | None, **kwargs) -> None: ...
@overload
def shared_paths(
    a: LineString | MultiLineString, b: LineString | MultiLineString, **kwargs
) -> GeometryCollection[MultiLineString]: ...
@overload
def shared_paths(a: LineString | MultiLineString | None, b: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def shared_paths(a: OptGeoArrayLikeSeq, b: OptGeoArrayLike, **kwargs) -> GeoArray: ...
@overload
def shortest_line(a: Geometry | None, b: None, **kwargs) -> None: ...
@overload
def shortest_line(a: None, b: Geometry | None, **kwargs) -> None: ...
@overload
def shortest_line(a: Geometry, b: Geometry, **kwargs) -> LineString: ...
@overload
def shortest_line(a: Geometry | None, b: OptGeoArrayLikeSeq, **kwargs) -> GeoArray: ...
@overload
def shortest_line(a: OptGeoArrayLikeSeq, b: OptGeoArrayLike, **kwargs) -> GeoArray: ...
