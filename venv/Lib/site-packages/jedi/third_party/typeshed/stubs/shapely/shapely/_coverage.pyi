from typing import overload

import numpy as np

from ._typing import ArrayLike, GeoArray, GeoArrayLikeSeq, OptGeoArrayLike
from .geometry import Polygon

__all__ = ["coverage_invalid_edges", "coverage_is_valid", "coverage_simplify"]

def coverage_is_valid(geometry: OptGeoArrayLike, gap_width: float = 0.0, **kwargs) -> np.bool_: ...
def coverage_invalid_edges(geometry: OptGeoArrayLike, gap_width: float = 0.0, **kwargs) -> GeoArray: ...
@overload
def coverage_simplify(geometry: Polygon, tolerance: ArrayLike[float], *, simplify_boundary: bool = True) -> Polygon: ...
@overload
def coverage_simplify(geometry: GeoArrayLikeSeq, tolerance: ArrayLike[float], *, simplify_boundary: bool = True) -> GeoArray: ...
