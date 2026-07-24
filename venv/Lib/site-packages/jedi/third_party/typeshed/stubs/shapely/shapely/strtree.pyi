from typing import Any, Literal, SupportsIndex, overload
from typing_extensions import TypeAlias

import numpy as np
from numpy.typing import NDArray

from ._enum import ParamEnum
from ._typing import ArrayLike, GeoArray, GeoArrayLikeSeq, OptGeoArrayLike
from .lib import Geometry

__all__ = ["STRtree"]

_BinaryPredicate: TypeAlias = Literal[
    "intersects", "within", "contains", "overlaps", "crosses", "touches", "covers", "covered_by", "contains_properly"
]

class BinaryPredicate(ParamEnum):
    intersects = 1
    within = 2
    contains = 3
    overlaps = 4
    crosses = 5
    touches = 6
    covers = 7
    covered_by = 8
    contains_properly = 9

class STRtree:
    def __init__(self, geoms: GeoArrayLikeSeq, node_capacity: SupportsIndex = 10) -> None: ...
    def __len__(self) -> int: ...
    @property
    def geometries(self) -> GeoArray: ...
    @overload
    def query(
        self, geometry: OptGeoArrayLike, predicate: Literal["dwithin"], distance: ArrayLike[float]
    ) -> NDArray[np.int64]: ...
    @overload
    def query(
        self, geometry: OptGeoArrayLike, predicate: _BinaryPredicate | None = None, distance: object = None
    ) -> NDArray[np.int64]: ...
    # nearest may return `None` if the tree is empty, use the "Any trick"
    @overload
    def nearest(self, geometry: Geometry) -> np.int64 | Any: ...
    @overload
    def nearest(self, geometry: GeoArrayLikeSeq) -> NDArray[np.int64] | Any: ...
    @overload  # return_distance=False
    def query_nearest(
        self,
        geometry: OptGeoArrayLike,
        max_distance: float | None = None,
        return_distance: Literal[False] = False,
        exclusive: bool = False,
        all_matches: bool = True,
    ) -> NDArray[np.int64]: ...
    @overload  # return_distance=True keyword
    def query_nearest(
        self,
        geometry: OptGeoArrayLike,
        max_distance: float | None = None,
        *,
        return_distance: Literal[True],
        exclusive: bool = False,
        all_matches: bool = True,
    ) -> tuple[NDArray[np.int64], NDArray[np.float64]]: ...
    @overload  # return_distance=True positional
    def query_nearest(
        self,
        geometry: OptGeoArrayLike,
        max_distance: float | None,
        return_distance: Literal[True],
        exclusive: bool = False,
        all_matches: bool = True,
    ) -> tuple[NDArray[np.int64], NDArray[np.float64]]: ...
    @overload  # return_distance=bool fallback
    def query_nearest(
        self,
        geometry: OptGeoArrayLike,
        max_distance: float | None = None,
        return_distance: bool = False,
        exclusive: bool = False,
        all_matches: bool = True,
    ) -> NDArray[np.int64] | tuple[NDArray[np.int64], NDArray[np.float64]]: ...
