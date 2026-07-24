from collections.abc import Iterable
from typing import Any, Final, Literal, overload

import numpy as np
from numpy.typing import ArrayLike, NDArray
from shapely import Geometry

from .array import _Array1D, _Array2D

PREDICATES: Final[set[str | None]]

class SpatialIndex:
    geometries: NDArray[np.object_]
    def __init__(self, geometry: NDArray[np.object_]) -> None: ...
    @property
    def valid_query_predicates(self) -> set[str | None]: ...
    @overload
    def query(
        self,
        geometry: Geometry | ArrayLike,
        predicate: str | None = None,
        sort: bool = False,
        distance: float | ArrayLike | None = None,
        output_format: Literal["indices"] = "indices",
    ) -> NDArray[np.int64]: ...
    @overload
    def query(
        self,
        geometry: Geometry | ArrayLike,
        predicate: str | None = None,
        sort: bool = False,
        distance: float | ArrayLike | None = None,
        *,
        output_format: Literal["dense"],
    ) -> NDArray[np.bool_]: ...
    @overload
    def query(
        self,
        geometry: Geometry | ArrayLike,
        predicate: str | None = None,
        sort: bool = False,
        distance: float | ArrayLike | None = None,
        *,
        output_format: Literal["sparse"],
    ) -> Any: ...  # returns scipy coo_array but we don't depend on scipy
    @overload
    def nearest(
        self,
        geometry,
        return_all: bool = True,
        max_distance: float | None = None,
        return_distance: Literal[False] = False,
        exclusive: bool = False,
    ) -> _Array2D[np.int64]: ...
    @overload
    def nearest(
        self,
        geometry,
        return_all: bool = True,
        max_distance: float | None = None,
        *,
        return_distance: Literal[True],
        exclusive: bool = False,
    ) -> tuple[_Array2D[np.int64], _Array1D[np.float64]]: ...
    @overload
    def nearest(
        self,
        geometry,
        return_all: bool = True,
        max_distance: float | None = None,
        return_distance: bool = False,
        exclusive: bool = False,
    ) -> _Array2D[np.int64] | tuple[_Array2D[np.int64], _Array1D[np.float64]]: ...
    def intersection(self, coordinates: Iterable[float]) -> _Array1D[np.int64]: ...
    @property
    def size(self) -> int: ...
    @property
    def is_empty(self) -> bool: ...
    def __len__(self) -> int: ...
