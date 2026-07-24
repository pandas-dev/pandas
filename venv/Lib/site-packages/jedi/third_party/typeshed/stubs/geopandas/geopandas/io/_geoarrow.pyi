from _typeshed import Incomplete
from collections.abc import Mapping
from typing import (
    # pyarrow types returned as Any to avoid depending on pyarrow (40 MB) in stubs
    Any as _PAArray,
    Any as _PAField,
    Any as _PATable,
    Literal,
    Protocol,
    type_check_only,
)
from typing_extensions import CapsuleType, TypeAlias

import numpy as np
from numpy.typing import NDArray

from ..array import GeometryArray
from ..geodataframe import GeoDataFrame

# Literal for language server completions and str because runtime normalizes to lowercase
_GeomEncoding: TypeAlias = Literal["WKB", "geoarrow"] | str  # noqa: Y051

@type_check_only
class _PyarrowTableLike(Protocol):
    def __arrow_c_stream__(self, requested_schema=None) -> CapsuleType: ...

@type_check_only
class _PyarrowFieldLike(Protocol):
    def __arrow_c_schema__(self) -> CapsuleType: ...

@type_check_only
class _PyarrowArrayLike(Protocol):
    def __arrow_c_array__(self) -> tuple[CapsuleType, CapsuleType]: ...

GEOARROW_ENCODINGS: list[str]

class ArrowTable:
    def __init__(self, pa_table: _PyarrowTableLike) -> None: ...
    def __arrow_c_stream__(self, requested_schema=None) -> CapsuleType: ...

class GeoArrowArray:
    def __init__(self, pa_field: _PyarrowFieldLike, pa_array: _PyarrowArrayLike) -> None: ...
    def __arrow_c_array__(self, requested_schema=None) -> tuple[CapsuleType, CapsuleType]: ...

def geopandas_to_arrow(
    df: GeoDataFrame,
    index: bool | None = None,
    geometry_encoding: _GeomEncoding = "WKB",
    interleaved: bool = True,
    include_z: bool | None = None,
) -> tuple[_PATable, dict[str, str]]: ...
def construct_wkb_array(
    shapely_arr: NDArray[np.object_], *, field_name: str = "geometry", crs: str | None = None
) -> tuple[_PAField, _PAArray]: ...
def construct_geometry_array(
    shapely_arr: NDArray[np.object_],
    include_z: bool | None = None,
    *,
    field_name: str = "geometry",
    crs: str | None = None,
    interleaved: bool = True,
) -> tuple[_PAField, _PAArray]: ...
def arrow_to_geopandas(
    table, geometry: str | None = None, to_pandas_kwargs: Mapping[str, Incomplete] | None = None
) -> GeoDataFrame: ...
def arrow_to_geometry_array(arr) -> GeometryArray: ...
def construct_shapely_array(arr: _PAArray, extension_name: str) -> NDArray[np.object_]: ...
