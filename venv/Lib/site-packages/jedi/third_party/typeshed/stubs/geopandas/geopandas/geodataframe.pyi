import io
import os
from _typeshed import Incomplete, SupportsGetItem, SupportsLenAndGetItem, SupportsRead, SupportsWrite
from collections.abc import Callable, Container, Hashable, Iterable, Iterator, Mapping, Sequence
from json import JSONEncoder
from typing import Any, Literal, overload
from typing_extensions import Self

import pandas as pd
from numpy.typing import ArrayLike
from pandas._typing import AggFuncTypeFrame, Axes, Axis, Dtype, GroupByObject, IndexLabel, Scalar
from pyproj import CRS

from ._decorator import doc
from .base import (
    GeoPandasBase,
    _BboxLike,
    _ClipMask,
    _ConvertibleToCRS,
    _ConvertibleToDataFrame,
    _GeomCol,
    _GeomSeq,
    _MaskLike,
    _SupportsGeoInterface,
)
from .explore import _explore
from .geoseries import GeoSeries
from .io._geoarrow import ArrowTable, _GeomEncoding
from .io.sql import _SQLConnection
from .plotting import GeoplotAccessor

crs_mismatch_error: str

class GeoDataFrame(GeoPandasBase, pd.DataFrame):  # type: ignore[misc]
    # Override the weird annotation of DataFrame.__new__ in pandas-stubs
    @overload
    def __new__(
        cls,
        data: _ConvertibleToDataFrame | None = None,
        index: Axes | None = None,
        columns: Axes | None = None,
        dtype: Dtype | None = None,
        copy: bool | None = None,
        *,
        geometry: _GeomCol | None = None,
        crs: _ConvertibleToCRS | None = None,
    ) -> Self: ...
    @overload
    def __new__(
        cls,
        data: Scalar,
        index: Axes,
        columns: Axes,
        dtype: Dtype | None = None,
        copy: bool | None = None,
        *,
        geometry: _GeomCol | None = None,
        crs: _ConvertibleToCRS | None = None,
    ) -> Self: ...
    def __init__(
        self,
        data: _ConvertibleToDataFrame | None = None,
        index: Axes | None = None,
        columns: Axes | None = None,
        dtype: Dtype | None = None,
        copy: bool | None = None,
        *,
        geometry: _GeomCol | None = None,
        crs: _ConvertibleToCRS | None = None,
    ) -> None: ...
    def __setattr__(self, attr: str, val: Any) -> None: ...  # type: ignore[misc]  # Can set arbitrary objects
    @property
    def geometry(self) -> GeoSeries: ...
    @geometry.setter
    def geometry(self, col: _GeomSeq) -> None: ...
    @overload
    def set_geometry(
        self, col: _GeomCol, drop: bool | None = None, inplace: Literal[False] = False, crs: _ConvertibleToCRS | None = None
    ) -> Self: ...
    @overload
    def set_geometry(
        self, col: _GeomCol, drop: bool | None = None, *, inplace: Literal[True], crs: _ConvertibleToCRS | None = None
    ) -> None: ...
    @overload
    def set_geometry(
        self, col: _GeomCol, drop: bool | None, inplace: Literal[True], crs: _ConvertibleToCRS | None = None
    ) -> None: ...
    @overload
    def rename_geometry(self, col: Hashable, inplace: Literal[False] = False) -> Self: ...
    @overload
    def rename_geometry(self, col: Hashable, inplace: Literal[True]) -> None: ...
    @property
    def active_geometry_name(self) -> str | None: ...
    @property
    def crs(self) -> CRS | None: ...
    @crs.setter
    def crs(self, value: _ConvertibleToCRS | None) -> None: ...
    @classmethod
    def from_dict(  # type: ignore[override]
        # Mapping[Any, Any] because of invariance keys and arbitrary values
        cls,
        data: Mapping[Any, Any],
        geometry: _GeomCol | None = None,
        crs: _ConvertibleToCRS | None = None,
        **kwargs,
    ) -> Self: ...
    # Keep inline with GeoSeries.from_file and geopandas.io.file._read_file
    @classmethod
    def from_file(
        cls,
        filename: str | os.PathLike[str] | SupportsRead[Incomplete],
        *,
        bbox: _BboxLike | None = None,
        mask: _MaskLike | None = None,
        rows: int | slice | None = None,
        engine: Literal["fiona", "pyogrio"] | None = None,
        ignore_geometry: Literal[False] = False,
        layer: int | str | None = None,
        encoding: str | None = None,
        **kwargs,  # engine dependent
    ) -> Self: ...
    @classmethod
    def from_features(
        cls,
        features: (
            _SupportsGeoInterface
            | Mapping[str, _SupportsGeoInterface | SupportsGetItem[str, Incomplete]]
            | Iterable[_SupportsGeoInterface | SupportsGetItem[str, Incomplete]]
        ),
        crs: _ConvertibleToCRS | None = None,
        columns: Axes | None = None,
    ) -> Self: ...
    @overload
    @classmethod
    def from_postgis(
        cls,
        sql: str,
        con: _SQLConnection,
        geom_col: str = "geom",
        crs: _ConvertibleToCRS | None = None,
        index_col: str | list[str] | None = None,
        coerce_float: bool = True,
        parse_dates: Container[str | Mapping[str, Incomplete]] | Mapping[str, str | Mapping[str, Incomplete]] | None = None,
        params: SupportsLenAndGetItem[Scalar] | Mapping[str, Scalar] | None = None,
        *,
        chunksize: int,
    ) -> Iterator[GeoDataFrame]: ...
    @overload
    @classmethod
    def from_postgis(
        cls,
        sql: str,
        con: _SQLConnection,
        geom_col: str = "geom",
        crs: _ConvertibleToCRS | None = None,
        index_col: str | list[str] | None = None,
        coerce_float: bool = True,
        parse_dates: Container[str | Mapping[str, Incomplete]] | Mapping[str, str | Mapping[str, Incomplete]] | None = None,
        params: SupportsLenAndGetItem[Scalar] | Mapping[str, Scalar] | None = None,
        chunksize: None = None,
    ) -> GeoDataFrame: ...
    @classmethod
    def from_arrow(
        cls, table, geometry: str | None = None, to_pandas_kwargs: Mapping[str, Incomplete] | None = None
    ) -> GeoDataFrame: ...  # TODO: `table: pyarrow.Table | table-like`
    def to_json(  # type: ignore[override]
        self,
        na: str = "null",
        show_bbox: bool = False,
        drop_id: bool = False,
        to_wgs84: bool = False,
        *,
        # json.dumps kwargs
        skipkeys: bool = False,
        ensure_ascii: bool = True,
        check_circular: bool = True,
        allow_nan: bool = True,
        cls: type[JSONEncoder] | None = None,
        indent: int | str | None = None,
        separators: tuple[str, str] | None = None,
        default: Callable[..., Any] | None = None,  # as typed in the json stdlib module
        sort_keys: bool = False,
        **kwargs,
    ) -> str: ...
    @property
    def __geo_interface__(self) -> dict[str, Any]: ...  # values are arbitrary
    def iterfeatures(
        self, na: str = "null", show_bbox: bool = False, drop_id: bool = False
    ) -> Iterator[dict[str, Incomplete]]: ...
    def to_geo_dict(self, na: str = "null", show_bbox: bool = False, drop_id: bool = False) -> dict[str, Incomplete]: ...
    def to_wkb(
        self,
        hex: bool = False,
        *,
        # shapely kwargs
        output_dimension: int = ...,
        byte_order: int = ...,
        include_srid: bool = ...,
        flavor: Literal["iso", "extended"] = ...,
        **kwargs,
    ) -> pd.DataFrame: ...
    def to_wkt(
        self,
        *,
        # shapely kwargs
        rounding_precision: int = ...,
        trim: bool = ...,
        output_dimension: int = ...,
        old_3d: bool = ...,
        **kwargs,
    ) -> pd.DataFrame: ...
    def to_arrow(
        self,
        *,
        index: bool | None = None,
        geometry_encoding: _GeomEncoding = "WKB",
        interleaved: bool = True,
        include_z: bool | None = None,
    ) -> ArrowTable: ...
    def to_parquet(  # type: ignore[override]
        self,
        path: str | os.PathLike[str] | SupportsWrite[Incomplete],
        index: bool | None = None,
        compression: Literal["snappy", "gzip", "brotli", "lz4", "zstd"] | None = "snappy",
        geometry_encoding: _GeomEncoding = "WKB",
        write_covering_bbox: bool = False,
        schema_version: str | None = None,
        *,
        engine: Literal["auto", "pyarrow"] = "auto",  # Only these engines are supported, unlike pandas
        **kwargs,
    ) -> None: ...
    def to_feather(
        self,
        path: str | os.PathLike[str] | SupportsWrite[Incomplete],
        index: bool | None = None,
        compression: Literal["zstd", "lz4", "uncompressed"] | None = None,
        schema_version: str | None = None,
        **kwargs,
    ) -> None: ...
    # Keep method to_file roughly in line with GeoSeries.to_file
    def to_file(
        self,
        filename: str | os.PathLike[str] | io.BytesIO,
        driver: str | None = None,
        schema: dict[str, Incomplete] | None = None,
        index: bool | None = None,
        *,
        # kwargs from `_to_file` function
        mode: Literal["w", "a"] = "w",
        crs: _ConvertibleToCRS | None = None,
        engine: Literal["fiona", "pyogrio"] | None = None,
        metadata: dict[str, str] | None = None,
        # kwargs extracted from engines
        layer: int | str | None = None,
        encoding: str | None = None,
        overwrite: bool | None = ...,
        **kwargs,  # engine and driver dependent
    ) -> None: ...
    @overload
    def set_crs(
        self, crs: _ConvertibleToCRS, epsg: int | None = None, inplace: bool = False, allow_override: bool = False
    ) -> Self: ...
    @overload
    def set_crs(
        self, crs: _ConvertibleToCRS | None = None, *, epsg: int, inplace: bool = False, allow_override: bool = False
    ) -> Self: ...
    @overload
    def set_crs(self, crs: _ConvertibleToCRS | None, epsg: int, inplace: bool = False, allow_override: bool = False) -> Self: ...
    @overload
    def to_crs(self, crs: _ConvertibleToCRS, epsg: int | None = None, inplace: Literal[False] = False) -> Self: ...
    @overload
    def to_crs(self, crs: _ConvertibleToCRS | None = None, *, epsg: int, inplace: Literal[False] = False) -> Self: ...
    @overload
    def to_crs(self, crs: _ConvertibleToCRS | None, epsg: int, inplace: Literal[False] = False) -> Self: ...
    @overload
    def to_crs(self, crs: _ConvertibleToCRS, epsg: int | None = None, *, inplace: Literal[True]) -> None: ...
    @overload
    def to_crs(self, crs: _ConvertibleToCRS, epsg: int | None, inplace: Literal[True]) -> None: ...
    @overload
    def to_crs(self, crs: _ConvertibleToCRS | None = None, *, epsg: int, inplace: Literal[True]) -> None: ...
    @overload
    def to_crs(self, crs: _ConvertibleToCRS | None, epsg: int, inplace: Literal[True]) -> None: ...
    def estimate_utm_crs(self, datum_name: str = "WGS 84") -> CRS: ...
    # def __getitem__(self, key): ...
    def __delitem__(self, key) -> None: ...  # type: ignore[misc]
    # def __setitem__(self, key, value) -> None: ...
    def copy(self, deep: bool = True) -> Self: ...  # type: ignore[misc]
    # def merge(self, *args, **kwargs) -> GeoDataFrame | pd.DataFrame: ...
    def apply(  # type: ignore[override]
        self,
        func: Callable[..., Incomplete],
        axis: Axis = 0,
        raw: bool = False,
        result_type: Literal["expand", "reduce", "broadcast"] | None = None,
        args: tuple[Any, ...] = (),  # type inexpressible in the typing system
        *,
        by_row: Literal[False, "compat"] = "compat",
        engine: Literal["python", "numba"] = "python",
        engine_kwargs: dict[str, bool] | None = None,
        **kwargs,
    ) -> pd.DataFrame | pd.Series[Incomplete]: ...
    def __finalize__(self, other, method: str | None = None, **kwargs) -> Self: ...  # type: ignore[misc]
    def dissolve(
        self,
        by: GroupByObject | None = None,
        aggfunc: AggFuncTypeFrame = "first",
        as_index: bool = True,
        level: IndexLabel | None = None,
        sort: bool = True,
        observed: bool = False,
        dropna: bool = True,
        method: Literal["coverage", "unary", "disjoint_subset"] = "unary",
        grid_size: float | None = None,
        **kwargs,
    ) -> GeoDataFrame: ...
    def explode(self, column: IndexLabel | None = None, ignore_index: bool = False, index_parts: bool = False) -> Self: ...
    def to_postgis(
        self,
        name: str,
        con: _SQLConnection,
        schema: str | None = None,
        if_exists: Literal["fail", "replace", "append"] = "fail",
        index: bool = False,
        index_label: IndexLabel | None = None,
        chunksize: int | None = None,
        dtype: dict[Any, Incomplete] | None = None,  # columns can be of "any" type
    ) -> None: ...
    @property
    def plot(self) -> GeoplotAccessor: ...
    @doc(_explore)  # pyright: ignore[reportUnknownArgumentType]
    def explore(self, *args, **kwargs): ...  # signature of `_explore` copied in `@doc`
    def sjoin(
        self,
        df: GeoDataFrame,
        how: Literal["left", "right", "inner"] = "inner",
        predicate: str = "intersects",
        lsuffix: str = "left",
        rsuffix: str = "right",
        *,
        # **kwargs passed to geopandas.sjoin
        distance: float | ArrayLike | None = None,
        on_attribute: str | Sequence[str] | None = None,
        **kwargs,
    ) -> GeoDataFrame: ...
    def sjoin_nearest(
        self,
        right: GeoDataFrame,
        how: Literal["left", "right", "inner"] = "inner",
        max_distance: float | None = None,
        lsuffix: str = "left",
        rsuffix: str = "right",
        distance_col: str | None = None,
        exclusive: bool = False,
    ) -> GeoDataFrame: ...
    def clip(self, mask: _ClipMask, keep_geom_type: bool = False, sort: bool = False) -> GeoDataFrame: ...  # type: ignore[override]
    def overlay(
        self, right: GeoDataFrame, how: str = "intersection", keep_geom_type: bool | None = None, make_valid: bool = True
    ) -> GeoDataFrame: ...
