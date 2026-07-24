import io
import json
import os
from _typeshed import Incomplete, SupportsRead
from collections.abc import Callable, Hashable
from typing import Any, Literal, final, overload
from typing_extensions import Self

import pandas as pd
from numpy.typing import ArrayLike
from pandas._typing import Axes, Dtype
from pyproj import CRS
from shapely.geometry.base import BaseGeometry

from ._decorator import doc
from .array import GeometryArray
from .base import GeoPandasBase, _BboxLike, _ClipMask, _ConvertibleToCRS, _ConvertibleToGeoSeries, _MaskLike
from .explore import _explore_geoseries
from .io._geoarrow import GeoArrowArray
from .plotting import plot_series

class GeoSeries(GeoPandasBase, pd.Series[BaseGeometry]):  # type: ignore[type-var,misc]  # pyright: ignore[reportInvalidTypeArguments]
    # Override the weird annotation of Series.__new__ in pandas-stubs
    def __new__(
        self,
        data: _ConvertibleToGeoSeries | None = None,
        index: Axes | None = None,
        crs: _ConvertibleToCRS | None = None,
        *,
        dtype: Dtype | None = None,
        name: Hashable = None,
        copy: bool | None = None,
        fastpath: bool = False,
    ) -> Self: ...
    def __init__(
        self,
        data: _ConvertibleToGeoSeries | None = None,
        index: Axes | None = None,
        crs: _ConvertibleToCRS | None = None,
        *,
        dtype: Dtype | None = None,
        name: Hashable = None,
        copy: bool | None = None,
        fastpath: bool = False,
    ) -> None: ...
    @final  # type: ignore[misc]
    def copy(self, deep: bool = True) -> Self: ...
    @property
    def values(self) -> GeometryArray: ...
    @property
    def geometry(self) -> Self: ...
    @property
    def x(self) -> pd.Series[float]: ...
    @property
    def y(self) -> pd.Series[float]: ...
    @property
    def z(self) -> pd.Series[float]: ...
    @property
    def m(self) -> pd.Series[float]: ...
    # Keep inline with GeoDataFrame.from_file and geopandas.io.file._read_file
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
    ) -> GeoSeries: ...
    @classmethod
    def from_wkb(
        cls,
        data: ArrayLike,  # array-like of bytes handled by shapely.from_wkb(data)
        index: Axes | None = None,
        crs: _ConvertibleToCRS | None = None,
        on_invalid: Literal["raise", "warn", "ignore", "fix"] = "raise",
        *,
        dtype: Dtype | None = None,
        name: Hashable = None,
        copy: bool | None = None,
        fastpath: bool = False,
    ) -> Self: ...
    @classmethod
    def from_wkt(
        cls,
        data: ArrayLike,  # array-like of str handled by shapely.from_wkt(data)
        index: Axes | None = None,
        crs: _ConvertibleToCRS | None = None,
        on_invalid: Literal["raise", "warn", "ignore", "fix"] = "raise",
        *,
        dtype: Dtype | None = None,
        name: Hashable = None,
        copy: bool | None = None,
        fastpath: bool = False,
    ) -> Self: ...
    @classmethod
    def from_xy(
        cls,
        # x, y, z: array-like of floats handled by np.asarray(..., dtype="float64")
        x: ArrayLike,
        y: ArrayLike,
        z: ArrayLike | None = None,
        index: Axes | None = None,
        crs: _ConvertibleToCRS | None = None,
        *,
        dtype: Dtype | None = None,
        name: Hashable = None,
        copy: bool | None = None,
        fastpath: bool = False,
    ) -> Self: ...
    @classmethod
    def from_arrow(
        cls,
        arr,
        *,
        # GeoSeries constructor kwargs
        index: Axes | None = None,
        crs: _ConvertibleToCRS | None = None,
        dtype: Dtype | None = None,
        name: Hashable = None,
        copy: bool | None = None,
        fastpath: bool = False,
    ) -> Self: ...
    @property
    def __geo_interface__(self) -> dict[str, Any]: ...  # values are arbitrary
    # Keep method to_file roughly in line with GeoDataFrame.to_file
    def to_file(
        self,
        filename: str | os.PathLike[str] | io.BytesIO,
        driver: str | None = None,
        index: bool | None = None,
        *,
        # kwargs from `_to_file` function
        schema: dict[str, Incomplete] | None = None,
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
    # *** `__getitem__`, `sort_index` and `take` are annotated with `-> Self` in pandas-stubs; no need to override them ***
    # *** `apply` annotation in pandas-stubs is compatible except for deprecated `convert_dtype` argument ***
    # def apply(self, func, convert_dtype: bool | None = None, args=(), **kwargs): ...
    def isna(self) -> pd.Series[bool]: ...
    def isnull(self) -> pd.Series[bool]: ...
    def notna(self) -> pd.Series[bool]: ...
    def notnull(self) -> pd.Series[bool]: ...
    # *** TODO: `fillna` annotation in pandas-stubs is NOT compatible; must `-> Self` ***
    # def fillna(self, value=None, method: FillnaOptions | None = None, inplace: bool = False, **kwargs): ...
    def __contains__(self, other: object) -> bool: ...  # type: ignore[misc]
    @doc(plot_series)
    def plot(self, *args, **kwargs): ...  # type: ignore[override]  # signature of `plot_series` copied in `@doc`
    @doc(_explore_geoseries)  # pyright: ignore[reportUnknownArgumentType]
    def explore(self, *args, **kwargs): ...  # signature of `_explore_geoseries` copied in `@doc`
    def explode(self, ignore_index: bool = False, index_parts: bool = False) -> GeoSeries: ...
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
    def to_crs(self, crs: _ConvertibleToCRS, epsg: int | None = None) -> GeoSeries: ...
    @overload
    def to_crs(self, crs: _ConvertibleToCRS | None = None, *, epsg: int) -> GeoSeries: ...
    @overload
    def to_crs(self, crs: _ConvertibleToCRS | None, epsg: int) -> GeoSeries: ...
    def estimate_utm_crs(self, datum_name: str = "WGS 84") -> CRS: ...
    def to_json(  # type: ignore[override]
        self,
        show_bbox: bool = True,
        drop_id: bool = False,
        to_wgs84: bool = False,
        *,
        # Keywords from json.dumps
        skipkeys: bool = False,
        ensure_ascii: bool = True,
        check_circular: bool = True,
        allow_nan: bool = True,
        cls: type[json.JSONEncoder] | None = None,
        indent: None | int | str = None,
        separators: tuple[str, str] | None = None,
        default: Callable[..., Any] | None = None,  # as typed in the json stdlib module
        sort_keys: bool = False,
        **kwds,
    ) -> str: ...
    @overload
    def to_wkb(self, hex: Literal[False] = False, **kwargs) -> pd.Series[bytes]: ...
    @overload
    def to_wkb(self, hex: Literal[True], **kwargs) -> pd.Series[str]: ...
    @overload
    def to_wkb(self, hex: bool = False, **kwargs) -> pd.Series[str] | pd.Series[bytes]: ...
    def to_wkt(self, **kwargs) -> pd.Series[str]: ...
    def to_arrow(
        self,
        geometry_encoding: Literal["WKB", "geoarrow", "wkb", "GeoArrow"] = "WKB",
        interleaved: bool | None = True,
        include_z: bool | None = None,
    ) -> GeoArrowArray: ...
    def clip(self, mask: _ClipMask, keep_geom_type: bool = False, sort: bool = False) -> GeoSeries: ...  # type: ignore[override]
