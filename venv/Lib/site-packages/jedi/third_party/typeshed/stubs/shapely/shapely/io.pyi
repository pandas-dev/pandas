from _typeshed import Incomplete
from typing import Literal, overload
from typing_extensions import TypeAlias

import numpy as np
from numpy.typing import NDArray

from ._enum import ParamEnum
from ._ragged_array import from_ragged_array as from_ragged_array, to_ragged_array as to_ragged_array
from ._typing import ArrayLikeSeq, GeoArray, OptGeoArrayLikeSeq
from .geometry.base import BaseGeometry
from .lib import Geometry

__all__ = ["from_geojson", "from_ragged_array", "from_wkb", "from_wkt", "to_geojson", "to_ragged_array", "to_wkb", "to_wkt"]

_OutputDimension: TypeAlias = Literal[2, 3, 4]

# Mypy and stubtest aren't happy with the following definition and
# raise is a reserved keyword, so we cannot use the class syntax of enums
# DecodingErrorOptions = ParamEnum("DecodingErrorOptions", {"ignore": 0, "warn": 1, "raise": 2, "fix": 3})
DecodingErrorOptions: Incomplete

class WKBFlavorOptions(ParamEnum):
    extended = 1
    iso = 2

@overload
def to_wkt(
    geometry: None,
    rounding_precision: int = 6,
    trim: bool = True,
    output_dimension: _OutputDimension | None = None,
    old_3d: bool = False,
    **kwargs,
) -> None: ...
@overload
def to_wkt(
    geometry: Geometry,
    rounding_precision: int = 6,
    trim: bool = True,
    output_dimension: _OutputDimension | None = None,
    old_3d: bool = False,
    **kwargs,
) -> str: ...
@overload
def to_wkt(
    geometry: OptGeoArrayLikeSeq,
    rounding_precision: int = 6,
    trim: bool = True,
    output_dimension: _OutputDimension | None = None,
    old_3d: bool = False,
    **kwargs,
) -> NDArray[np.str_]: ...
@overload
def to_wkb(
    geometry: None,
    hex: bool = False,
    output_dimension: _OutputDimension | None = None,
    byte_order: int = -1,
    include_srid: bool = False,
    flavor: Literal["iso", "extended"] = "extended",
    **kwargs,
) -> None: ...
@overload
def to_wkb(
    geometry: Geometry,
    hex: Literal[False] = False,
    output_dimension: _OutputDimension | None = None,
    byte_order: int = -1,
    include_srid: bool = False,
    flavor: Literal["iso", "extended"] = "extended",
    **kwargs,
) -> bytes: ...
@overload
def to_wkb(
    geometry: Geometry,
    hex: Literal[True],
    output_dimension: _OutputDimension | None = None,
    byte_order: int = -1,
    include_srid: bool = False,
    flavor: Literal["iso", "extended"] = "extended",
    **kwargs,
) -> str: ...
@overload
def to_wkb(
    geometry: Geometry,
    hex: bool,
    output_dimension: _OutputDimension | None = None,
    byte_order: int = -1,
    include_srid: bool = False,
    flavor: Literal["iso", "extended"] = "extended",
    **kwargs,
) -> bytes | str: ...
@overload
def to_wkb(
    geometry: OptGeoArrayLikeSeq,
    hex: Literal[False] = False,
    output_dimension: _OutputDimension | None = None,
    byte_order: int = -1,
    include_srid: bool = False,
    flavor: Literal["iso", "extended"] = "extended",
    **kwargs,
) -> NDArray[np.bytes_]: ...
@overload
def to_wkb(
    geometry: OptGeoArrayLikeSeq,
    hex: Literal[True],
    output_dimension: _OutputDimension | None = None,
    byte_order: int = -1,
    include_srid: bool = False,
    flavor: Literal["iso", "extended"] = "extended",
    **kwargs,
) -> NDArray[np.str_]: ...
@overload
def to_wkb(
    geometry: OptGeoArrayLikeSeq,
    hex: bool,
    output_dimension: _OutputDimension | None = None,
    byte_order: int = -1,
    include_srid: bool = False,
    flavor: Literal["iso", "extended"] = "extended",
    **kwargs,
) -> NDArray[np.bytes_] | NDArray[np.str_]: ...
@overload
def to_geojson(geometry: None, indent: int | None = None, **kwargs) -> None: ...
@overload
def to_geojson(geometry: Geometry, indent: int | None = None, **kwargs) -> str: ...
@overload
def to_geojson(geometry: OptGeoArrayLikeSeq, indent: int | None = None, **kwargs) -> NDArray[np.str_]: ...
@overload
def from_wkt(geometry: None, on_invalid: Literal["raise", "warn", "ignore", "fix"] = "raise", **kwargs) -> None: ...
@overload
def from_wkt(geometry: str, on_invalid: Literal["raise", "warn", "ignore", "fix"] = "raise", **kwargs) -> BaseGeometry: ...  # type: ignore[overload-overlap]
@overload
def from_wkt(
    geometry: ArrayLikeSeq[str | None], on_invalid: Literal["raise", "warn", "ignore", "fix"] = "raise", **kwargs
) -> GeoArray: ...
@overload
def from_wkb(geometry: None, on_invalid: Literal["raise", "warn", "ignore", "fix"] = "raise", **kwargs) -> None: ...
@overload
def from_wkb(
    geometry: str | bytes, on_invalid: Literal["raise", "warn", "ignore", "fix"] = "raise", **kwargs
) -> BaseGeometry: ...  # type: ignore[overload-overlap]
@overload
def from_wkb(
    geometry: ArrayLikeSeq[str | bytes | None], on_invalid: Literal["raise", "warn", "ignore", "fix"] = "raise", **kwargs
) -> GeoArray: ...
@overload
def from_geojson(geometry: None, on_invalid: Literal["raise", "warn", "ignore", "fix"] = "raise", **kwargs) -> None: ...
@overload
def from_geojson(
    geometry: str | bytes, on_invalid: Literal["raise", "warn", "ignore", "fix"] = "raise", **kwargs
) -> BaseGeometry: ...  # type: ignore[overload-overlap]
@overload
def from_geojson(
    geometry: ArrayLikeSeq[str | bytes | None], on_invalid: Literal["raise", "warn", "ignore", "fix"] = "raise", **kwargs
) -> GeoArray: ...
