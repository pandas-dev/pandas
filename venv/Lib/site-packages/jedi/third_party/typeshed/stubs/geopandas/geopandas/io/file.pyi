import os
from _typeshed import Incomplete, SupportsRead
from collections import OrderedDict
from typing import Literal, TypedDict, overload, type_check_only

import pandas as pd
from pandas._typing import Axes

from ..base import _BboxLike, _MaskLike
from ..geodataframe import GeoDataFrame

# Keep inline with GeoDataFrame.from_file and GeoSeries.from_file
@overload
def _read_file(
    filename: str | os.PathLike[str] | SupportsRead[Incomplete],
    bbox: _BboxLike | None = None,
    mask: _MaskLike | None = None,
    columns: Axes | None = None,
    rows: int | slice | None = None,
    engine: Literal["fiona", "pyogrio"] | None = None,
    *,
    ignore_geometry: Literal[False] = False,
    layer: int | str | None = None,
    encoding: str | None = None,
    **kwargs,  # depend on engine
) -> GeoDataFrame: ...
@overload
def _read_file(
    filename: str | os.PathLike[str] | SupportsRead[Incomplete],
    bbox: _BboxLike | None = None,
    mask: _MaskLike | None = None,
    columns: Axes | None = None,
    rows: int | slice | None = None,
    engine: Literal["fiona", "pyogrio"] | None = None,
    *,
    ignore_geometry: Literal[True],
    layer: int | str | None = None,
    encoding: str | None = None,
    **kwargs,  # depend on engine
) -> pd.DataFrame: ...
@type_check_only
class _Schema(TypedDict):
    geometry: str | list[str]
    properties: OrderedDict[str, str]

def infer_schema(df: GeoDataFrame) -> _Schema: ...
def _list_layers(filename: str | bytes | os.PathLike[str] | os.PathLike[bytes] | SupportsRead[Incomplete]) -> pd.DataFrame: ...
