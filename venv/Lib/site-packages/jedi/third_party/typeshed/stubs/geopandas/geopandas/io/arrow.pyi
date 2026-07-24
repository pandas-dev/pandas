import os
from _typeshed import Incomplete, SupportsGetItem, SupportsKeysAndGetItem
from collections.abc import Iterable, Mapping
from typing import Any, Final

from ..geodataframe import GeoDataFrame

METADATA_VERSION: Final[str]
SUPPORTED_VERSIONS_LITERAL = Incomplete
SUPPORTED_VERSIONS: Final[list[str]]
GEOARROW_ENCODINGS: Final[list[str]]
SUPPORTED_ENCODINGS: Final[list[str]]
PARQUET_GEOMETRY_ENCODINGS = Incomplete

def _read_parquet(
    path: str | os.PathLike[str],
    columns: Iterable[str] | None = None,
    storage_options: SupportsKeysAndGetItem[str, Any] | None = None,  # type depend on the connection
    bbox: SupportsGetItem[int, float] | None = None,
    to_pandas_kwargs: Mapping[str, Incomplete] | None = None,
    **kwargs,  # kwargs passed to pyarrow.parquet.read_table
) -> GeoDataFrame: ...
def _read_feather(
    path: str | os.PathLike[str],
    columns: Iterable[str] | None = None,
    to_pandas_kwargs: Mapping[str, Incomplete] | None = None,
    **kwargs,  # kwargs passed to pyarrow.feather.read_table
) -> GeoDataFrame: ...
