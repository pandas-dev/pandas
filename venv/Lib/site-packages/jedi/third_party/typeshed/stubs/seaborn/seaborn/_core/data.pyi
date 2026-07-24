from _typeshed import Incomplete
from collections.abc import Mapping
from typing import TypeVar, overload

from pandas import DataFrame
from seaborn._core.typing import DataSource, SupportsDataFrame, VariableSpec

_T = TypeVar("_T", Mapping[Incomplete, Incomplete], None)

class PlotData:
    frame: DataFrame
    frames: dict[tuple[str, str], DataFrame]
    names: dict[str, str | None]
    ids: dict[str, str | int]
    source_data: DataSource
    source_vars: dict[str, VariableSpec]
    def __init__(self, data: DataSource, variables: dict[str, VariableSpec]) -> None: ...
    def __contains__(self, key: str) -> bool: ...
    def join(self, data: DataSource, variables: dict[str, VariableSpec] | None) -> PlotData: ...

@overload
def handle_data_source(data: _T) -> _T: ...
@overload
def handle_data_source(data: SupportsDataFrame) -> DataFrame: ...
def convert_dataframe_to_pandas(data: object) -> DataFrame: ...
