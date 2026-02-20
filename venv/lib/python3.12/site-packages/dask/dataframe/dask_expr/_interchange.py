from __future__ import annotations

from pandas.core.interchange.dataframe_protocol import DtypeKind

from dask.dataframe._compat import is_string_dtype
from dask.dataframe.dask_expr._collection import DataFrame
from dask.dataframe.dispatch import is_categorical_dtype

_NP_KINDS = {
    "i": DtypeKind.INT,
    "u": DtypeKind.UINT,
    "f": DtypeKind.FLOAT,
    "b": DtypeKind.BOOL,
    "U": DtypeKind.STRING,
    "M": DtypeKind.DATETIME,
    "m": DtypeKind.DATETIME,
}


class DaskDataFrameInterchange:
    def __init__(
        self, df: DataFrame, nan_as_null: bool = False, allow_copy: bool = True
    ) -> None:
        self._df = df
        self._nan_as_null = nan_as_null
        self._allow_copy = allow_copy

    def get_columns(self):
        return [DaskColumn(self._df[name]) for name in self._df.columns]

    def column_names(self):
        return self._df.columns

    def num_columns(self) -> int:
        return len(self._df.columns)

    def num_rows(self) -> int:
        return len(self._df)


class DaskColumn:
    def __init__(self, column, allow_copy: bool = True) -> None:
        self._col = column
        self._allow_copy = allow_copy

    def dtype(self) -> tuple[DtypeKind, None, None, None]:
        dtype = self._col.dtype

        if is_categorical_dtype(dtype):
            return (
                DtypeKind.CATEGORICAL,
                None,
                None,
                None,
            )
        elif is_string_dtype(dtype):
            return (
                DtypeKind.STRING,
                None,
                None,
                None,
            )
        else:
            return _NP_KINDS.get(dtype.kind, None), None, None, None
