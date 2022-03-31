import collections.abc
import ctypes

from typing import Tuple, Any

from pandas.core.exchange.dataframe_protocol import (
    DataFrame as DataFrameXchg,
    DtypeKind,
)

import pandas as pd

from pandas.core.exchange.column import (
    convert_column_to_ndarray,
    convert_categorical_column,
    convert_string_column,
    PandasColumn,
)


def from_dataframe(df: DataFrameXchg, allow_copy: bool = True) -> pd.DataFrame:
    """
    Construct a pandas DataFrame from ``df`` if it supports ``__dataframe__``
    """
    if isinstance(df, pd.DataFrame):
        return df

    if not hasattr(df, "__dataframe__"):
        raise ValueError("`df` does not support __dataframe__")

    return _from_dataframe(df.__dataframe__(allow_copy=allow_copy))


def _from_dataframe(df: DataFrameXchg) -> pd.DataFrame:
    """
    Note: not all cases are handled yet, only ones that can be implemented with
    only Pandas. Later, we need to implement/test support for categoricals,
    bit/byte masks, chunk handling, etc.
    """
    _buffers = []  # hold on to buffers, keeps memory alive
    result = []
    for chunk in df.get_chunks():
        # We need a dict of columns here, with each column being a numpy array (at
        # least for now, deal with non-numpy dtypes later).
        chunk_cols = {}
        for name in chunk.column_names():
            if not isinstance(name, str):
                raise ValueError(f"Column {name} is not a string")
            if name in chunk_cols:
                raise ValueError(f"Column {name} is not unique")
            col = chunk.get_column_by_name(name)
            if col.dtype[0] in (
                DtypeKind.INT,
                DtypeKind.UINT,
                DtypeKind.FLOAT,
                DtypeKind.BOOL,
            ):
                # Simple numerical or bool dtype, turn into numpy array
                chunk_cols[name], _buf = convert_column_to_ndarray(col)
            elif col.dtype[0] == DtypeKind.CATEGORICAL:
                chunk_cols[name], _buf = convert_categorical_column(col)
            elif col.dtype[0] == DtypeKind.STRING:
                chunk_cols[name], _buf = convert_string_column(col)
            else:
                raise NotImplementedError(f"Data type {col.dtype[0]} not handled yet")

            _buffers.append(_buf)

        df_new = pd.DataFrame(chunk_cols)
        result.append(df_new)

    df_new = pd.concat(result)
    df_new._buffers = _buffers
    return df_new


class PandasDataFrameXchg(DataFrameXchg):
    """
    A data frame class, with only the methods required by the interchange
    protocol defined.
    Instances of this (private) class are returned from
    ``pd.DataFrame.__dataframe__`` as objects with the methods and
    attributes defined on this class.
    """

    def __init__(
        self, df: pd.DataFrame, nan_as_null: bool = False, allow_copy: bool = True
    ) -> None:
        """
        Constructor - an instance of this (private) class is returned from
        `pd.DataFrame.__dataframe__`.
        """
        self._df = df
        # ``nan_as_null`` is a keyword intended for the consumer to tell the
        # producer to overwrite null values in the data with ``NaN`` (or ``NaT``).
        # This currently has no effect; once support for nullable extension
        # dtypes is added, this value should be propagated to columns.
        self._nan_as_null = nan_as_null
        self._allow_copy = allow_copy

    @property
    def metadata(self):
        # `index` isn't a regular column, and the protocol doesn't support row
        # labels - so we export it as Pandas-specific metadata here.
        return {"pandas.index": self._df.index}

    def num_columns(self) -> int:
        return len(self._df.columns)

    def num_rows(self) -> int:
        return len(self._df)

    def num_chunks(self) -> int:
        return 1

    def column_names(self):
        return self._df.columns.tolist()

    def get_column(self, i: int) -> PandasColumn:
        return PandasColumn(self._df.iloc[:, i], allow_copy=self._allow_copy)

    def get_column_by_name(self, name: str) -> PandasColumn:
        return PandasColumn(self._df[name], allow_copy=self._allow_copy)

    def get_columns(self):
        return [
            PandasColumn(self._df[name], allow_copy=self._allow_copy)
            for name in self._df.columns
        ]

    def select_columns(self, indices):
        if not isinstance(indices, collections.abc.Sequence):
            raise ValueError("`indices` is not a sequence")
        if not isinstance(indices, list):
            indices = list(indices)

        return PandasDataFrameXchg(
            self._df.iloc[:, indices], self._nan_as_null, self._allow_copy
        )

    def select_columns_by_name(self, names):
        if not isinstance(names, collections.abc.Sequence):
            raise ValueError("`names` is not a sequence")
        if not isinstance(names, list):
            names = list(names)

        return PandasDataFrameXchg(
            self._df.loc[:, names], self._nan_as_null, self._allow_copy
        )

    def get_chunks(self, n_chunks=None):
        """
        Return an iterator yielding the chunks.
        """
        if n_chunks and n_chunks > 1:
            size = len(self._df)
            step = size // n_chunks
            if size % n_chunks != 0:
                step += 1
            for start in range(0, step * n_chunks, step):
                yield PandasDataFrameXchg(
                    self._df.iloc[start : start + step, :],
                    self._nan_as_null,
                    self._allow_copy,
                )
        else:
            yield self
