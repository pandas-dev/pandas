from typing import TYPE_CHECKING, Any, Hashable, Iterable, Optional, Sequence

from pandas.wesm import dataframe as dataframe_protocol

if TYPE_CHECKING:
    import pandas as pd
    import numpy as np


class Column(dataframe_protocol.Column):
    """
    Construct generic column from pandas Series

    Parameters
    ----------
    ser : pd.Series
    """

    _ser: "pd.Series"

    def __init__(self, ser: "pd.Series"):
        self._ser = ser

    @property
    def name(self) -> Optional[Hashable]:
        return self._ser.name

    @property
    def type(self) -> dataframe_protocol.DataType:
        """
        Return the logical type of each column cell value
        """
        raise NotImplementedError

    def to_numpy(self) -> "np.ndarray":
        """
        Access column's data as a NumPy array. Recommended to return a view if
        able but not required
        """
        return self._ser.to_numpy()

    def to_arrow(self, **kwargs):
        """
        Access column's data in the Apache Arrow format as pyarrow.Array or
        ChunkedArray. Recommended to return a view if able but not required
        """
        raise NotImplementedError("Conversion to Arrow not available")


class DataFrame(dataframe_protocol.DataFrame):
    """
    Construct generic data frame from pandas DataFrame

    Parameters
    ----------
    df : pd.DataFrame
    """

    _df: "pd.DataFrame"

    def __init__(self, df: "pd.DataFrame"):
        self._df = df

    def __str__(self) -> str:
        return str(self._df)

    def __repr__(self) -> str:
        return repr(self._df)

    def column_by_index(self, i: int) -> dataframe_protocol.Column:
        """
        Return the column at the indicated position.
        """
        return Column(self._df.iloc[:, i])

    def column_by_name(self, key: Hashable) -> dataframe_protocol.Column:
        """
        Return the column whose name is the indicated key.
        """
        return Column(self._df[key])

    @property
    def column_names(self) -> Sequence[Any]:
        """
        Return the column names as a materialized sequence.
        """
        return self._df.columns.to_list()

    @property
    def row_names(self) -> Sequence[Any]:
        """
        Return the row names (if any) as a materialized sequence. It is not
        necessary to implement this method
        """
        return self._df.index.to_list()

    def iter_column_names(self) -> Iterable[Any]:
        """
        Return the column names as an iterable.
        """
        return self.column_names

    @property
    def num_columns(self) -> int:
        """
        Return the number of columns in the DataFrame.
        """
        return self._df.shape[1]

    @property
    def num_rows(self) -> int:
        """
        Return the number of rows in the DataFrame.
        """
        return len(self._df)
