from typing import TYPE_CHECKING, Any, Hashable, Iterable, Sequence

from pandas.wesm import dataframe as dataframe_protocol
from pandas.wesm.example_dict_of_ndarray import NumPyColumn

if TYPE_CHECKING:
    import pandas as pd


class Column(NumPyColumn):
    """
    Construct generic column from pandas Series

    Parameters
    ----------
    ser : pd.Series
    """

    _ser: "pd.Series"

    def __init__(self, ser: "pd.Series"):
        self._ser = ser
        super().__init__(ser.name, ser.to_numpy())


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
