from typing import TYPE_CHECKING, Any, Hashable, Iterable, Sequence

from pandas.wesm import dataframe as dataframe_protocol

if TYPE_CHECKING:
    from pandas import DataFrame


class PandasDataFrame(dataframe_protocol.DataFrame):
    """
    Construct generic data frame from pandas DataFrame

    Parameters
    ----------
    df : pd.DataFrame
    """

    def __init__(self, df: "DataFrame"):
        self._df = df

    def column_by_index(self, i: int) -> dataframe_protocol.Column:
        """
        Return the column at the indicated position.
        """
        pass

    def column_by_name(self, key: Hashable) -> dataframe_protocol.Column:
        """
        Return the column whose name is the indicated key.
        """
        pass

    @property
    def column_names(self) -> Sequence[Any]:
        """
        Return the column names as a materialized sequence.
        """
        pass

    def iter_column_names(self) -> Iterable[Any]:
        """
        Return the column names as an iterable.
        """
        pass

    @property
    def num_columns(self) -> int:
        """
        Return the number of columns in the DataFrame.
        """
        pass

    @property
    def num_rows(self) -> int:
        """
        Return the number of rows in the DataFrame.
        """
        pass
