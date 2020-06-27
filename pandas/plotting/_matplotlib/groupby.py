from typing import Dict, List, Optional, Union

import numpy as np

from pandas._typing import FrameOrSeriesUnion, Label

from pandas.core.dtypes.missing import isna

from pandas import DataFrame, MultiIndex, Series, concat


def create_iter_data_given_by(
    data: DataFrame, by: Optional[List[Label]]
) -> Union[DataFrame, Dict[str, FrameOrSeriesUnion]]:
    """
    Create data for iteration given `by` is assigned or not, and it is only
    used in both hist and boxplot.

    If `by` is assigned, return a dictionary of DataFrames in which the key of
    dictionary is the values in groups.
    If `by` is not assigned, return input as is, and this preserves current
    status of iter_data.

    Parameters
    ----------
    data: reformatted grouped data from `_compute_plot_data` method
    by: list or None, value assigned to `by`.

    Returns
    -------
    iter_data: DataFrame or Dictionary of DataFrames

    Examples
    --------
    If `by` is assigned:

    >>> import numpy as np
    >>> tuples = [('h1', 'a'), ('h1', 'b'), ('h2', 'a'), ('h2', 'b')]
    >>> mi = MultiIndex.from_tuples(tuples)
    >>> value = [[1, 3, np.nan, np.nan],
    ...          [3, 4, np.nan, np.nan], [np.nan, np.nan, 5, 6]]
    >>> data = DataFrame(value, columns=mi)
    >>> create_iter_data_given_by(data, by=["col1"])
    {'h1': DataFrame({'a': [1, 3, np.nan], 'b': [3, 4, np.nan]}),
     'h2': DataFrame({'a': [np.nan, np.nan, 5], 'b': [np.nan, np.nan, 6]})}
    """
    iter_data: Union[DataFrame, Dict[str, FrameOrSeriesUnion]]
    if not by:
        iter_data = data
    else:
        # Select sub-columns based on the value of first level of MI
        assert isinstance(data.columns, MultiIndex)
        cols = data.columns.levels[0]
        iter_data = {
            col: data.loc[:, data.columns.get_level_values(0) == col] for col in cols
        }
    return iter_data


def reconstruct_data_with_by(
    data: DataFrame, by: Union[Label, List[Label]], cols: List[Label]
) -> DataFrame:
    """
    Internal function to group data, and reassign multiindex column names onto the
    result in order to let grouped data be used in _compute_plot_data method.

    Parameters
    ----------
    data: Original DataFrame to plot
    by: grouped `by` parameter selected by users
    cols: columns of data set (excluding columns used in `by`)

    Returns
    -------
    Output is the reconstructed DataFrame with MultiIndex columns. The first level
    of MI is unique values of groups, and second level of MI is the columns
    selected by users.

    Examples
    --------
    >>> d = {'h': ['h1', 'h1', 'h2'], 'a': [1, 3, 5], 'b': [3, 4, 6]}
    >>> df = DataFrame(d)
    >>> reconstruct_data_with_by(df, by='h', cols=['a', 'b'])
       h1      h2
       a   b   a   b
    0  1   3   NaN NaN
    1  3   4   NaN NaN
    2  NaN NaN 5   6
    """
    grouped = data.groupby(by)

    data_list = []
    for key, group in grouped:
        columns = MultiIndex.from_product([[key], cols])
        sub_group = group[cols]
        sub_group.columns = columns
        data_list.append(sub_group)

    data = concat(data_list, axis=1)
    return data


def reformat_hist_y_given_by(
    y: Union[Series, np.array], by: Optional[Union[Label, List[Label]]]
) -> Union[Series, np.array]:
    """Internal function to reformat y given `by` is applied or not for hist plot.

    If by is None, input y is 1-d with NaN removed; and if by is not None, groupby
    will take place and input y is multi-dimensional array.
    """
    if by is not None and len(y.shape) > 1:
        notna = [col[~isna(col)] for col in y.T]
        y = np.array(np.array(notna).T)
    else:
        y = y[~isna(y)]
    return y
