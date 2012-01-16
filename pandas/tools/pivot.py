# pylint: disable=E1103

from pandas import Series, DataFrame
from pandas.tools.merge import concat
import numpy as np

def pivot_table(data, values=None, rows=None, cols=None, aggfunc=np.mean,
                fill_value=None, margins=False):
    """
    Create a spreadsheet-style pivot table as a DataFrame. The levels in the
    pivot table will be stored in MultiIndex objects (hierarchical indexes) on
    the index and columns of the result DataFrame

    Parameters
    ----------
    data : DataFrame
    values : column to aggregate, optional
    rows : list
        Columns to group on the x-axis of the pivot table
    cols : list
        Columns to group on the x-axis of the pivot table
    aggfunc : function, default numpy.mean, or list of functions
        If list of functions passed, the resulting pivot table will have
        hierarchical columns whose top level are the function names (inferred
        from the function objects themselves)
    fill_value : scalar, default None
        Value to replace missing values with
    margins : boolean, default False
        Add all row / columns (e.g. for subtotal / grand totals)

    Examples
    --------
    >>> df
       A   B   C      D
    0  foo one small  1
    1  foo one large  2
    2  foo one large  2
    3  foo two small  3
    4  foo two small  3
    5  bar one large  4
    6  bar one small  5
    7  bar two small  6
    8  bar two large  7

    >>> table = pivot_table(df, values='D', rows=['A', 'B'],
    ...                     cols=['C'], aggfunc=np.sum)
    >>> table
              small  large
    foo  one  1      4
         two  6      NaN
    bar  one  5      4
         two  6      7

    Returns
    -------
    table : DataFrame
    """
    rows = _convert_by(rows)
    cols = _convert_by(cols)

    if isinstance(aggfunc, list):
        pieces = []
        keys = []
        for func in aggfunc:
            table = pivot_table(data, values=values, rows=rows, cols=cols,
                                fill_value=fill_value, aggfunc=func,
                                margins=margins)
            pieces.append(table)
            keys.append(func.__name__)
        return concat(pieces, keys=keys, axis=1)

    keys = rows + cols

    values_passed = values is not None
    if values_passed:
        if isinstance(values, (list, tuple)):
            values_multi = True
        else:
            values_multi = False
            values = [values]
    else:
        values = list(data.columns.drop(keys))

    if values_passed:
        data = data[keys + values]

    grouped = data.groupby(keys)
    agged = grouped.agg(aggfunc)

    table = agged
    for k in cols:
        table = table.unstack(level=k)

    if fill_value is not None:
        table = table.fillna(value=fill_value)

    if margins:
        table = _add_margins(table, data, values, rows=rows,
                             cols=cols, aggfunc=aggfunc)

    # discard the top level
    if values_passed and not values_multi:
        table = table[values[0]]

    return table

DataFrame.pivot_table = pivot_table

def _add_margins(table, data, values, rows=None, cols=None, aggfunc=np.mean):
    if len(cols) > 0:
        col_margin = data[rows + values].groupby(rows).agg(aggfunc)

        # need to "interleave" the margins
        table_pieces = []
        margin_keys = []
        for key, piece in table.groupby(level=0, axis=1):
            all_key = (key, 'All') + ('',) * (len(cols) - 1)
            piece[all_key] = col_margin[key]
            table_pieces.append(piece)
            margin_keys.append(all_key)

        result = concat(table_pieces, axis=1)
    else:
        result = table
        margin_keys = table.columns

    grand_margin = {}
    for k, v in data[values].iteritems():
        try:
            grand_margin[k] = aggfunc(v)
        except TypeError:
            pass

    if len(cols) > 0:
        row_margin = data[cols + values].groupby(cols).agg(aggfunc)
        row_margin = row_margin.stack()

        # slight hack
        new_order = [len(cols)] + range(len(cols))
        row_margin.index = row_margin.index.reorder_levels(new_order)
    else:
        row_margin = Series(np.nan, index=result.columns)

    key = ('All',) + ('',) * (len(rows) - 1)

    row_margin = row_margin.reindex(result.columns)
    # populate grand margin
    for k in margin_keys:
        if len(cols) > 0:
            row_margin[k] = grand_margin[k[0]]
        else:
            row_margin[k] = grand_margin[k]

    margin_dummy = DataFrame(row_margin, columns=[key]).T
    result = result.append(margin_dummy)

    return result

def _convert_by(by):
    if by is None:
        by = []
    elif np.isscalar(by):
        by = [by]
    else:
        by = list(by)
    return by

def crosstab(rows, columns):
    """
    Compute a simple cross-tabulation of two (or more) factors

    Parameters
    ----------
    rows :
    columns :

    Returns
    -------
    crosstab : DataFrame
    """
    rname = cname = None
    if isinstance(rows, Series):
        rname = rows.name

    if isinstance(columns, Series):
        cname = columns.name

    df = DataFrame({'rows' : rows, 'columns' : columns})
    table = df.groupby(['rows', 'columns']).size()

    result = table.unstack()
    result.columns.name = cname
    result.index.name = rname

    return result
