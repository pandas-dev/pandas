# pylint: disable=E1103

from pandas import Series, DataFrame
from pandas.core.index import MultiIndex
from pandas.tools.merge import concat
from pandas.tools.util import cartesian_product
from pandas.compat import range, lrange, zip
from pandas import compat
import pandas.core.common as com
import numpy as np


def pivot_table(data, values=None, rows=None, cols=None, aggfunc='mean',
                fill_value=None, margins=False, dropna=True):
    """
    Create a spreadsheet-style pivot table as a DataFrame. The levels in the
    pivot table will be stored in MultiIndex objects (hierarchical indexes) on
    the index and columns of the result DataFrame

    Parameters
    ----------
    data : DataFrame
    values : column to aggregate, optional
    rows : list of column names or arrays to group on
        Keys to group on the x-axis of the pivot table
    cols : list of column names or arrays to group on
        Keys to group on the y-axis of the pivot table
    aggfunc : function, default numpy.mean, or list of functions
        If list of functions passed, the resulting pivot table will have
        hierarchical columns whose top level are the function names (inferred
        from the function objects themselves)
    fill_value : scalar, default None
        Value to replace missing values with
    margins : boolean, default False
        Add all row / columns (e.g. for subtotal / grand totals)
    dropna : boolean, default True
        Do not include columns whose entries are all NaN

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
        to_filter = []
        for x in keys + values:
            try:
                if x in data:
                    to_filter.append(x)
            except TypeError:
                pass
        if len(to_filter) < len(data.columns):
            data = data[to_filter]

    grouped = data.groupby(keys)
    agged = grouped.agg(aggfunc)

    table = agged
    if table.index.nlevels > 1:
        to_unstack = [agged.index.names[i]
                      for i in range(len(rows), len(keys))]
        table = agged.unstack(to_unstack)

    if not dropna:
        try:
            m = MultiIndex.from_arrays(cartesian_product(table.index.levels))
            table = table.reindex_axis(m, axis=0)
        except AttributeError:
            pass # it's a single level

        try:
            m = MultiIndex.from_arrays(cartesian_product(table.columns.levels))
            table = table.reindex_axis(m, axis=1)
        except AttributeError:
            pass # it's a single level or a series

    if isinstance(table, DataFrame):
        if isinstance(table.columns, MultiIndex):
            table = table.sortlevel(axis=1)
        else:
            table = table.sort_index(axis=1)

    if fill_value is not None:
        table = table.fillna(value=fill_value, downcast='infer')

    if margins:
        table = _add_margins(table, data, values, rows=rows,
                             cols=cols, aggfunc=aggfunc)

    # discard the top level
    if values_passed and not values_multi:
        table = table[values[0]]

    if len(rows) == 0 and len(cols) > 0:
        table = table.T

    return table


DataFrame.pivot_table = pivot_table


def _add_margins(table, data, values, rows, cols, aggfunc):

    grand_margin = _compute_grand_margin(data, values, aggfunc)

    if not values and isinstance(table, Series):
        # If there are no values and the table is a series, then there is only
        # one column in the data. Compute grand margin and return it.
        row_key = ('All',) + ('',) * (len(rows) - 1) if len(rows) > 1 else 'All'
        return table.append(Series({row_key: grand_margin['All']}))

    if values:
        marginal_result_set = _generate_marginal_results(table, data, values, rows, cols, aggfunc, grand_margin)
        if not isinstance(marginal_result_set, tuple):
            return marginal_result_set
        result, margin_keys, row_margin = marginal_result_set
    else:
        marginal_result_set = _generate_marginal_results_without_values(table, data, rows, cols, aggfunc)
        if not isinstance(marginal_result_set, tuple):
            return marginal_result_set
        result, margin_keys, row_margin = marginal_result_set

    key = ('All',) + ('',) * (len(rows) - 1) if len(rows) > 1 else 'All'

    row_margin = row_margin.reindex(result.columns)
    # populate grand margin
    for k in margin_keys:
        if isinstance(k, compat.string_types):
            row_margin[k] = grand_margin[k]
        else:
            row_margin[k] = grand_margin[k[0]]

    margin_dummy = DataFrame(row_margin, columns=[key]).T

    row_names = result.index.names
    result = result.append(margin_dummy)
    result.index.names = row_names

    return result


def _compute_grand_margin(data, values, aggfunc):

    if values:
        grand_margin = {}
        for k, v in data[values].iteritems():
            try:
                if isinstance(aggfunc, compat.string_types):
                    grand_margin[k] = getattr(v, aggfunc)()
                else:
                    grand_margin[k] = aggfunc(v)
            except TypeError:
                pass
        return grand_margin
    else:
        return {'All': aggfunc(data.index)}


def _generate_marginal_results(table, data, values, rows, cols, aggfunc, grand_margin):
    if len(cols) > 0:
        # need to "interleave" the margins
        table_pieces = []
        margin_keys = []

        def _all_key(key):
            return (key, 'All') + ('',) * (len(cols) - 1)

        if len(rows) > 0:
            margin = data[rows + values].groupby(rows).agg(aggfunc)
            cat_axis = 1
            for key, piece in table.groupby(level=0, axis=cat_axis):
                all_key = _all_key(key)
                piece[all_key] = margin[key]
                table_pieces.append(piece)
                margin_keys.append(all_key)
        else:
            margin = grand_margin
            cat_axis = 0
            for key, piece in table.groupby(level=0, axis=cat_axis):
                all_key = _all_key(key)
                table_pieces.append(piece)
                table_pieces.append(Series(margin[key], index=[all_key]))
                margin_keys.append(all_key)

        result = concat(table_pieces, axis=cat_axis)

        if len(rows) == 0:
            return result
    else:
        result = table
        margin_keys = table.columns

    if len(cols) > 0:
        row_margin = data[cols + values].groupby(cols).agg(aggfunc)
        row_margin = row_margin.stack()

        # slight hack
        new_order = [len(cols)] + lrange(len(cols))
        row_margin.index = row_margin.index.reorder_levels(new_order)
    else:
        row_margin = Series(np.nan, index=result.columns)

    return result, margin_keys, row_margin


def _generate_marginal_results_without_values(table, data, rows, cols, aggfunc):
    if len(cols) > 0:
        # need to "interleave" the margins
        margin_keys = []

        def _all_key():
            if len(cols) == 1:
                return 'All'
            return ('All', ) + ('', ) * (len(cols) - 1)

        if len(rows) > 0:
            margin = data[rows].groupby(rows).apply(aggfunc)
            all_key = _all_key()
            table[all_key] = margin
            result = table
            margin_keys.append(all_key)

        else:
            margin = data.groupby(level=0, axis=0).apply(aggfunc)
            all_key = _all_key()
            table[all_key] = margin
            result = table
            margin_keys.append(all_key)
            return result
    else:
        result = table
        margin_keys = table.columns

    if len(cols):
        row_margin = data[cols].groupby(cols).apply(aggfunc)
    else:
        row_margin = Series(np.nan, index=result.columns)

    return result, margin_keys, row_margin


def _convert_by(by):
    if by is None:
        by = []
    elif (np.isscalar(by) or isinstance(by, (np.ndarray, Series))
          or hasattr(by, '__call__')):
        by = [by]
    else:
        by = list(by)
    return by


def crosstab(rows, cols, values=None, rownames=None, colnames=None,
             aggfunc=None, margins=False, dropna=True):
    """
    Compute a simple cross-tabulation of two (or more) factors. By default
    computes a frequency table of the factors unless an array of values and an
    aggregation function are passed

    Parameters
    ----------
    rows : array-like, Series, or list of arrays/Series
        Values to group by in the rows
    cols : array-like, Series, or list of arrays/Series
        Values to group by in the columns
    values : array-like, optional
        Array of values to aggregate according to the factors
    aggfunc : function, optional
        If no values array is passed, computes a frequency table
    rownames : sequence, default None
        If passed, must match number of row arrays passed
    colnames : sequence, default None
        If passed, must match number of column arrays passed
    margins : boolean, default False
        Add row/column margins (subtotals)
    dropna : boolean, default True
        Do not include columns whose entries are all NaN

    Notes
    -----
    Any Series passed will have their name attributes used unless row or column
    names for the cross-tabulation are specified

    Examples
    --------
    >>> a
    array([foo, foo, foo, foo, bar, bar,
           bar, bar, foo, foo, foo], dtype=object)
    >>> b
    array([one, one, one, two, one, one,
           one, two, two, two, one], dtype=object)
    >>> c
    array([dull, dull, shiny, dull, dull, shiny,
           shiny, dull, shiny, shiny, shiny], dtype=object)

    >>> crosstab(a, [b, c], rownames=['a'], colnames=['b', 'c'])
    b    one          two
    c    dull  shiny  dull  shiny
    a
    bar  1     2      1     0
    foo  2     2      1     2

    Returns
    -------
    crosstab : DataFrame
    """
    rows = com._maybe_make_list(rows)
    cols = com._maybe_make_list(cols)

    rownames = _get_names(rows, rownames, prefix='row')
    colnames = _get_names(cols, colnames, prefix='col')

    data = {}
    data.update(zip(rownames, rows))
    data.update(zip(colnames, cols))

    if values is None:
        df = DataFrame(data)
        df['__dummy__'] = 0
        table = df.pivot_table('__dummy__', rows=rownames, cols=colnames,
                               aggfunc=len, margins=margins, dropna=dropna)
        return table.fillna(0).astype(np.int64)
    else:
        data['__dummy__'] = values
        df = DataFrame(data)
        table = df.pivot_table('__dummy__', rows=rownames, cols=colnames,
                               aggfunc=aggfunc, margins=margins, dropna=dropna)
        return table


def _get_names(arrs, names, prefix='row'):
    if names is None:
        names = []
        for i, arr in enumerate(arrs):
            if isinstance(arr, Series) and arr.name is not None:
                names.append(arr.name)
            else:
                names.append('%s_%d' % (prefix, i))
    else:
        if len(names) != len(arrs):
            raise AssertionError('arrays and names must have the same length')
        if not isinstance(names, list):
            names = list(names)

    return names
