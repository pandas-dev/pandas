from pandas import DataFrame
import numpy as np

def pivot_table(data, values=None, rows=None, cols=None, aggfunc=np.mean,
                fill_value=None):
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
    aggfunc : function, default numpy.mean
    fill_value : scalar, default None
        Value to replace missing values with

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

    keys = rows + cols
    grouped = data.groupby(keys)

    if values is not None:
        grouped = grouped[values]

    agged = grouped.agg(aggfunc)

    table = agged
    for k in cols:
        table = table.unstack(level=k)

    if fill_value is not None:
        table = table.fillna(value=fill_value)

    return table

def _convert_by(by):
    if by is None:
        by = []
    elif np.isscalar(by):
        by = [by]
    else:
        by = list(by)
    return by

def pprint_table(table):
    pass

if __name__ == '__main__':
    def _sample(values, n):
        indexer = np.random.randint(0, len(values), n)
        return np.asarray(values).take(indexer)

    levels = [['a', 'b', 'c', 'd'],
              ['foo', 'bar', 'baz'],
              ['one', 'two'],
              ['US', 'JP', 'UK']]
    names = ['k1', 'k2', 'k3', 'k4']

    n = 100000

    data = {}
    for name, level in zip(names, levels):
        data[name] = _sample(level, n)

    data['values'] = np.random.randn(n)
    data = DataFrame(data)

    table = pivot_table(data, values='values',
                        rows=['k1', 'k2'], cols=['k3', 'k4'])

