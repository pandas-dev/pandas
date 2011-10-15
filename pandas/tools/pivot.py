from pandas import DataFrame
import numpy as np

def pivot_table(data, values=None, xby=None, yby=None, aggfunc=np.mean,
                fill_value=None):
    """

    """

    xby = [] if xby is None else list(xby)
    yby = [] if yby is None else list(yby)

    keys = xby + yby
    grouped = data.groupby(keys)

    if values is not None:
        grouped = grouped[values]

    agged = grouped.agg(aggfunc)

    table = agged
    for k in yby:
        table = table.unstack(level=k)

    if fill_value is not None:
        table = table.fillna(value=fill_value)

    return table

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

    table = pivot_table(data, values='values', xby=['k1', 'k2'], yby=['k3', 'k4'])

