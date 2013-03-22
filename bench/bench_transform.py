import numpy as np
import pandas as pd
from pandas import Index, MultiIndex, DataFrame
from pandas.core.groupby import SeriesGroupBy, DataFrameGroupBy

def apply_by_group(grouped, f):
    """
    Applies a function to each Series or DataFrame in a GroupBy object, concatenates the results
    and returns the resulting Series or DataFrame.

    Parameters
    ----------
    grouped: SeriesGroupBy or DataFrameGroupBy
    f: callable
        Function to apply to each Series or DataFrame in the grouped object.

    Returns
    -------
    Series or DataFrame that results from applying the function to each Series or DataFrame in the
    GroupBy object and concatenating the results.

    """
    assert isinstance(grouped, (SeriesGroupBy, DataFrameGroupBy))
    assert hasattr(f, '__call__')

    groups = []
    for key, group in grouped:
        groups.append(f(group))
    c = pd.concat(groups)
    c.sort_index(inplace=True)
    return c

n_dates = 1000
n_securities = 2000
n_columns = 3
share_na = 0.1

dates = pd.date_range('1997-12-31', periods=n_dates, freq='B')
dates = Index(map(lambda x: x.year * 10000 + x.month * 100 + x.day, dates))

secid_min = int('10000000', 16)
secid_max = int('F0000000', 16)
step = (secid_max - secid_min) // (n_securities - 1)
security_ids = map(lambda x: hex(x)[2:10].upper(), range(secid_min, secid_max + 1, step))

data_index = MultiIndex(levels=[dates.values, security_ids],
    labels=[[i for i in xrange(n_dates) for _ in xrange(n_securities)], range(n_securities) * n_dates],
    names=['date', 'security_id'])
n_data = len(data_index)

columns = Index(['factor{}'.format(i) for i in xrange(1, n_columns + 1)])

data = DataFrame(np.random.randn(n_data, n_columns), index=data_index, columns=columns)

step = int(n_data * share_na)
for column_index in xrange(n_columns):
    index = column_index
    while index < n_data:
        data.set_value(data_index[index], columns[column_index], np.nan)
        index += step

grouped = data.groupby(level='security_id')
f_fillna = lambda x: x.fillna(method='pad')

#%timeit grouped.transform(f_fillna)
#%timeit apply_by_group(grouped, f_fillna)
