import numpy as np
import pytest

from pandas.compat import lrange, lzip

from pandas import DataFrame, MultiIndex, Series
from pandas.core import common as com
import pandas.util.testing as tm


def test_detect_chained_assignment():
    # Inplace ops, originally from:
    # http://stackoverflow.com/questions/20508968/series-fillna-in-a-multiindex-dataframe-does-not-fill-is-this-a-bug
    a = [12, 23]
    b = [123, None]
    c = [1234, 2345]
    d = [12345, 23456]
    tuples = [('eyes', 'left'), ('eyes', 'right'), ('ears', 'left'),
              ('ears', 'right')]
    events = {('eyes', 'left'): a,
              ('eyes', 'right'): b,
              ('ears', 'left'): c,
              ('ears', 'right'): d}
    multiind = MultiIndex.from_tuples(tuples, names=['part', 'side'])
    zed = DataFrame(events, index=['a', 'b'], columns=multiind)

    with pytest.raises(com.SettingWithCopyError):
        zed['eyes']['right'].fillna(value=555, inplace=True)


def test_cache_updating():
    # 5216
    # make sure that we don't try to set a dead cache
    a = np.random.rand(10, 3)
    df = DataFrame(a, columns=['x', 'y', 'z'])
    tuples = [(i, j) for i in range(5) for j in range(2)]
    index = MultiIndex.from_tuples(tuples)
    df.index = index

    # setting via chained assignment
    # but actually works, since everything is a view
    df.loc[0]['z'].iloc[0] = 1.
    result = df.loc[(0, 0), 'z']
    assert result == 1

    # correct setting
    df.loc[(0, 0), 'z'] = 2
    result = df.loc[(0, 0), 'z']
    assert result == 2


def test_indexer_caching():
    # GH5727
    # make sure that indexers are in the _internal_names_set
    n = 1000001
    arrays = [lrange(n), lrange(n)]
    index = MultiIndex.from_tuples(lzip(*arrays))
    s = Series(np.zeros(n), index=index)
    str(s)

    # setitem
    expected = Series(np.ones(n), index=index)
    s = Series(np.zeros(n), index=index)
    s[s == 0] = 1
    tm.assert_series_equal(s, expected)
