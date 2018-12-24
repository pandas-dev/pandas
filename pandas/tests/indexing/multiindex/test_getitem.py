import numpy as np
import pytest

from pandas.compat import StringIO, lrange, range, u, zip

import pandas as pd
from pandas import DataFrame, Index, MultiIndex, Series, date_range
import pandas.core.common as com
from pandas.util import testing as tm


@pytest.mark.parametrize('access_method', [lambda s, x: s[:, x],
                                           lambda s, x: s.loc[:, x],
                                           lambda s, x: s.xs(x, level=1)])
@pytest.mark.parametrize('level1_value, expected', [
    (0, Series([1], index=[0])),
    (1, Series([2, 3], index=[1, 2]))
])
def test_series_getitem_multiindex(access_method, level1_value, expected):

    # GH 6018
    # series regression getitem with a multi-index

    s = Series([1, 2, 3])
    s.index = MultiIndex.from_tuples([(0, 0), (1, 1), (2, 1)])
    result = access_method(s, level1_value)
    tm.assert_series_equal(result, expected)


def test_series_getitem_multiindex_xs():
    # GH6258
    dt = list(date_range('20130903', periods=3))
    idx = MultiIndex.from_product([list('AB'), dt])
    s = Series([1, 3, 4, 1, 3, 4], index=idx)

    result = s.xs('20130903', level=1)
    expected = Series([1, 1], index=list('AB'))
    tm.assert_series_equal(result, expected)


def test_series_getitem_multiindex_xs_by_label():
    # GH5684
    idx = MultiIndex.from_tuples([('a', 'one'), ('a', 'two'), ('b', 'one'),
                                  ('b', 'two')])
    s = Series([1, 2, 3, 4], index=idx)
    s.index.set_names(['L1', 'L2'], inplace=True)
    result = s.xs('one', level='L2')
    expected = Series([1, 3], index=['a', 'b'])
    expected.index.set_names(['L1'], inplace=True)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize('level0_value', ['D', 'A'])
def test_getitem_duplicates_multiindex(level0_value):
    # GH 5725 the 'A' happens to be a valid Timestamp so the doesn't raise
    # the appropriate error, only in PY3 of course!

    index = MultiIndex(levels=[[level0_value, 'B', 'C'],
                               [0, 26, 27, 37, 57, 67, 75, 82]],
                       codes=[[0, 0, 0, 1, 2, 2, 2, 2, 2, 2],
                              [1, 3, 4, 6, 0, 2, 2, 3, 5, 7]],
                       names=['tag', 'day'])
    arr = np.random.randn(len(index), 1)
    df = DataFrame(arr, index=index, columns=['val'])

    # confirm indexing on missing value raises KeyError
    if level0_value != 'A':
        msg = "'A'"
        with pytest.raises(KeyError, match=msg):
            df.val['A']

    msg = "'X'"
    with pytest.raises(KeyError, match=msg):
        df.val['X']

    result = df.val[level0_value]
    expected = Series(arr.ravel()[0:3], name='val', index=Index(
        [26, 37, 57], name='day'))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize('indexer, is_level1, expected_error', [
    ([], False, None),  # empty ok
    (['A'], False, None),
    (['A', 'D'], False, None),
    (['D'], False, r"\['D'\] not in index"),  # not any values found
    (pd.IndexSlice[:, ['foo']], True, None),
    (pd.IndexSlice[:, ['foo', 'bah']], True, None)
])
def test_getitem_duplicates_multiindex_missing_indexers(indexer, is_level1,
                                                        expected_error):
    # GH 7866
    # multi-index slicing with missing indexers
    idx = MultiIndex.from_product([['A', 'B', 'C'],
                                   ['foo', 'bar', 'baz']],
                                  names=['one', 'two'])
    s = Series(np.arange(9, dtype='int64'), index=idx).sort_index()

    if indexer == []:
        expected = s.iloc[[]]
    elif is_level1:
        expected = Series([0, 3, 6], index=MultiIndex.from_product(
            [['A', 'B', 'C'], ['foo']], names=['one', 'two'])).sort_index()
    else:
        exp_idx = MultiIndex.from_product([['A'], ['foo', 'bar', 'baz']],
                                          names=['one', 'two'])
        expected = Series(np.arange(3, dtype='int64'),
                          index=exp_idx).sort_index()

    if expected_error is not None:
        with pytest.raises(KeyError, match=expected_error):
            s.loc[indexer]
    else:
        result = s.loc[indexer]
        tm.assert_series_equal(result, expected)


@pytest.mark.parametrize('columns_indexer', [
    ([], slice(None)),
    (['foo'], [])
])
def test_getitem_duplicates_multiindex_empty_indexer(columns_indexer):
    # GH 8737
    # empty indexer
    multi_index = MultiIndex.from_product((['foo', 'bar', 'baz'],
                                           ['alpha', 'beta']))
    df = DataFrame(np.random.randn(5, 6), index=range(5), columns=multi_index)
    df = df.sort_index(level=0, axis=1)

    expected = DataFrame(index=range(5), columns=multi_index.reindex([])[0])
    result = df.loc[:, columns_indexer]
    tm.assert_frame_equal(result, expected)


def test_getitem_duplicates_multiindex_non_scalar_type_object():
    # regression from < 0.14.0
    # GH 7914
    df = DataFrame([[np.mean, np.median], ['mean', 'median']],
                   columns=MultiIndex.from_tuples([('functs', 'mean'),
                                                   ('functs', 'median')]),
                   index=['function', 'name'])
    result = df.loc['function', ('functs', 'mean')]
    expected = np.mean
    assert result == expected


def test_getitem_simple(multiindex_dataframe_random_data):
    frame = multiindex_dataframe_random_data
    df = frame.T

    col = df['foo', 'one']
    tm.assert_almost_equal(col.values, df.values[:, 0])
    msg = r"\('foo', 'four'\)"
    with pytest.raises(KeyError, match=msg):
        df[('foo', 'four')]
    msg = "'foobar'"
    with pytest.raises(KeyError, match=msg):
        df['foobar']


@pytest.mark.filterwarnings("ignore:\\n.ix:DeprecationWarning")
def test_series_getitem(multiindex_year_month_day_dataframe_random_data):
    ymd = multiindex_year_month_day_dataframe_random_data
    s = ymd['A']

    result = s[2000, 3]

    # TODO(wesm): unused?
    # result2 = s.loc[2000, 3]

    expected = s.reindex(s.index[42:65])
    expected.index = expected.index.droplevel(0).droplevel(0)
    tm.assert_series_equal(result, expected)

    result = s[2000, 3, 10]
    expected = s[49]
    assert result == expected

    # fancy
    expected = s.reindex(s.index[49:51])
    result = s.loc[[(2000, 3, 10), (2000, 3, 13)]]
    tm.assert_series_equal(result, expected)

    result = s.ix[[(2000, 3, 10), (2000, 3, 13)]]
    tm.assert_series_equal(result, expected)

    # key error
    msg = "356"
    with pytest.raises(KeyError, match=msg):
        s.__getitem__((2000, 3, 4))


def test_series_getitem_corner(
        multiindex_year_month_day_dataframe_random_data):
    ymd = multiindex_year_month_day_dataframe_random_data
    s = ymd['A']

    # don't segfault, GH #495
    # out of bounds access
    msg = "index out of bounds"
    with pytest.raises(IndexError, match=msg):
        s.__getitem__(len(ymd))

    # generator
    result = s[(x > 0 for x in s)]
    expected = s[s > 0]
    tm.assert_series_equal(result, expected)


def test_frame_getitem_multicolumn_empty_level():
    f = DataFrame({'a': ['1', '2', '3'], 'b': ['2', '3', '4']})
    f.columns = [['level1 item1', 'level1 item2'], ['', 'level2 item2'],
                 ['level3 item1', 'level3 item2']]

    result = f['level1 item1']
    expected = DataFrame([['1'], ['2'], ['3']], index=f.index,
                         columns=['level3 item1'])
    tm.assert_frame_equal(result, expected)


@pytest.mark.filterwarnings("ignore:\\n.ix:DeprecationWarning")
def test_getitem_tuple_plus_slice():
    # GH #671
    df = DataFrame({'a': lrange(10),
                    'b': lrange(10),
                    'c': np.random.randn(10),
                    'd': np.random.randn(10)})

    idf = df.set_index(['a', 'b'])

    result = idf.loc[(0, 0), :]
    expected = idf.loc[0, 0]
    expected2 = idf.xs((0, 0))
    expected3 = idf.ix[0, 0]

    tm.assert_series_equal(result, expected)
    tm.assert_series_equal(result, expected2)
    tm.assert_series_equal(result, expected3)


def test_getitem_toplevel(multiindex_dataframe_random_data):
    frame = multiindex_dataframe_random_data
    df = frame.T

    result = df['foo']
    expected = df.reindex(columns=df.columns[:3])
    expected.columns = expected.columns.droplevel(0)
    tm.assert_frame_equal(result, expected)

    result = df['bar']
    result2 = df.loc[:, 'bar']

    expected = df.reindex(columns=df.columns[3:5])
    expected.columns = expected.columns.droplevel(0)
    tm.assert_frame_equal(result, expected)
    tm.assert_frame_equal(result, result2)


def test_getitem_int(multiindex_dataframe_random_data):
    levels = [[0, 1], [0, 1, 2]]
    codes = [[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]]
    index = MultiIndex(levels=levels, codes=codes)

    frame = DataFrame(np.random.randn(6, 2), index=index)

    result = frame.loc[1]
    expected = frame[-3:]
    expected.index = expected.index.droplevel(0)
    tm.assert_frame_equal(result, expected)

    # raises exception
    msg = "3"
    with pytest.raises(KeyError, match=msg):
        frame.loc.__getitem__(3)

    # however this will work
    frame = multiindex_dataframe_random_data
    result = frame.iloc[2]
    expected = frame.xs(frame.index[2])
    tm.assert_series_equal(result, expected)


def test_frame_getitem_view(multiindex_dataframe_random_data):
    frame = multiindex_dataframe_random_data
    df = frame.T.copy()

    # this works because we are modifying the underlying array
    # really a no-no
    df['foo'].values[:] = 0
    assert (df['foo'].values == 0).all()

    # but not if it's mixed-type
    df['foo', 'four'] = 'foo'
    df = df.sort_index(level=0, axis=1)

    # this will work, but will raise/warn as its chained assignment
    def f():
        df['foo']['one'] = 2
        return df

    msg = "A value is trying to be set on a copy of a slice from a DataFrame"
    with pytest.raises(com.SettingWithCopyError, match=msg):
        df['foo']['one'] = 2

    try:
        df = f()
    except ValueError:
        pass
    assert (df['foo', 'one'] == 0).all()


def test_getitem_lowerdim_corner(multiindex_dataframe_random_data):
    frame = multiindex_dataframe_random_data
    msg = "11"
    with pytest.raises(KeyError, match=msg):
        frame.loc.__getitem__((('bar', 'three'), 'B'))

    # in theory should be inserting in a sorted space????
    frame.loc[('bar', 'three'), 'B'] = 0
    assert frame.sort_index().loc[('bar', 'three'), 'B'] == 0


@pytest.mark.parametrize('unicode_strings', [True, False])
def test_mixed_depth_get(unicode_strings):
    # If unicode_strings is True, the column labels in dataframe
    # construction will use unicode strings in Python 2 (pull request
    # #17099).

    arrays = [['a', 'top', 'top', 'routine1', 'routine1', 'routine2'],
              ['', 'OD', 'OD', 'result1', 'result2', 'result1'],
              ['', 'wx', 'wy', '', '', '']]

    if unicode_strings:
        arrays = [[u(s) for s in arr] for arr in arrays]

    tuples = sorted(zip(*arrays))
    index = MultiIndex.from_tuples(tuples)
    df = DataFrame(np.random.randn(4, 6), columns=index)

    result = df['a']
    expected = df['a', '', ''].rename('a')
    tm.assert_series_equal(result, expected)

    result = df['routine1', 'result1']
    expected = df['routine1', 'result1', '']
    expected = expected.rename(('routine1', 'result1'))
    tm.assert_series_equal(result, expected)


def test_mi_access():

    # GH 4145
    data = """h1 main  h3 sub  h5
0  a    A   1  A1   1
1  b    B   2  B1   2
2  c    B   3  A1   3
3  d    A   4  B2   4
4  e    A   5  B2   5
5  f    B   6  A2   6
"""

    df = pd.read_csv(StringIO(data), sep=r'\s+', index_col=0)
    df2 = df.set_index(['main', 'sub']).T.sort_index(1)
    index = Index(['h1', 'h3', 'h5'])
    columns = MultiIndex.from_tuples([('A', 'A1')], names=['main', 'sub'])
    expected = DataFrame([['a', 1, 1]], index=columns, columns=index).T

    result = df2.loc[:, ('A', 'A1')]
    tm.assert_frame_equal(result, expected)

    result = df2[('A', 'A1')]
    tm.assert_frame_equal(result, expected)

    # GH 4146, not returning a block manager when selecting a unique index
    # from a duplicate index
    # as of 4879, this returns a Series (which is similar to what happens
    # with a non-unique)
    expected = Series(['a', 1, 1], index=['h1', 'h3', 'h5'], name='A1')
    result = df2['A']['A1']
    tm.assert_series_equal(result, expected)

    # selecting a non_unique from the 2nd level
    expected = DataFrame([['d', 4, 4], ['e', 5, 5]],
                         index=Index(['B2', 'B2'], name='sub'),
                         columns=['h1', 'h3', 'h5'], ).T
    result = df2['A']['B2']
    tm.assert_frame_equal(result, expected)
