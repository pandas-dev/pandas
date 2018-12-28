import numpy as np
import pytest

from pandas.compat import range, u, zip

import pandas as pd
from pandas import DataFrame, Index, MultiIndex, Series
import pandas.core.common as com
from pandas.core.indexing import IndexingError
from pandas.util import testing as tm


@pytest.fixture
def frame_random_data_integer_multi_index():
    levels = [[0, 1], [0, 1, 2]]
    codes = [[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]]
    index = MultiIndex(levels=levels, codes=codes)
    return DataFrame(np.random.randn(6, 2), index=index)


@pytest.fixture
def dataframe_with_duplicate_index():
    """Fixture for DataFrame used in tests for gh-4145 and gh-4146"""
    data = [['a', 'd', 'e', 'c', 'f', 'b'],
            [1, 4, 5, 3, 6, 2],
            [1, 4, 5, 3, 6, 2]]
    index = ['h1', 'h3', 'h5']
    columns = MultiIndex(
        levels=[['A', 'B'], ['A1', 'A2', 'B1', 'B2']],
        codes=[[0, 0, 0, 1, 1, 1], [0, 3, 3, 0, 1, 2]],
        names=['main', 'sub'])
    return DataFrame(data, index=index, columns=columns)


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
    df = multiindex_dataframe_random_data.T
    expected = df.values[:, 0]
    result = df['foo', 'one'].values
    tm.assert_almost_equal(result, expected)


@pytest.mark.parametrize('indexer,msg', [
    (lambda df: df[('foo', 'four')], r"\('foo', 'four'\)"),
    (lambda df: df['foobar'], "'foobar'")
])
def test_getitem_simple_key_error(
        multiindex_dataframe_random_data, indexer, msg):
    df = multiindex_dataframe_random_data.T
    with pytest.raises(KeyError, match=msg):
        indexer(df)


@pytest.mark.parametrize('indexer', [
    lambda s: s[2000, 3],
    lambda s: s.loc[2000, 3]
])
def test_series_getitem(
        multiindex_year_month_day_dataframe_random_data, indexer):
    s = multiindex_year_month_day_dataframe_random_data['A']
    expected = s.reindex(s.index[42:65])
    expected.index = expected.index.droplevel(0).droplevel(0)

    result = indexer(s)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize('indexer', [
    lambda s: s[2000, 3, 10],
    lambda s: s.loc[2000, 3, 10]
])
def test_series_getitem_returns_scalar(
        multiindex_year_month_day_dataframe_random_data, indexer):
    s = multiindex_year_month_day_dataframe_random_data['A']
    expected = s.iloc[49]

    result = indexer(s)
    assert result == expected


@pytest.mark.filterwarnings("ignore:\\n.ix:DeprecationWarning")
@pytest.mark.parametrize('indexer', [
    lambda s: s.loc[[(2000, 3, 10), (2000, 3, 13)]],
    lambda s: s.ix[[(2000, 3, 10), (2000, 3, 13)]]
])
def test_series_getitem_fancy(
        multiindex_year_month_day_dataframe_random_data, indexer):
    s = multiindex_year_month_day_dataframe_random_data['A']
    expected = s.reindex(s.index[49:51])

    result = indexer(s)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize('indexer,error,msg', [
    (lambda s: s.__getitem__((2000, 3, 4)), KeyError, '356'),
    (lambda s: s[(2000, 3, 4)], KeyError, '356'),
    (lambda s: s.loc[(2000, 3, 4)], IndexingError, 'Too many indexers'),
    (lambda s: s.__getitem__(len(s)), IndexError, 'index out of bounds'),
    (lambda s: s[len(s)], IndexError, 'index out of bounds'),
    (lambda s: s.iloc[len(s)], IndexError,
     'single positional indexer is out-of-bounds')
])
def test_series_getitem_indexing_errors(
        multiindex_year_month_day_dataframe_random_data, indexer, error, msg):
    s = multiindex_year_month_day_dataframe_random_data['A']
    with pytest.raises(error, match=msg):
        indexer(s)


def test_series_getitem_corner_generator(
        multiindex_year_month_day_dataframe_random_data):
    s = multiindex_year_month_day_dataframe_random_data['A']
    result = s[(x > 0 for x in s)]
    expected = s[s > 0]
    tm.assert_series_equal(result, expected)


def test_frame_getitem_multicolumn_empty_level():
    df = DataFrame({'a': ['1', '2', '3'], 'b': ['2', '3', '4']})
    df.columns = [['level1 item1', 'level1 item2'], ['', 'level2 item2'],
                  ['level3 item1', 'level3 item2']]

    result = df['level1 item1']
    expected = DataFrame([['1'], ['2'], ['3']], index=df.index,
                         columns=['level3 item1'])
    tm.assert_frame_equal(result, expected)


def test_getitem_tuple_plus_slice():
    # GH 671
    df = DataFrame({'a': np.arange(10),
                    'b': np.arange(10),
                    'c': np.random.randn(10),
                    'd': np.random.randn(10)}
                   ).set_index(['a', 'b'])
    expected = df.loc[0, 0]
    result = df.loc[(0, 0), :]
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize('indexer,expected_slice', [
    (lambda df: df['foo'], slice(3)),
    (lambda df: df['bar'], slice(3, 5)),
    (lambda df: df.loc[:, 'bar'], slice(3, 5))
])
def test_getitem_toplevel(
        multiindex_dataframe_random_data, indexer, expected_slice):
    df = multiindex_dataframe_random_data.T
    expected = df.reindex(columns=df.columns[expected_slice])
    expected.columns = expected.columns.droplevel(0)
    result = indexer(df)
    tm.assert_frame_equal(result, expected)


def test_getitem_int(frame_random_data_integer_multi_index):
    df = frame_random_data_integer_multi_index
    result = df.loc[1]
    expected = df[-3:]
    expected.index = expected.index.droplevel(0)
    tm.assert_frame_equal(result, expected)


def test_getitem_int_raises_exception(frame_random_data_integer_multi_index):
    df = frame_random_data_integer_multi_index
    msg = "3"
    with pytest.raises(KeyError, match=msg):
        df.loc.__getitem__(3)


def test_getitem_iloc(multiindex_dataframe_random_data):
    df = multiindex_dataframe_random_data
    result = df.iloc[2]
    expected = df.xs(df.index[2])
    tm.assert_series_equal(result, expected)


def test_frame_setitem_view_direct(multiindex_dataframe_random_data):
    # this works because we are modifying the underlying array
    # really a no-no
    df = multiindex_dataframe_random_data.T
    df['foo'].values[:] = 0
    assert (df['foo'].values == 0).all()


def test_frame_setitem_copy_raises(multiindex_dataframe_random_data):
    # will raise/warn as its chained assignment
    df = multiindex_dataframe_random_data.T
    msg = "A value is trying to be set on a copy of a slice from a DataFrame"
    with pytest.raises(com.SettingWithCopyError, match=msg):
        df['foo']['one'] = 2


def test_frame_setitem_copy_no_write(multiindex_dataframe_random_data):
    frame = multiindex_dataframe_random_data.T
    expected = frame
    df = frame.copy()
    msg = "A value is trying to be set on a copy of a slice from a DataFrame"
    with pytest.raises(com.SettingWithCopyError, match=msg):
        df['foo']['one'] = 2

    result = df
    tm.assert_frame_equal(result, expected)


def test_getitem_lowerdim_corner(multiindex_dataframe_random_data):
    df = multiindex_dataframe_random_data

    # test setup - check key not in dataframe
    with pytest.raises(KeyError, match="11"):
        df.loc[('bar', 'three'), 'B']

    # in theory should be inserting in a sorted space????
    df.loc[('bar', 'three'), 'B'] = 0
    expected = 0
    result = df.sort_index().loc[('bar', 'three'), 'B']
    assert result == expected


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


@pytest.mark.parametrize('indexer', [
    lambda df: df.loc[:, ('A', 'A1')],
    lambda df: df[('A', 'A1')]
])
def test_mi_access(dataframe_with_duplicate_index, indexer):
    # GH 4145
    df = dataframe_with_duplicate_index
    index = Index(['h1', 'h3', 'h5'])
    columns = MultiIndex.from_tuples([('A', 'A1')], names=['main', 'sub'])
    expected = DataFrame([['a', 1, 1]], index=columns, columns=index).T

    result = indexer(df)
    tm.assert_frame_equal(result, expected)


def test_mi_access_returns_series(dataframe_with_duplicate_index):
    # GH 4146, not returning a block manager when selecting a unique index
    # from a duplicate index
    # as of 4879, this returns a Series (which is similar to what happens
    # with a non-unique)
    df = dataframe_with_duplicate_index
    expected = Series(['a', 1, 1], index=['h1', 'h3', 'h5'], name='A1')
    result = df['A']['A1']
    tm.assert_series_equal(result, expected)


def test_mi_access_returns_frame(dataframe_with_duplicate_index):
    # selecting a non_unique from the 2nd level
    df = dataframe_with_duplicate_index
    expected = DataFrame([['d', 4, 4], ['e', 5, 5]],
                         index=Index(['B2', 'B2'], name='sub'),
                         columns=['h1', 'h3', 'h5'], ).T
    result = df['A']['B2']
    tm.assert_frame_equal(result, expected)
