import numpy as np
import pytest

from pandas import DataFrame, Index, MultiIndex, Series
from pandas.core.indexing import IndexingError
from pandas.util import testing as tm

# ----------------------------------------------------------------------------
# test indexing of Series with multi-level Index
# ----------------------------------------------------------------------------


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
def test_series_getitem_duplicates_multiindex(level0_value):
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
        with pytest.raises(KeyError, match=r"^'A'$"):
            df.val['A']

    with pytest.raises(KeyError, match=r"^'X'$"):
        df.val['X']

    result = df.val[level0_value]
    expected = Series(arr.ravel()[0:3], name='val', index=Index(
        [26, 37, 57], name='day'))
    tm.assert_series_equal(result, expected)


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


@pytest.mark.parametrize('indexer,expected_error,expected_error_msg', [
    (lambda s: s.__getitem__((2000, 3, 4)), KeyError, r"^356L?$"),
    (lambda s: s[(2000, 3, 4)], KeyError, r"^356L?$"),
    (lambda s: s.loc[(2000, 3, 4)], IndexingError, 'Too many indexers'),
    (lambda s: s.__getitem__(len(s)), IndexError, 'index out of bounds'),
    (lambda s: s[len(s)], IndexError, 'index out of bounds'),
    (lambda s: s.iloc[len(s)], IndexError,
     'single positional indexer is out-of-bounds')
])
def test_series_getitem_indexing_errors(
        multiindex_year_month_day_dataframe_random_data, indexer,
        expected_error, expected_error_msg):
    s = multiindex_year_month_day_dataframe_random_data['A']
    with pytest.raises(expected_error, match=expected_error_msg):
        indexer(s)


def test_series_getitem_corner_generator(
        multiindex_year_month_day_dataframe_random_data):
    s = multiindex_year_month_day_dataframe_random_data['A']
    result = s[(x > 0 for x in s)]
    expected = s[s > 0]
    tm.assert_series_equal(result, expected)


# ----------------------------------------------------------------------------
# test indexing of DataFrame with multi-level Index
# ----------------------------------------------------------------------------

def test_getitem_simple(multiindex_dataframe_random_data):
    df = multiindex_dataframe_random_data.T
    expected = df.values[:, 0]
    result = df['foo', 'one'].values
    tm.assert_almost_equal(result, expected)


@pytest.mark.parametrize('indexer,expected_error_msg', [
    (lambda df: df[('foo', 'four')], r"^\('foo', 'four'\)$"),
    (lambda df: df['foobar'], r"^'foobar'$")
])
def test_frame_getitem_simple_key_error(
        multiindex_dataframe_random_data, indexer, expected_error_msg):
    df = multiindex_dataframe_random_data.T
    with pytest.raises(KeyError, match=expected_error_msg):
        indexer(df)


def test_frame_getitem_multicolumn_empty_level():
    df = DataFrame({'a': ['1', '2', '3'], 'b': ['2', '3', '4']})
    df.columns = [['level1 item1', 'level1 item2'], ['', 'level2 item2'],
                  ['level3 item1', 'level3 item2']]

    result = df['level1 item1']
    expected = DataFrame([['1'], ['2'], ['3']], index=df.index,
                         columns=['level3 item1'])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('indexer,expected_slice', [
    (lambda df: df['foo'], slice(3)),
    (lambda df: df['bar'], slice(3, 5)),
    (lambda df: df.loc[:, 'bar'], slice(3, 5))
])
def test_frame_getitem_toplevel(
        multiindex_dataframe_random_data, indexer, expected_slice):
    df = multiindex_dataframe_random_data.T
    expected = df.reindex(columns=df.columns[expected_slice])
    expected.columns = expected.columns.droplevel(0)
    result = indexer(df)
    tm.assert_frame_equal(result, expected)


def test_frame_mixed_depth_get():
    arrays = [['a', 'top', 'top', 'routine1', 'routine1', 'routine2'],
              ['', 'OD', 'OD', 'result1', 'result2', 'result1'],
              ['', 'wx', 'wy', '', '', '']]

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


# ----------------------------------------------------------------------------
# test indexing of DataFrame with multi-level Index with duplicates
# ----------------------------------------------------------------------------

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


@pytest.mark.parametrize('indexer', [
    lambda df: df[('A', 'A1')],
    lambda df: df.loc[:, ('A', 'A1')]
])
def test_frame_mi_access(dataframe_with_duplicate_index, indexer):
    # GH 4145
    df = dataframe_with_duplicate_index
    index = Index(['h1', 'h3', 'h5'])
    columns = MultiIndex.from_tuples([('A', 'A1')], names=['main', 'sub'])
    expected = DataFrame([['a', 1, 1]], index=columns, columns=index).T

    result = indexer(df)
    tm.assert_frame_equal(result, expected)


def test_frame_mi_access_returns_series(dataframe_with_duplicate_index):
    # GH 4146, not returning a block manager when selecting a unique index
    # from a duplicate index
    # as of 4879, this returns a Series (which is similar to what happens
    # with a non-unique)
    df = dataframe_with_duplicate_index
    expected = Series(['a', 1, 1], index=['h1', 'h3', 'h5'], name='A1')
    result = df['A']['A1']
    tm.assert_series_equal(result, expected)


def test_frame_mi_access_returns_frame(dataframe_with_duplicate_index):
    # selecting a non_unique from the 2nd level
    df = dataframe_with_duplicate_index
    expected = DataFrame([['d', 4, 4], ['e', 5, 5]],
                         index=Index(['B2', 'B2'], name='sub'),
                         columns=['h1', 'h3', 'h5'], ).T
    result = df['A']['B2']
    tm.assert_frame_equal(result, expected)
