from warnings import catch_warnings, simplefilter

import numpy as np
import pytest

from pandas.compat import lrange, range, u, zip

import pandas as pd
from pandas import DataFrame, Index, MultiIndex, Series, date_range
import pandas.core.common as com
from pandas.util import testing as tm


@pytest.mark.filterwarnings("ignore:\\n.ix:DeprecationWarning")
class TestMultiIndexGetItem(object):

    def test_series_getitem_multiindex(self):

        # GH 6018
        # series regression getitem with a multi-index

        s = Series([1, 2, 3])
        s.index = MultiIndex.from_tuples([(0, 0), (1, 1), (2, 1)])

        result = s[:, 0]
        expected = Series([1], index=[0])
        tm.assert_series_equal(result, expected)

        result = s.loc[:, 1]
        expected = Series([2, 3], index=[1, 2])
        tm.assert_series_equal(result, expected)

        # xs
        result = s.xs(0, level=0)
        expected = Series([1], index=[0])
        tm.assert_series_equal(result, expected)

        result = s.xs(1, level=1)
        expected = Series([2, 3], index=[1, 2])
        tm.assert_series_equal(result, expected)

        # GH6258
        dt = list(date_range('20130903', periods=3))
        idx = MultiIndex.from_product([list('AB'), dt])
        s = Series([1, 3, 4, 1, 3, 4], index=idx)

        result = s.xs('20130903', level=1)
        expected = Series([1, 1], index=list('AB'))
        tm.assert_series_equal(result, expected)

        # GH5684
        idx = MultiIndex.from_tuples([('a', 'one'), ('a', 'two'), ('b', 'one'),
                                      ('b', 'two')])
        s = Series([1, 2, 3, 4], index=idx)
        s.index.set_names(['L1', 'L2'], inplace=True)
        result = s.xs('one', level='L2')
        expected = Series([1, 3], index=['a', 'b'])
        expected.index.set_names(['L1'], inplace=True)
        tm.assert_series_equal(result, expected)

    def test_getitem_duplicates_multiindex(self):
        # GH 5725 the 'A' happens to be a valid Timestamp so the doesn't raise
        # the appropriate error, only in PY3 of course!

        index = MultiIndex(levels=[['D', 'B', 'C'],
                                   [0, 26, 27, 37, 57, 67, 75, 82]],
                           labels=[[0, 0, 0, 1, 2, 2, 2, 2, 2, 2],
                                   [1, 3, 4, 6, 0, 2, 2, 3, 5, 7]],
                           names=['tag', 'day'])
        arr = np.random.randn(len(index), 1)
        df = DataFrame(arr, index=index, columns=['val'])
        result = df.val['D']
        expected = Series(arr.ravel()[0:3], name='val', index=Index(
            [26, 37, 57], name='day'))
        tm.assert_series_equal(result, expected)

        def f():
            df.val['A']

        pytest.raises(KeyError, f)

        def f():
            df.val['X']

        pytest.raises(KeyError, f)

        # A is treated as a special Timestamp
        index = MultiIndex(levels=[['A', 'B', 'C'],
                                   [0, 26, 27, 37, 57, 67, 75, 82]],
                           labels=[[0, 0, 0, 1, 2, 2, 2, 2, 2, 2],
                                   [1, 3, 4, 6, 0, 2, 2, 3, 5, 7]],
                           names=['tag', 'day'])
        df = DataFrame(arr, index=index, columns=['val'])
        result = df.val['A']
        expected = Series(arr.ravel()[0:3], name='val', index=Index(
            [26, 37, 57], name='day'))
        tm.assert_series_equal(result, expected)

        def f():
            df.val['X']

        pytest.raises(KeyError, f)

        # GH 7866
        # multi-index slicing with missing indexers
        idx = MultiIndex.from_product([['A', 'B', 'C'],
                                       ['foo', 'bar', 'baz']],
                                      names=['one', 'two'])
        s = Series(np.arange(9, dtype='int64'), index=idx).sort_index()

        exp_idx = MultiIndex.from_product([['A'], ['foo', 'bar', 'baz']],
                                          names=['one', 'two'])
        expected = Series(np.arange(3, dtype='int64'),
                          index=exp_idx).sort_index()

        result = s.loc[['A']]
        tm.assert_series_equal(result, expected)
        result = s.loc[['A', 'D']]
        tm.assert_series_equal(result, expected)

        # not any values found
        pytest.raises(KeyError, lambda: s.loc[['D']])

        # empty ok
        result = s.loc[[]]
        expected = s.iloc[[]]
        tm.assert_series_equal(result, expected)

        idx = pd.IndexSlice
        expected = Series([0, 3, 6], index=MultiIndex.from_product(
            [['A', 'B', 'C'], ['foo']], names=['one', 'two'])).sort_index()

        result = s.loc[idx[:, ['foo']]]
        tm.assert_series_equal(result, expected)
        result = s.loc[idx[:, ['foo', 'bah']]]
        tm.assert_series_equal(result, expected)

        # GH 8737
        # empty indexer
        multi_index = MultiIndex.from_product((['foo', 'bar', 'baz'],
                                               ['alpha', 'beta']))
        df = DataFrame(
            np.random.randn(5, 6), index=range(5), columns=multi_index)
        df = df.sort_index(level=0, axis=1)

        expected = DataFrame(index=range(5),
                             columns=multi_index.reindex([])[0])
        result1 = df.loc[:, ([], slice(None))]
        result2 = df.loc[:, (['foo'], [])]
        tm.assert_frame_equal(result1, expected)
        tm.assert_frame_equal(result2, expected)

        # regression from < 0.14.0
        # GH 7914
        df = DataFrame([[np.mean, np.median], ['mean', 'median']],
                       columns=MultiIndex.from_tuples([('functs', 'mean'),
                                                       ('functs', 'median')]),
                       index=['function', 'name'])
        result = df.loc['function', ('functs', 'mean')]
        assert result == np.mean

    def test_getitem_simple(self, multiindex_dataframe_random_data):
        frame = multiindex_dataframe_random_data
        df = frame.T

        col = df['foo', 'one']
        tm.assert_almost_equal(col.values, df.values[:, 0])
        with pytest.raises(KeyError):
            df[('foo', 'four')]
        with pytest.raises(KeyError):
            df['foobar']

    def test_series_getitem(
            self, multiindex_year_month_day_dataframe_random_data):
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

        with catch_warnings(record=True):
            simplefilter("ignore", DeprecationWarning)
            result = s.ix[[(2000, 3, 10), (2000, 3, 13)]]
        tm.assert_series_equal(result, expected)

        # key error
        pytest.raises(KeyError, s.__getitem__, (2000, 3, 4))

    def test_series_getitem_corner(
            self, multiindex_year_month_day_dataframe_random_data):
        ymd = multiindex_year_month_day_dataframe_random_data
        s = ymd['A']

        # don't segfault, GH #495
        # out of bounds access
        pytest.raises(IndexError, s.__getitem__, len(ymd))

        # generator
        result = s[(x > 0 for x in s)]
        expected = s[s > 0]
        tm.assert_series_equal(result, expected)

    def test_frame_getitem_multicolumn_empty_level(self):
        f = DataFrame({'a': ['1', '2', '3'], 'b': ['2', '3', '4']})
        f.columns = [['level1 item1', 'level1 item2'], ['', 'level2 item2'],
                     ['level3 item1', 'level3 item2']]

        result = f['level1 item1']
        expected = DataFrame([['1'], ['2'], ['3']], index=f.index,
                             columns=['level3 item1'])
        tm.assert_frame_equal(result, expected)

    def test_getitem_tuple_plus_slice(self):
        # GH #671
        df = DataFrame({'a': lrange(10),
                        'b': lrange(10),
                        'c': np.random.randn(10),
                        'd': np.random.randn(10)})

        idf = df.set_index(['a', 'b'])

        result = idf.loc[(0, 0), :]
        expected = idf.loc[0, 0]
        expected2 = idf.xs((0, 0))
        with catch_warnings(record=True):
            simplefilter("ignore", DeprecationWarning)
            expected3 = idf.ix[0, 0]

        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result, expected2)
        tm.assert_series_equal(result, expected3)

    def test_getitem_toplevel(self, multiindex_dataframe_random_data):
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

    def test_getitem_int(self, multiindex_dataframe_random_data):
        levels = [[0, 1], [0, 1, 2]]
        labels = [[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]]
        index = MultiIndex(levels=levels, labels=labels)

        frame = DataFrame(np.random.randn(6, 2), index=index)

        result = frame.loc[1]
        expected = frame[-3:]
        expected.index = expected.index.droplevel(0)
        tm.assert_frame_equal(result, expected)

        # raises exception
        pytest.raises(KeyError, frame.loc.__getitem__, 3)

        # however this will work
        frame = multiindex_dataframe_random_data
        result = frame.iloc[2]
        expected = frame.xs(frame.index[2])
        tm.assert_series_equal(result, expected)

    def test_frame_getitem_view(self, multiindex_dataframe_random_data):
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

        pytest.raises(com.SettingWithCopyError, f)

        try:
            df = f()
        except ValueError:
            pass
        assert (df['foo', 'one'] == 0).all()

    def test_getitem_lowerdim_corner(self, multiindex_dataframe_random_data):
        frame = multiindex_dataframe_random_data
        pytest.raises(KeyError, frame.loc.__getitem__,
                      (('bar', 'three'), 'B'))

        # in theory should be inserting in a sorted space????
        frame.loc[('bar', 'three'), 'B'] = 0
        assert frame.sort_index().loc[('bar', 'three'), 'B'] == 0

    @pytest.mark.parametrize('unicode_strings', [True, False])
    def test_mixed_depth_get(self, unicode_strings):
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
