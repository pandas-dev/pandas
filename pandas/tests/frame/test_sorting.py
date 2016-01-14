# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np

from pandas.compat import lrange
from pandas import (DataFrame, Series, MultiIndex, Timestamp,
                    date_range)

from pandas.util.testing import (assert_series_equal,
                                 assert_frame_equal,
                                 assertRaisesRegexp)

import pandas.util.testing as tm

from pandas.tests.frame.common import TestData


class TestDataFrameSorting(tm.TestCase, TestData):

    _multiprocess_can_split_ = True

    def test_sort_values(self):
        # API for 9816

        # sort_index
        frame = DataFrame(np.arange(16).reshape(4, 4), index=[1, 2, 3, 4],
                          columns=['A', 'B', 'C', 'D'])

        # 9816 deprecated
        with tm.assert_produces_warning(FutureWarning):
            frame.sort(columns='A')
        with tm.assert_produces_warning(FutureWarning):
            frame.sort()

        unordered = frame.ix[[3, 2, 4, 1]]
        expected = unordered.sort_index()

        result = unordered.sort_index(axis=0)
        assert_frame_equal(result, expected)

        unordered = frame.ix[:, [2, 1, 3, 0]]
        expected = unordered.sort_index(axis=1)

        result = unordered.sort_index(axis=1)
        assert_frame_equal(result, expected)
        assert_frame_equal(result, expected)

        # sortlevel
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        df = DataFrame([[1, 2], [3, 4]], mi)

        result = df.sort_index(level='A', sort_remaining=False)
        expected = df.sortlevel('A', sort_remaining=False)
        assert_frame_equal(result, expected)

        df = df.T
        result = df.sort_index(level='A', axis=1, sort_remaining=False)
        expected = df.sortlevel('A', axis=1, sort_remaining=False)
        assert_frame_equal(result, expected)

        # MI sort, but no by
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        df = DataFrame([[1, 2], [3, 4]], mi)
        result = df.sort_index(sort_remaining=False)
        expected = df.sort_index()
        assert_frame_equal(result, expected)

    def test_sort_index(self):
        frame = DataFrame(np.arange(16).reshape(4, 4), index=[1, 2, 3, 4],
                          columns=['A', 'B', 'C', 'D'])

        # axis=0
        unordered = frame.ix[[3, 2, 4, 1]]
        sorted_df = unordered.sort_index(axis=0)
        expected = frame
        assert_frame_equal(sorted_df, expected)

        sorted_df = unordered.sort_index(ascending=False)
        expected = frame[::-1]
        assert_frame_equal(sorted_df, expected)

        # axis=1
        unordered = frame.ix[:, ['D', 'B', 'C', 'A']]
        sorted_df = unordered.sort_index(axis=1)
        expected = frame
        assert_frame_equal(sorted_df, expected)

        sorted_df = unordered.sort_index(axis=1, ascending=False)
        expected = frame.ix[:, ::-1]
        assert_frame_equal(sorted_df, expected)

        # by column
        sorted_df = frame.sort_values(by='A')
        indexer = frame['A'].argsort().values
        expected = frame.ix[frame.index[indexer]]
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.sort_values(by='A', ascending=False)
        indexer = indexer[::-1]
        expected = frame.ix[frame.index[indexer]]
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.sort_values(by='A', ascending=False)
        assert_frame_equal(sorted_df, expected)

        # GH4839
        sorted_df = frame.sort_values(by=['A'], ascending=[False])
        assert_frame_equal(sorted_df, expected)

        # check for now
        sorted_df = frame.sort_values(by='A')
        assert_frame_equal(sorted_df, expected[::-1])
        expected = frame.sort_values(by='A')
        assert_frame_equal(sorted_df, expected)

        expected = frame.sort_values(by=['A', 'B'], ascending=False)
        sorted_df = frame.sort_values(by=['A', 'B'])
        assert_frame_equal(sorted_df, expected[::-1])

        self.assertRaises(ValueError, lambda: frame.sort_values(
            by=['A', 'B'], axis=2, inplace=True))

        msg = 'When sorting by column, axis must be 0'
        with assertRaisesRegexp(ValueError, msg):
            frame.sort_values(by='A', axis=1)

        msg = r'Length of ascending \(5\) != length of by \(2\)'
        with assertRaisesRegexp(ValueError, msg):
            frame.sort_values(by=['A', 'B'], axis=0, ascending=[True] * 5)

    def test_sort_index_categorical_index(self):

        df = (DataFrame({'A': np.arange(6, dtype='int64'),
                         'B': Series(list('aabbca'))
                         .astype('category', categories=list('cab'))})
              .set_index('B'))

        result = df.sort_index()
        expected = df.iloc[[4, 0, 1, 5, 2, 3]]
        assert_frame_equal(result, expected)

        result = df.sort_index(ascending=False)
        expected = df.iloc[[3, 2, 5, 1, 0, 4]]
        assert_frame_equal(result, expected)

    def test_sort_nan(self):
        # GH3917
        nan = np.nan
        df = DataFrame({'A': [1, 2, nan, 1, 6, 8, 4],
                        'B': [9, nan, 5, 2, 5, 4, 5]})

        # sort one column only
        expected = DataFrame(
            {'A': [nan, 1, 1, 2, 4, 6, 8],
             'B': [5, 9, 2, nan, 5, 5, 4]},
            index=[2, 0, 3, 1, 6, 4, 5])
        sorted_df = df.sort_values(['A'], na_position='first')
        assert_frame_equal(sorted_df, expected)

        expected = DataFrame(
            {'A': [nan, 8, 6, 4, 2, 1, 1],
             'B': [5, 4, 5, 5, nan, 9, 2]},
            index=[2, 5, 4, 6, 1, 0, 3])
        sorted_df = df.sort_values(['A'], na_position='first', ascending=False)
        assert_frame_equal(sorted_df, expected)

        # na_position='last', order
        expected = DataFrame(
            {'A': [1, 1, 2, 4, 6, 8, nan],
             'B': [2, 9, nan, 5, 5, 4, 5]},
            index=[3, 0, 1, 6, 4, 5, 2])
        sorted_df = df.sort_values(['A', 'B'])
        assert_frame_equal(sorted_df, expected)

        # na_position='first', order
        expected = DataFrame(
            {'A': [nan, 1, 1, 2, 4, 6, 8],
             'B': [5, 2, 9, nan, 5, 5, 4]},
            index=[2, 3, 0, 1, 6, 4, 5])
        sorted_df = df.sort_values(['A', 'B'], na_position='first')
        assert_frame_equal(sorted_df, expected)

        # na_position='first', not order
        expected = DataFrame(
            {'A': [nan, 1, 1, 2, 4, 6, 8],
             'B': [5, 9, 2, nan, 5, 5, 4]},
            index=[2, 0, 3, 1, 6, 4, 5])
        sorted_df = df.sort_values(['A', 'B'], ascending=[
                                   1, 0], na_position='first')
        assert_frame_equal(sorted_df, expected)

        # na_position='last', not order
        expected = DataFrame(
            {'A': [8, 6, 4, 2, 1, 1, nan],
             'B': [4, 5, 5, nan, 2, 9, 5]},
            index=[5, 4, 6, 1, 3, 0, 2])
        sorted_df = df.sort_values(['A', 'B'], ascending=[
                                   0, 1], na_position='last')
        assert_frame_equal(sorted_df, expected)

        # Test DataFrame with nan label
        df = DataFrame({'A': [1, 2, nan, 1, 6, 8, 4],
                        'B': [9, nan, 5, 2, 5, 4, 5]},
                       index=[1, 2, 3, 4, 5, 6, nan])

        # NaN label, ascending=True, na_position='last'
        sorted_df = df.sort_index(
            kind='quicksort', ascending=True, na_position='last')
        expected = DataFrame({'A': [1, 2, nan, 1, 6, 8, 4],
                              'B': [9, nan, 5, 2, 5, 4, 5]},
                             index=[1, 2, 3, 4, 5, 6, nan])
        assert_frame_equal(sorted_df, expected)

        # NaN label, ascending=True, na_position='first'
        sorted_df = df.sort_index(na_position='first')
        expected = DataFrame({'A': [4, 1, 2, nan, 1, 6, 8],
                              'B': [5, 9, nan, 5, 2, 5, 4]},
                             index=[nan, 1, 2, 3, 4, 5, 6])
        assert_frame_equal(sorted_df, expected)

        # NaN label, ascending=False, na_position='last'
        sorted_df = df.sort_index(kind='quicksort', ascending=False)
        expected = DataFrame({'A': [8, 6, 1, nan, 2, 1, 4],
                              'B': [4, 5, 2, 5, nan, 9, 5]},
                             index=[6, 5, 4, 3, 2, 1, nan])
        assert_frame_equal(sorted_df, expected)

        # NaN label, ascending=False, na_position='first'
        sorted_df = df.sort_index(
            kind='quicksort', ascending=False, na_position='first')
        expected = DataFrame({'A': [4, 8, 6, 1, nan, 2, 1],
                              'B': [5, 4, 5, 2, 5, nan, 9]},
                             index=[nan, 6, 5, 4, 3, 2, 1])
        assert_frame_equal(sorted_df, expected)

    def test_stable_descending_sort(self):
        # GH #6399
        df = DataFrame([[2, 'first'], [2, 'second'], [1, 'a'], [1, 'b']],
                       columns=['sort_col', 'order'])
        sorted_df = df.sort_values(by='sort_col', kind='mergesort',
                                   ascending=False)
        assert_frame_equal(df, sorted_df)

    def test_stable_descending_multicolumn_sort(self):
        nan = np.nan
        df = DataFrame({'A': [1, 2, nan, 1, 6, 8, 4],
                        'B': [9, nan, 5, 2, 5, 4, 5]})
        # test stable mergesort
        expected = DataFrame(
            {'A': [nan, 8, 6, 4, 2, 1, 1],
             'B': [5, 4, 5, 5, nan, 2, 9]},
            index=[2, 5, 4, 6, 1, 3, 0])
        sorted_df = df.sort_values(['A', 'B'], ascending=[0, 1],
                                   na_position='first',
                                   kind='mergesort')
        assert_frame_equal(sorted_df, expected)

        expected = DataFrame(
            {'A': [nan, 8, 6, 4, 2, 1, 1],
             'B': [5, 4, 5, 5, nan, 9, 2]},
            index=[2, 5, 4, 6, 1, 0, 3])
        sorted_df = df.sort_values(['A', 'B'], ascending=[0, 0],
                                   na_position='first',
                                   kind='mergesort')
        assert_frame_equal(sorted_df, expected)

    def test_sort_index_multicolumn(self):
        import random
        A = np.arange(5).repeat(20)
        B = np.tile(np.arange(5), 20)
        random.shuffle(A)
        random.shuffle(B)
        frame = DataFrame({'A': A, 'B': B,
                           'C': np.random.randn(100)})

        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            frame.sort_index(by=['A', 'B'])
        result = frame.sort_values(by=['A', 'B'])
        indexer = np.lexsort((frame['B'], frame['A']))
        expected = frame.take(indexer)
        assert_frame_equal(result, expected)

        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            frame.sort_index(by=['A', 'B'], ascending=False)
        result = frame.sort_values(by=['A', 'B'], ascending=False)
        indexer = np.lexsort((frame['B'].rank(ascending=False),
                              frame['A'].rank(ascending=False)))
        expected = frame.take(indexer)
        assert_frame_equal(result, expected)

        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            frame.sort_index(by=['B', 'A'])
        result = frame.sort_values(by=['B', 'A'])
        indexer = np.lexsort((frame['A'], frame['B']))
        expected = frame.take(indexer)
        assert_frame_equal(result, expected)

    def test_sort_index_inplace(self):
        frame = DataFrame(np.random.randn(4, 4), index=[1, 2, 3, 4],
                          columns=['A', 'B', 'C', 'D'])

        # axis=0
        unordered = frame.ix[[3, 2, 4, 1]]
        a_id = id(unordered['A'])
        df = unordered.copy()
        df.sort_index(inplace=True)
        expected = frame
        assert_frame_equal(df, expected)
        self.assertNotEqual(a_id, id(df['A']))

        df = unordered.copy()
        df.sort_index(ascending=False, inplace=True)
        expected = frame[::-1]
        assert_frame_equal(df, expected)

        # axis=1
        unordered = frame.ix[:, ['D', 'B', 'C', 'A']]
        df = unordered.copy()
        df.sort_index(axis=1, inplace=True)
        expected = frame
        assert_frame_equal(df, expected)

        df = unordered.copy()
        df.sort_index(axis=1, ascending=False, inplace=True)
        expected = frame.ix[:, ::-1]
        assert_frame_equal(df, expected)

    def test_sort_index_different_sortorder(self):
        A = np.arange(20).repeat(5)
        B = np.tile(np.arange(5), 20)

        indexer = np.random.permutation(100)
        A = A.take(indexer)
        B = B.take(indexer)

        df = DataFrame({'A': A, 'B': B,
                        'C': np.random.randn(100)})

        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            df.sort_index(by=['A', 'B'], ascending=[1, 0])
        result = df.sort_values(by=['A', 'B'], ascending=[1, 0])

        ex_indexer = np.lexsort((df.B.max() - df.B, df.A))
        expected = df.take(ex_indexer)
        assert_frame_equal(result, expected)

        # test with multiindex, too
        idf = df.set_index(['A', 'B'])

        result = idf.sort_index(ascending=[1, 0])
        expected = idf.take(ex_indexer)
        assert_frame_equal(result, expected)

        # also, Series!
        result = idf['C'].sort_index(ascending=[1, 0])
        assert_series_equal(result, expected['C'])

    def test_sort_inplace(self):
        frame = DataFrame(np.random.randn(4, 4), index=[1, 2, 3, 4],
                          columns=['A', 'B', 'C', 'D'])

        sorted_df = frame.copy()
        sorted_df.sort_values(by='A', inplace=True)
        expected = frame.sort_values(by='A')
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.copy()
        sorted_df.sort_values(by='A', ascending=False, inplace=True)
        expected = frame.sort_values(by='A', ascending=False)
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.copy()
        sorted_df.sort_values(by=['A', 'B'], ascending=False, inplace=True)
        expected = frame.sort_values(by=['A', 'B'], ascending=False)
        assert_frame_equal(sorted_df, expected)

    def test_sort_index_duplicates(self):

        # with 9816, these are all translated to .sort_values

        df = DataFrame([lrange(5, 9), lrange(4)],
                       columns=['a', 'a', 'b', 'b'])

        with assertRaisesRegexp(ValueError, 'duplicate'):
            # use .sort_values #9816
            with tm.assert_produces_warning(FutureWarning):
                df.sort_index(by='a')
        with assertRaisesRegexp(ValueError, 'duplicate'):
            df.sort_values(by='a')

        with assertRaisesRegexp(ValueError, 'duplicate'):
            # use .sort_values #9816
            with tm.assert_produces_warning(FutureWarning):
                df.sort_index(by=['a'])
        with assertRaisesRegexp(ValueError, 'duplicate'):
            df.sort_values(by=['a'])

        with assertRaisesRegexp(ValueError, 'duplicate'):
            # use .sort_values #9816
            with tm.assert_produces_warning(FutureWarning):
                # multi-column 'by' is separate codepath
                df.sort_index(by=['a', 'b'])
        with assertRaisesRegexp(ValueError, 'duplicate'):
            # multi-column 'by' is separate codepath
            df.sort_values(by=['a', 'b'])

        # with multi-index
        # GH4370
        df = DataFrame(np.random.randn(4, 2),
                       columns=MultiIndex.from_tuples([('a', 0), ('a', 1)]))
        with assertRaisesRegexp(ValueError, 'levels'):
            # use .sort_values #9816
            with tm.assert_produces_warning(FutureWarning):
                df.sort_index(by='a')
        with assertRaisesRegexp(ValueError, 'levels'):
            df.sort_values(by='a')

        # convert tuples to a list of tuples
        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            df.sort_index(by=[('a', 1)])
        expected = df.sort_values(by=[('a', 1)])

        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            df.sort_index(by=('a', 1))
        result = df.sort_values(by=('a', 1))
        assert_frame_equal(result, expected)

    def test_sortlevel(self):
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        df = DataFrame([[1, 2], [3, 4]], mi)
        res = df.sortlevel('A', sort_remaining=False)
        assert_frame_equal(df, res)

        res = df.sortlevel(['A', 'B'], sort_remaining=False)
        assert_frame_equal(df, res)

    def test_sort_datetimes(self):

        # GH 3461, argsort / lexsort differences for a datetime column
        df = DataFrame(['a', 'a', 'a', 'b', 'c', 'd', 'e', 'f', 'g'],
                       columns=['A'],
                       index=date_range('20130101', periods=9))
        dts = [Timestamp(x)
               for x in ['2004-02-11', '2004-01-21', '2004-01-26',
                         '2005-09-20', '2010-10-04', '2009-05-12',
                         '2008-11-12', '2010-09-28', '2010-09-28']]
        df['B'] = dts[::2] + dts[1::2]
        df['C'] = 2.
        df['A1'] = 3.

        df1 = df.sort_values(by='A')
        df2 = df.sort_values(by=['A'])
        assert_frame_equal(df1, df2)

        df1 = df.sort_values(by='B')
        df2 = df.sort_values(by=['B'])
        assert_frame_equal(df1, df2)

    def test_frame_column_inplace_sort_exception(self):
        s = self.frame['A']
        with assertRaisesRegexp(ValueError, "This Series is a view"):
            s.sort_values(inplace=True)

        cp = s.copy()
        cp.sort_values()  # it works!
