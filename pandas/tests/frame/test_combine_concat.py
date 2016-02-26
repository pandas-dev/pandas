# -*- coding: utf-8 -*-

from __future__ import print_function

from datetime import datetime

from numpy import nan
import numpy as np

from pandas.compat import lrange
from pandas import DataFrame, Series, Index, Timestamp
import pandas as pd

from pandas.util.testing import (assert_series_equal,
                                 assert_frame_equal,
                                 assertRaisesRegexp)

import pandas.util.testing as tm

from pandas.tests.frame.common import TestData


class TestDataFrameCombineConcat(tm.TestCase, TestData):

    _multiprocess_can_split_ = True

    def test_combine_first_mixed(self):
        a = Series(['a', 'b'], index=lrange(2))
        b = Series(lrange(2), index=lrange(2))
        f = DataFrame({'A': a, 'B': b})

        a = Series(['a', 'b'], index=lrange(5, 7))
        b = Series(lrange(2), index=lrange(5, 7))
        g = DataFrame({'A': a, 'B': b})

        # TODO(wesm): no verification?
        combined = f.combine_first(g)  # noqa

    def test_combine_multiple_frames_dtypes(self):

        # GH 2759
        A = DataFrame(data=np.ones((10, 2)), columns=[
                      'foo', 'bar'], dtype=np.float64)
        B = DataFrame(data=np.ones((10, 2)), dtype=np.float32)
        results = pd.concat((A, B), axis=1).get_dtype_counts()
        expected = Series(dict(float64=2, float32=2))
        assert_series_equal(results, expected)

    def test_combine_multiple_tzs(self):
        # GH 12467
        # combining datetime tz-aware and naive DataFrames
        ts1 = Timestamp('2015-01-01', tz=None)
        ts2 = Timestamp('2015-01-01', tz='UTC')
        ts3 = Timestamp('2015-01-01', tz='EST')

        df1 = DataFrame(dict(time=[ts1]))
        df2 = DataFrame(dict(time=[ts2]))
        df3 = DataFrame(dict(time=[ts3]))

        results = pd.concat([df1, df2]).reset_index(drop=True)
        expected = DataFrame(dict(time=[ts1, ts2]), dtype=object)
        assert_frame_equal(results, expected)

        results = pd.concat([df1, df3]).reset_index(drop=True)
        expected = DataFrame(dict(time=[ts1, ts3]), dtype=object)
        assert_frame_equal(results, expected)

        results = pd.concat([df2, df3]).reset_index(drop=True)
        expected = DataFrame(dict(time=[ts2, ts3]))
        assert_frame_equal(results, expected)

    def test_append_series_dict(self):
        df = DataFrame(np.random.randn(5, 4),
                       columns=['foo', 'bar', 'baz', 'qux'])

        series = df.ix[4]
        with assertRaisesRegexp(ValueError, 'Indexes have overlapping values'):
            df.append(series, verify_integrity=True)
        series.name = None
        with assertRaisesRegexp(TypeError, 'Can only append a Series if '
                                'ignore_index=True'):
            df.append(series, verify_integrity=True)

        result = df.append(series[::-1], ignore_index=True)
        expected = df.append(DataFrame({0: series[::-1]}, index=df.columns).T,
                             ignore_index=True)
        assert_frame_equal(result, expected)

        # dict
        result = df.append(series.to_dict(), ignore_index=True)
        assert_frame_equal(result, expected)

        result = df.append(series[::-1][:3], ignore_index=True)
        expected = df.append(DataFrame({0: series[::-1][:3]}).T,
                             ignore_index=True)
        assert_frame_equal(result, expected.ix[:, result.columns])

        # can append when name set
        row = df.ix[4]
        row.name = 5
        result = df.append(row)
        expected = df.append(df[-1:], ignore_index=True)
        assert_frame_equal(result, expected)

    def test_append_list_of_series_dicts(self):
        df = DataFrame(np.random.randn(5, 4),
                       columns=['foo', 'bar', 'baz', 'qux'])

        dicts = [x.to_dict() for idx, x in df.iterrows()]

        result = df.append(dicts, ignore_index=True)
        expected = df.append(df, ignore_index=True)
        assert_frame_equal(result, expected)

        # different columns
        dicts = [{'foo': 1, 'bar': 2, 'baz': 3, 'peekaboo': 4},
                 {'foo': 5, 'bar': 6, 'baz': 7, 'peekaboo': 8}]
        result = df.append(dicts, ignore_index=True)
        expected = df.append(DataFrame(dicts), ignore_index=True)
        assert_frame_equal(result, expected)

    def test_append_empty_dataframe(self):

        # Empty df append empty df
        df1 = DataFrame([])
        df2 = DataFrame([])
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        # Non-empty df append empty df
        df1 = DataFrame(np.random.randn(5, 2))
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        # Empty df with columns append empty df
        df1 = DataFrame(columns=['bar', 'foo'])
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        # Non-Empty df with columns append empty df
        df1 = DataFrame(np.random.randn(5, 2), columns=['bar', 'foo'])
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

    def test_append_dtypes(self):

        # GH 5754
        # row appends of different dtypes (so need to do by-item)
        # can sometimes infer the correct type

        df1 = DataFrame({'bar': Timestamp('20130101')}, index=lrange(5))
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        df1 = DataFrame({'bar': Timestamp('20130101')}, index=lrange(1))
        df2 = DataFrame({'bar': 'foo'}, index=lrange(1, 2))
        result = df1.append(df2)
        expected = DataFrame({'bar': [Timestamp('20130101'), 'foo']})
        assert_frame_equal(result, expected)

        df1 = DataFrame({'bar': Timestamp('20130101')}, index=lrange(1))
        df2 = DataFrame({'bar': np.nan}, index=lrange(1, 2))
        result = df1.append(df2)
        expected = DataFrame(
            {'bar': Series([Timestamp('20130101'), np.nan], dtype='M8[ns]')})
        assert_frame_equal(result, expected)

        df1 = DataFrame({'bar': Timestamp('20130101')}, index=lrange(1))
        df2 = DataFrame({'bar': np.nan}, index=lrange(1, 2), dtype=object)
        result = df1.append(df2)
        expected = DataFrame(
            {'bar': Series([Timestamp('20130101'), np.nan], dtype='M8[ns]')})
        assert_frame_equal(result, expected)

        df1 = DataFrame({'bar': np.nan}, index=lrange(1))
        df2 = DataFrame({'bar': Timestamp('20130101')}, index=lrange(1, 2))
        result = df1.append(df2)
        expected = DataFrame(
            {'bar': Series([np.nan, Timestamp('20130101')], dtype='M8[ns]')})
        assert_frame_equal(result, expected)

        df1 = DataFrame({'bar': Timestamp('20130101')}, index=lrange(1))
        df2 = DataFrame({'bar': 1}, index=lrange(1, 2), dtype=object)
        result = df1.append(df2)
        expected = DataFrame({'bar': Series([Timestamp('20130101'), 1])})
        assert_frame_equal(result, expected)

    def test_combine_first(self):
        # disjoint
        head, tail = self.frame[:5], self.frame[5:]

        combined = head.combine_first(tail)
        reordered_frame = self.frame.reindex(combined.index)
        assert_frame_equal(combined, reordered_frame)
        self.assertTrue(tm.equalContents(combined.columns, self.frame.columns))
        assert_series_equal(combined['A'], reordered_frame['A'])

        # same index
        fcopy = self.frame.copy()
        fcopy['A'] = 1
        del fcopy['C']

        fcopy2 = self.frame.copy()
        fcopy2['B'] = 0
        del fcopy2['D']

        combined = fcopy.combine_first(fcopy2)

        self.assertTrue((combined['A'] == 1).all())
        assert_series_equal(combined['B'], fcopy['B'])
        assert_series_equal(combined['C'], fcopy2['C'])
        assert_series_equal(combined['D'], fcopy['D'])

        # overlap
        head, tail = reordered_frame[:10].copy(), reordered_frame
        head['A'] = 1

        combined = head.combine_first(tail)
        self.assertTrue((combined['A'][:10] == 1).all())

        # reverse overlap
        tail['A'][:10] = 0
        combined = tail.combine_first(head)
        self.assertTrue((combined['A'][:10] == 0).all())

        # no overlap
        f = self.frame[:10]
        g = self.frame[10:]
        combined = f.combine_first(g)
        assert_series_equal(combined['A'].reindex(f.index), f['A'])
        assert_series_equal(combined['A'].reindex(g.index), g['A'])

        # corner cases
        comb = self.frame.combine_first(self.empty)
        assert_frame_equal(comb, self.frame)

        comb = self.empty.combine_first(self.frame)
        assert_frame_equal(comb, self.frame)

        comb = self.frame.combine_first(DataFrame(index=["faz", "boo"]))
        self.assertTrue("faz" in comb.index)

        # #2525
        df = DataFrame({'a': [1]}, index=[datetime(2012, 1, 1)])
        df2 = DataFrame({}, columns=['b'])
        result = df.combine_first(df2)
        self.assertTrue('b' in result)

    def test_combine_first_mixed_bug(self):
        idx = Index(['a', 'b', 'c', 'e'])
        ser1 = Series([5.0, -9.0, 4.0, 100.], index=idx)
        ser2 = Series(['a', 'b', 'c', 'e'], index=idx)
        ser3 = Series([12, 4, 5, 97], index=idx)

        frame1 = DataFrame({"col0": ser1,
                            "col2": ser2,
                            "col3": ser3})

        idx = Index(['a', 'b', 'c', 'f'])
        ser1 = Series([5.0, -9.0, 4.0, 100.], index=idx)
        ser2 = Series(['a', 'b', 'c', 'f'], index=idx)
        ser3 = Series([12, 4, 5, 97], index=idx)

        frame2 = DataFrame({"col1": ser1,
                            "col2": ser2,
                            "col5": ser3})

        combined = frame1.combine_first(frame2)
        self.assertEqual(len(combined.columns), 5)

        # gh 3016 (same as in update)
        df = DataFrame([[1., 2., False, True], [4., 5., True, False]],
                       columns=['A', 'B', 'bool1', 'bool2'])

        other = DataFrame([[45, 45]], index=[0], columns=['A', 'B'])
        result = df.combine_first(other)
        assert_frame_equal(result, df)

        df.ix[0, 'A'] = np.nan
        result = df.combine_first(other)
        df.ix[0, 'A'] = 45
        assert_frame_equal(result, df)

        # doc example
        df1 = DataFrame({'A': [1., np.nan, 3., 5., np.nan],
                         'B': [np.nan, 2., 3., np.nan, 6.]})

        df2 = DataFrame({'A': [5., 2., 4., np.nan, 3., 7.],
                         'B': [np.nan, np.nan, 3., 4., 6., 8.]})

        result = df1.combine_first(df2)
        expected = DataFrame(
            {'A': [1, 2, 3, 5, 3, 7.], 'B': [np.nan, 2, 3, 4, 6, 8]})
        assert_frame_equal(result, expected)

        # GH3552, return object dtype with bools
        df1 = DataFrame(
            [[np.nan, 3., True], [-4.6, np.nan, True], [np.nan, 7., False]])
        df2 = DataFrame(
            [[-42.6, np.nan, True], [-5., 1.6, False]], index=[1, 2])

        result = df1.combine_first(df2)[2]
        expected = Series([True, True, False], name=2)
        assert_series_equal(result, expected)

        # GH 3593, converting datetime64[ns] incorrecly
        df0 = DataFrame({"a": [datetime(2000, 1, 1),
                               datetime(2000, 1, 2),
                               datetime(2000, 1, 3)]})
        df1 = DataFrame({"a": [None, None, None]})
        df2 = df1.combine_first(df0)
        assert_frame_equal(df2, df0)

        df2 = df0.combine_first(df1)
        assert_frame_equal(df2, df0)

        df0 = DataFrame({"a": [datetime(2000, 1, 1),
                               datetime(2000, 1, 2),
                               datetime(2000, 1, 3)]})
        df1 = DataFrame({"a": [datetime(2000, 1, 2), None, None]})
        df2 = df1.combine_first(df0)
        result = df0.copy()
        result.iloc[0, :] = df1.iloc[0, :]
        assert_frame_equal(df2, result)

        df2 = df0.combine_first(df1)
        assert_frame_equal(df2, df0)

    def test_update(self):
        df = DataFrame([[1.5, nan, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[3.6, 2., np.nan],
                           [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other)

        expected = DataFrame([[1.5, nan, 3],
                              [3.6, 2, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 7.]])
        assert_frame_equal(df, expected)

    def test_update_dtypes(self):

        # gh 3016
        df = DataFrame([[1., 2., False, True], [4., 5., True, False]],
                       columns=['A', 'B', 'bool1', 'bool2'])

        other = DataFrame([[45, 45]], index=[0], columns=['A', 'B'])
        df.update(other)

        expected = DataFrame([[45., 45., False, True], [4., 5., True, False]],
                             columns=['A', 'B', 'bool1', 'bool2'])
        assert_frame_equal(df, expected)

    def test_update_nooverwrite(self):
        df = DataFrame([[1.5, nan, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[3.6, 2., np.nan],
                           [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other, overwrite=False)

        expected = DataFrame([[1.5, nan, 3],
                              [1.5, 2, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 3.]])
        assert_frame_equal(df, expected)

    def test_update_filtered(self):
        df = DataFrame([[1.5, nan, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[3.6, 2., np.nan],
                           [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other, filter_func=lambda x: x > 2)

        expected = DataFrame([[1.5, nan, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 7.]])
        assert_frame_equal(df, expected)

    def test_update_raise(self):
        df = DataFrame([[1.5, 1, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[2., nan],
                           [nan, 7]], index=[1, 3], columns=[1, 2])
        with assertRaisesRegexp(ValueError, "Data overlaps"):
            df.update(other, raise_conflict=True)

    def test_update_from_non_df(self):
        d = {'a': Series([1, 2, 3, 4]), 'b': Series([5, 6, 7, 8])}
        df = DataFrame(d)

        d['a'] = Series([5, 6, 7, 8])
        df.update(d)

        expected = DataFrame(d)

        assert_frame_equal(df, expected)

        d = {'a': [1, 2, 3, 4], 'b': [5, 6, 7, 8]}
        df = DataFrame(d)

        d['a'] = [5, 6, 7, 8]
        df.update(d)

        expected = DataFrame(d)

        assert_frame_equal(df, expected)

    def test_join_str_datetime(self):
        str_dates = ['20120209', '20120222']
        dt_dates = [datetime(2012, 2, 9), datetime(2012, 2, 22)]

        A = DataFrame(str_dates, index=lrange(2), columns=['aa'])
        C = DataFrame([[1, 2], [3, 4]], index=str_dates, columns=dt_dates)

        tst = A.join(C, on='aa')

        self.assertEqual(len(tst.columns), 3)

    def test_join_multiindex_leftright(self):
        # GH 10741
        df1 = (pd.DataFrame([['a', 'x', 0.471780], ['a', 'y', 0.774908],
                             ['a', 'z', 0.563634], ['b', 'x', -0.353756],
                             ['b', 'y', 0.368062], ['b', 'z', -1.721840],
                             ['c', 'x', 1], ['c', 'y', 2], ['c', 'z', 3]],
                            columns=['first', 'second', 'value1'])
               .set_index(['first', 'second']))

        df2 = (pd.DataFrame([['a', 10], ['b', 20]],
                            columns=['first', 'value2'])
               .set_index(['first']))

        exp = pd.DataFrame([[0.471780, 10], [0.774908, 10], [0.563634, 10],
                            [-0.353756, 20], [0.368062, 20],
                            [-1.721840, 20],
                            [1.000000, np.nan], [2.000000, np.nan],
                            [3.000000, np.nan]],
                           index=df1.index, columns=['value1', 'value2'])

        # these must be the same results (but columns are flipped)
        assert_frame_equal(df1.join(df2, how='left'), exp)
        assert_frame_equal(df2.join(df1, how='right'),
                           exp[['value2', 'value1']])

        exp_idx = pd.MultiIndex.from_product([['a', 'b'], ['x', 'y', 'z']],
                                             names=['first', 'second'])
        exp = pd.DataFrame([[0.471780, 10], [0.774908, 10], [0.563634, 10],
                            [-0.353756, 20], [0.368062, 20], [-1.721840, 20]],
                           index=exp_idx, columns=['value1', 'value2'])

        assert_frame_equal(df1.join(df2, how='right'), exp)
        assert_frame_equal(df2.join(df1, how='left'),
                           exp[['value2', 'value1']])
