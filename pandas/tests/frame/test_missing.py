# -*- coding: utf-8 -*-

from __future__ import print_function

from numpy import nan, random
import numpy as np

from pandas.compat import lrange
from pandas import (DataFrame, Series, Timestamp,
                    date_range)
import pandas as pd

from pandas.util.testing import (assert_series_equal,
                                 assert_frame_equal,
                                 assertRaisesRegexp)

import pandas.util.testing as tm
from pandas.tests.frame.common import TestData, _check_mixed_float


class TestDataFrameMissingData(tm.TestCase, TestData):

    _multiprocess_can_split_ = True

    def test_dropEmptyRows(self):
        N = len(self.frame.index)
        mat = random.randn(N)
        mat[:5] = nan

        frame = DataFrame({'foo': mat}, index=self.frame.index)
        original = Series(mat, index=self.frame.index, name='foo')
        expected = original.dropna()
        inplace_frame1, inplace_frame2 = frame.copy(), frame.copy()

        smaller_frame = frame.dropna(how='all')
        # check that original was preserved
        assert_series_equal(frame['foo'], original)
        inplace_frame1.dropna(how='all', inplace=True)
        assert_series_equal(smaller_frame['foo'], expected)
        assert_series_equal(inplace_frame1['foo'], expected)

        smaller_frame = frame.dropna(how='all', subset=['foo'])
        inplace_frame2.dropna(how='all', subset=['foo'], inplace=True)
        assert_series_equal(smaller_frame['foo'], expected)
        assert_series_equal(inplace_frame2['foo'], expected)

    def test_dropIncompleteRows(self):
        N = len(self.frame.index)
        mat = random.randn(N)
        mat[:5] = nan

        frame = DataFrame({'foo': mat}, index=self.frame.index)
        frame['bar'] = 5
        original = Series(mat, index=self.frame.index, name='foo')
        inp_frame1, inp_frame2 = frame.copy(), frame.copy()

        smaller_frame = frame.dropna()
        assert_series_equal(frame['foo'], original)
        inp_frame1.dropna(inplace=True)
        self.assert_numpy_array_equal(smaller_frame['foo'], mat[5:])
        self.assert_numpy_array_equal(inp_frame1['foo'], mat[5:])

        samesize_frame = frame.dropna(subset=['bar'])
        assert_series_equal(frame['foo'], original)
        self.assertTrue((frame['bar'] == 5).all())
        inp_frame2.dropna(subset=['bar'], inplace=True)
        self.assertTrue(samesize_frame.index.equals(self.frame.index))
        self.assertTrue(inp_frame2.index.equals(self.frame.index))

    def test_dropna(self):
        df = DataFrame(np.random.randn(6, 4))
        df[2][:2] = nan

        dropped = df.dropna(axis=1)
        expected = df.ix[:, [0, 1, 3]]
        inp = df.copy()
        inp.dropna(axis=1, inplace=True)
        assert_frame_equal(dropped, expected)
        assert_frame_equal(inp, expected)

        dropped = df.dropna(axis=0)
        expected = df.ix[lrange(2, 6)]
        inp = df.copy()
        inp.dropna(axis=0, inplace=True)
        assert_frame_equal(dropped, expected)
        assert_frame_equal(inp, expected)

        # threshold
        dropped = df.dropna(axis=1, thresh=5)
        expected = df.ix[:, [0, 1, 3]]
        inp = df.copy()
        inp.dropna(axis=1, thresh=5, inplace=True)
        assert_frame_equal(dropped, expected)
        assert_frame_equal(inp, expected)

        dropped = df.dropna(axis=0, thresh=4)
        expected = df.ix[lrange(2, 6)]
        inp = df.copy()
        inp.dropna(axis=0, thresh=4, inplace=True)
        assert_frame_equal(dropped, expected)
        assert_frame_equal(inp, expected)

        dropped = df.dropna(axis=1, thresh=4)
        assert_frame_equal(dropped, df)

        dropped = df.dropna(axis=1, thresh=3)
        assert_frame_equal(dropped, df)

        # subset
        dropped = df.dropna(axis=0, subset=[0, 1, 3])
        inp = df.copy()
        inp.dropna(axis=0, subset=[0, 1, 3], inplace=True)
        assert_frame_equal(dropped, df)
        assert_frame_equal(inp, df)

        # all
        dropped = df.dropna(axis=1, how='all')
        assert_frame_equal(dropped, df)

        df[2] = nan
        dropped = df.dropna(axis=1, how='all')
        expected = df.ix[:, [0, 1, 3]]
        assert_frame_equal(dropped, expected)

        # bad input
        self.assertRaises(ValueError, df.dropna, axis=3)

    def test_drop_and_dropna_caching(self):
        # tst that cacher updates
        original = Series([1, 2, np.nan], name='A')
        expected = Series([1, 2], dtype=original.dtype, name='A')
        df = pd.DataFrame({'A': original.values.copy()})
        df2 = df.copy()
        df['A'].dropna()
        assert_series_equal(df['A'], original)
        df['A'].dropna(inplace=True)
        assert_series_equal(df['A'], expected)
        df2['A'].drop([1])
        assert_series_equal(df2['A'], original)
        df2['A'].drop([1], inplace=True)
        assert_series_equal(df2['A'], original.drop([1]))

    def test_dropna_corner(self):
        # bad input
        self.assertRaises(ValueError, self.frame.dropna, how='foo')
        self.assertRaises(TypeError, self.frame.dropna, how=None)
        # non-existent column - 8303
        self.assertRaises(KeyError, self.frame.dropna, subset=['A', 'X'])

    def test_dropna_multiple_axes(self):
        df = DataFrame([[1, np.nan, 2, 3],
                        [4, np.nan, 5, 6],
                        [np.nan, np.nan, np.nan, np.nan],
                        [7, np.nan, 8, 9]])
        cp = df.copy()
        result = df.dropna(how='all', axis=[0, 1])
        result2 = df.dropna(how='all', axis=(0, 1))
        expected = df.dropna(how='all').dropna(how='all', axis=1)

        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)
        assert_frame_equal(df, cp)

        inp = df.copy()
        inp.dropna(how='all', axis=(0, 1), inplace=True)
        assert_frame_equal(inp, expected)

    def test_fillna(self):
        self.tsframe.ix[:5, 'A'] = nan
        self.tsframe.ix[-5:, 'A'] = nan

        zero_filled = self.tsframe.fillna(0)
        self.assertTrue((zero_filled.ix[:5, 'A'] == 0).all())

        padded = self.tsframe.fillna(method='pad')
        self.assertTrue(np.isnan(padded.ix[:5, 'A']).all())
        self.assertTrue((padded.ix[-5:, 'A'] == padded.ix[-5, 'A']).all())

        # mixed type
        self.mixed_frame.ix[5:20, 'foo'] = nan
        self.mixed_frame.ix[-10:, 'A'] = nan
        result = self.mixed_frame.fillna(value=0)
        result = self.mixed_frame.fillna(method='pad')

        self.assertRaises(ValueError, self.tsframe.fillna)
        self.assertRaises(ValueError, self.tsframe.fillna, 5, method='ffill')

        # mixed numeric (but no float16)
        mf = self.mixed_float.reindex(columns=['A', 'B', 'D'])
        mf.ix[-10:, 'A'] = nan
        result = mf.fillna(value=0)
        _check_mixed_float(result, dtype=dict(C=None))

        result = mf.fillna(method='pad')
        _check_mixed_float(result, dtype=dict(C=None))

        # empty frame (GH #2778)
        df = DataFrame(columns=['x'])
        for m in ['pad', 'backfill']:
            df.x.fillna(method=m, inplace=1)
            df.x.fillna(method=m)

        # with different dtype (GH3386)
        df = DataFrame([['a', 'a', np.nan, 'a'], [
                       'b', 'b', np.nan, 'b'], ['c', 'c', np.nan, 'c']])

        result = df.fillna({2: 'foo'})
        expected = DataFrame([['a', 'a', 'foo', 'a'],
                              ['b', 'b', 'foo', 'b'],
                              ['c', 'c', 'foo', 'c']])
        assert_frame_equal(result, expected)

        df.fillna({2: 'foo'}, inplace=True)
        assert_frame_equal(df, expected)

        # limit and value
        df = DataFrame(np.random.randn(10, 3))
        df.iloc[2:7, 0] = np.nan
        df.iloc[3:5, 2] = np.nan

        expected = df.copy()
        expected.iloc[2, 0] = 999
        expected.iloc[3, 2] = 999
        result = df.fillna(999, limit=1)
        assert_frame_equal(result, expected)

        # with datelike
        # GH 6344
        df = DataFrame({
            'Date': [pd.NaT, Timestamp("2014-1-1")],
            'Date2': [Timestamp("2013-1-1"), pd.NaT]
        })

        expected = df.copy()
        expected['Date'] = expected['Date'].fillna(df.ix[0, 'Date2'])
        result = df.fillna(value={'Date': df['Date2']})
        assert_frame_equal(result, expected)

    def test_fillna_dtype_conversion(self):
        # make sure that fillna on an empty frame works
        df = DataFrame(index=["A", "B", "C"], columns=[1, 2, 3, 4, 5])
        result = df.get_dtype_counts().sort_values()
        expected = Series({'object': 5})
        assert_series_equal(result, expected)

        result = df.fillna(1)
        expected = DataFrame(1, index=["A", "B", "C"], columns=[1, 2, 3, 4, 5])
        result = result.get_dtype_counts().sort_values()
        expected = Series({'int64': 5})
        assert_series_equal(result, expected)

        # empty block
        df = DataFrame(index=lrange(3), columns=['A', 'B'], dtype='float64')
        result = df.fillna('nan')
        expected = DataFrame('nan', index=lrange(3), columns=['A', 'B'])
        assert_frame_equal(result, expected)

        # equiv of replace
        df = DataFrame(dict(A=[1, np.nan], B=[1., 2.]))
        for v in ['', 1, np.nan, 1.0]:
            expected = df.replace(np.nan, v)
            result = df.fillna(v)
            assert_frame_equal(result, expected)

    def test_fillna_datetime_columns(self):
        # GH 7095
        df = pd.DataFrame({'A': [-1, -2, np.nan],
                           'B': date_range('20130101', periods=3),
                           'C': ['foo', 'bar', None],
                           'D': ['foo2', 'bar2', None]},
                          index=date_range('20130110', periods=3))
        result = df.fillna('?')
        expected = pd.DataFrame({'A': [-1, -2, '?'],
                                 'B': date_range('20130101', periods=3),
                                 'C': ['foo', 'bar', '?'],
                                 'D': ['foo2', 'bar2', '?']},
                                index=date_range('20130110', periods=3))
        self.assert_frame_equal(result, expected)

        df = pd.DataFrame({'A': [-1, -2, np.nan],
                           'B': [pd.Timestamp('2013-01-01'),
                                 pd.Timestamp('2013-01-02'), pd.NaT],
                           'C': ['foo', 'bar', None],
                           'D': ['foo2', 'bar2', None]},
                          index=date_range('20130110', periods=3))
        result = df.fillna('?')
        expected = pd.DataFrame({'A': [-1, -2, '?'],
                                 'B': [pd.Timestamp('2013-01-01'),
                                       pd.Timestamp('2013-01-02'), '?'],
                                 'C': ['foo', 'bar', '?'],
                                 'D': ['foo2', 'bar2', '?']},
                                index=pd.date_range('20130110', periods=3))
        self.assert_frame_equal(result, expected)

    def test_ffill(self):
        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan

        assert_frame_equal(self.tsframe.ffill(),
                           self.tsframe.fillna(method='ffill'))

    def test_bfill(self):
        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan

        assert_frame_equal(self.tsframe.bfill(),
                           self.tsframe.fillna(method='bfill'))

    def test_fillna_skip_certain_blocks(self):
        # don't try to fill boolean, int blocks

        df = DataFrame(np.random.randn(10, 4).astype(int))

        # it works!
        df.fillna(np.nan)

    def test_fillna_inplace(self):
        df = DataFrame(np.random.randn(10, 4))
        df[1][:4] = np.nan
        df[3][-4:] = np.nan

        expected = df.fillna(value=0)
        self.assertIsNot(expected, df)

        df.fillna(value=0, inplace=True)
        assert_frame_equal(df, expected)

        df[1][:4] = np.nan
        df[3][-4:] = np.nan
        expected = df.fillna(method='ffill')
        self.assertIsNot(expected, df)

        df.fillna(method='ffill', inplace=True)
        assert_frame_equal(df, expected)

    def test_fillna_dict_series(self):
        df = DataFrame({'a': [nan, 1, 2, nan, nan],
                        'b': [1, 2, 3, nan, nan],
                        'c': [nan, 1, 2, 3, 4]})

        result = df.fillna({'a': 0, 'b': 5})

        expected = df.copy()
        expected['a'] = expected['a'].fillna(0)
        expected['b'] = expected['b'].fillna(5)
        assert_frame_equal(result, expected)

        # it works
        result = df.fillna({'a': 0, 'b': 5, 'd': 7})

        # Series treated same as dict
        result = df.fillna(df.max())
        expected = df.fillna(df.max().to_dict())
        assert_frame_equal(result, expected)

        # disable this for now
        with assertRaisesRegexp(NotImplementedError, 'column by column'):
            df.fillna(df.max(1), axis=1)

    def test_fillna_dataframe(self):
        # GH 8377
        df = DataFrame({'a': [nan, 1, 2, nan, nan],
                        'b': [1, 2, 3, nan, nan],
                        'c': [nan, 1, 2, 3, 4]},
                       index=list('VWXYZ'))

        # df2 may have different index and columns
        df2 = DataFrame({'a': [nan, 10, 20, 30, 40],
                         'b': [50, 60, 70, 80, 90],
                         'foo': ['bar'] * 5},
                        index=list('VWXuZ'))

        result = df.fillna(df2)

        # only those columns and indices which are shared get filled
        expected = DataFrame({'a': [nan, 1, 2, nan, 40],
                              'b': [1, 2, 3, nan, 90],
                              'c': [nan, 1, 2, 3, 4]},
                             index=list('VWXYZ'))

        assert_frame_equal(result, expected)

    def test_fillna_columns(self):
        df = DataFrame(np.random.randn(10, 10))
        df.values[:, ::2] = np.nan

        result = df.fillna(method='ffill', axis=1)
        expected = df.T.fillna(method='pad').T
        assert_frame_equal(result, expected)

        df.insert(6, 'foo', 5)
        result = df.fillna(method='ffill', axis=1)
        expected = df.astype(float).fillna(method='ffill', axis=1)
        assert_frame_equal(result, expected)

    def test_fillna_invalid_method(self):
        with assertRaisesRegexp(ValueError, 'ffil'):
            self.frame.fillna(method='ffil')

    def test_fillna_invalid_value(self):
        # list
        self.assertRaises(TypeError, self.frame.fillna, [1, 2])
        # tuple
        self.assertRaises(TypeError, self.frame.fillna, (1, 2))
        # frame with series
        self.assertRaises(ValueError, self.frame.iloc[:, 0].fillna,
                          self.frame)

    def test_fillna_col_reordering(self):
        cols = ["COL." + str(i) for i in range(5, 0, -1)]
        data = np.random.rand(20, 5)
        df = DataFrame(index=lrange(20), columns=cols, data=data)
        filled = df.fillna(method='ffill')
        self.assertEqual(df.columns.tolist(), filled.columns.tolist())

    def test_fill_corner(self):
        self.mixed_frame.ix[5:20, 'foo'] = nan
        self.mixed_frame.ix[-10:, 'A'] = nan

        filled = self.mixed_frame.fillna(value=0)
        self.assertTrue((filled.ix[5:20, 'foo'] == 0).all())
        del self.mixed_frame['foo']

        empty_float = self.frame.reindex(columns=[])

        # TODO(wesm): unused?
        result = empty_float.fillna(value=0)  # noqa
