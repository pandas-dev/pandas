# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np

from pandas import DataFrame, Series, MultiIndex, Panel
import pandas as pd
import pandas.util.testing as tm

from pandas.tests.frame.common import TestData


class TestDataFrameSubclassing(tm.TestCase, TestData):

    _multiprocess_can_split_ = True

    def test_frame_subclassing_and_slicing(self):
        # Subclass frame and ensure it returns the right class on slicing it
        # In reference to PR 9632

        class CustomSeries(Series):

            @property
            def _constructor(self):
                return CustomSeries

            def custom_series_function(self):
                return 'OK'

        class CustomDataFrame(DataFrame):
            """
            Subclasses pandas DF, fills DF with simulation results, adds some
            custom plotting functions.
            """

            def __init__(self, *args, **kw):
                super(CustomDataFrame, self).__init__(*args, **kw)

            @property
            def _constructor(self):
                return CustomDataFrame

            _constructor_sliced = CustomSeries

            def custom_frame_function(self):
                return 'OK'

        data = {'col1': range(10),
                'col2': range(10)}
        cdf = CustomDataFrame(data)

        # Did we get back our own DF class?
        self.assertTrue(isinstance(cdf, CustomDataFrame))

        # Do we get back our own Series class after selecting a column?
        cdf_series = cdf.col1
        self.assertTrue(isinstance(cdf_series, CustomSeries))
        self.assertEqual(cdf_series.custom_series_function(), 'OK')

        # Do we get back our own DF class after slicing row-wise?
        cdf_rows = cdf[1:5]
        self.assertTrue(isinstance(cdf_rows, CustomDataFrame))
        self.assertEqual(cdf_rows.custom_frame_function(), 'OK')

        # Make sure sliced part of multi-index frame is custom class
        mcol = pd.MultiIndex.from_tuples([('A', 'A'), ('A', 'B')])
        cdf_multi = CustomDataFrame([[0, 1], [2, 3]], columns=mcol)
        self.assertTrue(isinstance(cdf_multi['A'], CustomDataFrame))

        mcol = pd.MultiIndex.from_tuples([('A', ''), ('B', '')])
        cdf_multi2 = CustomDataFrame([[0, 1], [2, 3]], columns=mcol)
        self.assertTrue(isinstance(cdf_multi2['A'], CustomSeries))

    def test_dataframe_metadata(self):
        df = tm.SubclassedDataFrame({'X': [1, 2, 3], 'Y': [1, 2, 3]},
                                    index=['a', 'b', 'c'])
        df.testattr = 'XXX'

        self.assertEqual(df.testattr, 'XXX')
        self.assertEqual(df[['X']].testattr, 'XXX')
        self.assertEqual(df.loc[['a', 'b'], :].testattr, 'XXX')
        self.assertEqual(df.iloc[[0, 1], :].testattr, 'XXX')

        # GH9776
        self.assertEqual(df.iloc[0:1, :].testattr, 'XXX')

        # GH10553
        unpickled = self.round_trip_pickle(df)
        tm.assert_frame_equal(df, unpickled)
        self.assertEqual(df._metadata, unpickled._metadata)
        self.assertEqual(df.testattr, unpickled.testattr)

    def test_indexing_sliced(self):
        # GH 11559
        df = tm.SubclassedDataFrame({'X': [1, 2, 3],
                                     'Y': [4, 5, 6],
                                     'Z': [7, 8, 9]},
                                    index=['a', 'b', 'c'])
        res = df.loc[:, 'X']
        exp = tm.SubclassedSeries([1, 2, 3], index=list('abc'), name='X')
        tm.assert_series_equal(res, exp)
        tm.assertIsInstance(res, tm.SubclassedSeries)

        res = df.iloc[:, 1]
        exp = tm.SubclassedSeries([4, 5, 6], index=list('abc'), name='Y')
        tm.assert_series_equal(res, exp)
        tm.assertIsInstance(res, tm.SubclassedSeries)

        res = df.ix[:, 'Z']
        exp = tm.SubclassedSeries([7, 8, 9], index=list('abc'), name='Z')
        tm.assert_series_equal(res, exp)
        tm.assertIsInstance(res, tm.SubclassedSeries)

        res = df.loc['a', :]
        exp = tm.SubclassedSeries([1, 4, 7], index=list('XYZ'), name='a')
        tm.assert_series_equal(res, exp)
        tm.assertIsInstance(res, tm.SubclassedSeries)

        res = df.iloc[1, :]
        exp = tm.SubclassedSeries([2, 5, 8], index=list('XYZ'), name='b')
        tm.assert_series_equal(res, exp)
        tm.assertIsInstance(res, tm.SubclassedSeries)

        res = df.ix['c', :]
        exp = tm.SubclassedSeries([3, 6, 9], index=list('XYZ'), name='c')
        tm.assert_series_equal(res, exp)
        tm.assertIsInstance(res, tm.SubclassedSeries)

    def test_to_panel_expanddim(self):
        # GH 9762

        class SubclassedFrame(DataFrame):

            @property
            def _constructor_expanddim(self):
                return SubclassedPanel

        class SubclassedPanel(Panel):
            pass

        index = MultiIndex.from_tuples([(0, 0), (0, 1), (0, 2)])
        df = SubclassedFrame({'X': [1, 2, 3], 'Y': [4, 5, 6]}, index=index)
        result = df.to_panel()
        self.assertTrue(isinstance(result, SubclassedPanel))
        expected = SubclassedPanel([[[1, 2, 3]], [[4, 5, 6]]],
                                   items=['X', 'Y'], major_axis=[0],
                                   minor_axis=[0, 1, 2],
                                   dtype='int64')
        tm.assert_panel_equal(result, expected)

    def test_subclass_attr_err_propagation(self):
        # GH 11808
        class A(DataFrame):

            @property
            def bar(self):
                return self.i_dont_exist
        with tm.assertRaisesRegexp(AttributeError, '.*i_dont_exist.*'):
            A().bar

    def test_subclass_align(self):
        # GH 12983
        df1 = tm.SubclassedDataFrame({'a': [1, 3, 5],
                                      'b': [1, 3, 5]}, index=list('ACE'))
        df2 = tm.SubclassedDataFrame({'c': [1, 2, 4],
                                      'd': [1, 2, 4]}, index=list('ABD'))

        res1, res2 = df1.align(df2, axis=0)
        exp1 = tm.SubclassedDataFrame({'a': [1, np.nan, 3, np.nan, 5],
                                       'b': [1, np.nan, 3, np.nan, 5]},
                                      index=list('ABCDE'))
        exp2 = tm.SubclassedDataFrame({'c': [1, 2, np.nan, 4, np.nan],
                                       'd': [1, 2, np.nan, 4, np.nan]},
                                      index=list('ABCDE'))
        tm.assertIsInstance(res1, tm.SubclassedDataFrame)
        tm.assert_frame_equal(res1, exp1)
        tm.assertIsInstance(res2, tm.SubclassedDataFrame)
        tm.assert_frame_equal(res2, exp2)

        res1, res2 = df1.a.align(df2.c)
        tm.assertIsInstance(res1, tm.SubclassedSeries)
        tm.assert_series_equal(res1, exp1.a)
        tm.assertIsInstance(res2, tm.SubclassedSeries)
        tm.assert_series_equal(res2, exp2.c)

    def test_subclass_align_combinations(self):
        # GH 12983
        df = tm.SubclassedDataFrame({'a': [1, 3, 5],
                                     'b': [1, 3, 5]}, index=list('ACE'))
        s = tm.SubclassedSeries([1, 2, 4], index=list('ABD'), name='x')

        # frame + series
        res1, res2 = df.align(s, axis=0)
        exp1 = pd.DataFrame({'a': [1, np.nan, 3, np.nan, 5],
                             'b': [1, np.nan, 3, np.nan, 5]},
                            index=list('ABCDE'))
        # name is lost when
        exp2 = pd.Series([1, 2, np.nan, 4, np.nan],
                         index=list('ABCDE'), name='x')

        tm.assertIsInstance(res1, tm.SubclassedDataFrame)
        tm.assert_frame_equal(res1, exp1)
        tm.assertIsInstance(res2, tm.SubclassedSeries)
        tm.assert_series_equal(res2, exp2)

        # series + frame
        res1, res2 = s.align(df)
        tm.assertIsInstance(res1, tm.SubclassedSeries)
        tm.assert_series_equal(res1, exp2)
        tm.assertIsInstance(res2, tm.SubclassedDataFrame)
        tm.assert_frame_equal(res2, exp1)

    def test_subclass_iterrows(self):
        # GH 13977
        df = tm.SubclassedDataFrame({'a': [1]})
        for i, row in df.iterrows():
            tm.assertIsInstance(row, tm.SubclassedSeries)
            tm.assert_series_equal(row, df.loc[i])

    def test_subclass_sparse_slice(self):
        rows = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
        ssdf = tm.SubclassedSparseDataFrame(rows)
        ssdf.testattr = "testattr"

        tm.assert_sp_frame_equal(ssdf.loc[:2],
                                 tm.SubclassedSparseDataFrame(rows[:3]))
        tm.assert_sp_frame_equal(ssdf.iloc[:2],
                                 tm.SubclassedSparseDataFrame(rows[:2]))
        tm.assert_sp_frame_equal(ssdf[:2],
                                 tm.SubclassedSparseDataFrame(rows[:2]))
        tm.assert_equal(ssdf.loc[:2].testattr, "testattr")
        tm.assert_equal(ssdf.iloc[:2].testattr, "testattr")
        tm.assert_equal(ssdf[:2].testattr, "testattr")

        tm.assert_sp_series_equal(ssdf.loc[1],
                                  tm.SubclassedSparseSeries(rows[1]),
                                  check_names=False)
        tm.assert_sp_series_equal(ssdf.iloc[1],
                                  tm.SubclassedSparseSeries(rows[1]),
                                  check_names=False)

    def test_subclass_sparse_transpose(self):
        ossdf = tm.SubclassedSparseDataFrame([[1, 2, 3],
                                              [4, 5, 6]])
        essdf = tm.SubclassedSparseDataFrame([[1, 4],
                                              [2, 5],
                                              [3, 6]])
        tm.assert_sp_frame_equal(ossdf.T, essdf)
