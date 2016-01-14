# -*- coding: utf-8 -*-

from __future__ import print_function

from pandas import DataFrame, Series, MultiIndex, Panel
import pandas as pd

from pandas.util.testing import (assert_frame_equal,
                                 SubclassedDataFrame)

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
        df = SubclassedDataFrame({'X': [1, 2, 3], 'Y': [1, 2, 3]},
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
        assert_frame_equal(df, unpickled)
        self.assertEqual(df._metadata, unpickled._metadata)
        self.assertEqual(df.testattr, unpickled.testattr)

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
