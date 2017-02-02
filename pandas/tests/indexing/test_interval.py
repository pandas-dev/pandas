import pytest
import numpy as np
import pandas as pd

from pandas import Series, DataFrame, IntervalIndex, Interval
import pandas.util.testing as tm


class TestIntervalIndex(tm.TestCase):

    def setUp(self):
        self.s = Series(np.arange(5), IntervalIndex.from_breaks(np.arange(6)))

    def test_loc_getitem_series(self):

        s = self.s
        expected = 0
        self.assertEqual(expected, s.loc[0.5])
        self.assertEqual(expected, s.loc[1])
        self.assertEqual(expected, s.loc[Interval(0, 1)])
        self.assertRaises(KeyError, s.loc.__getitem__, 0)

        expected = s.iloc[:3]
        tm.assert_series_equal(expected, s.loc[:3])
        tm.assert_series_equal(expected, s.loc[:2.5])
        tm.assert_series_equal(expected, s.loc[0.1:2.5])
        tm.assert_series_equal(expected, s.loc[-1:3])

        expected = s.iloc[1:4]
        tm.assert_series_equal(expected, s.loc[[1.5, 2.5, 3.5]])
        tm.assert_series_equal(expected, s.loc[[2, 3, 4]])
        tm.assert_series_equal(expected, s.loc[[1.5, 3, 4]])

        expected = s.iloc[2:5]
        tm.assert_series_equal(expected, s.loc[s >= 2])

        expected = s.iloc[2:5]
        result = s.loc[[pd.Interval(3, 6)]]
        tm.assert_series_equal(expected, result)

        expected = s.iloc[2:4]
        result = s.loc[[pd.Interval(3, 5)]]
        tm.assert_series_equal(expected, result)

        expected = s.iloc[[2, 3, 4, 2, 3, 4]]
        result = s.loc[[pd.Interval(3, 6), pd.Interval(3, 6)]]
        tm.assert_series_equal(expected, result)

        # slice of interval
        with pytest.raises(NotImplementedError):
            result = s.loc[pd.Interval(3, 6):]

    def test_loc_non_matching(self):
        s = self.s

        # TODO: We are getting at least 1 matching
        # interval so this meets our current semantics
        expected = s.iloc[[2, 3, 4]]
        result = s.loc[[-1, 3, 4, 5]]
        tm.assert_series_equal(expected, result)

    def test_getitem_series(self):

        s = self.s
        expected = 0
        self.assertEqual(expected, s[0.5])
        self.assertEqual(expected, s[1])
        self.assertEqual(expected, s[Interval(0, 1)])
        self.assertRaises(KeyError, s.__getitem__, 0)

        expected = s.iloc[:3]
        tm.assert_series_equal(expected, s[:3])
        tm.assert_series_equal(expected, s[:2.5])
        tm.assert_series_equal(expected, s[0.1:2.5])
        tm.assert_series_equal(expected, s[-1:3])

        expected = s.iloc[1:4]
        tm.assert_series_equal(expected, s[[1.5, 2.5, 3.5]])
        tm.assert_series_equal(expected, s[[2, 3, 4]])
        tm.assert_series_equal(expected, s[[1.5, 3, 4]])

        expected = s.iloc[2:5]
        tm.assert_series_equal(expected, s[s >= 2])

        expected = s.iloc[2:5]
        result = s[[pd.Interval(3, 6)]]
        tm.assert_series_equal(expected, result)

        # slice of interval
        with pytest.raises(NotImplementedError):
            result = s[pd.Interval(3, 6):]

        # slice of scalar
        with pytest.raises(NotImplementedError):
            s[0:4:2]

    def test_large_series(self):
        s = Series(np.arange(1000000),
                   index=IntervalIndex.from_breaks(np.arange(1000001)))

        result1 = s.loc[:80000]
        result2 = s.loc[0:80000]
        result3 = s.loc[0:80000:1]
        tm.assert_series_equal(result1, result2)
        tm.assert_series_equal(result1, result3)

    def test_loc_getitem_frame(self):

        df = DataFrame({'A': range(10)})
        s = pd.cut(df.A, 5)
        df['B'] = s
        df = df.set_index('B')

        result = df.loc[4]
        expected = df.iloc[4:6]
        tm.assert_frame_equal(result, expected)

        def f():
            df.loc[10]

        self.assertRaises(KeyError, f)

        # single list-like
        result = df.loc[[4]]
        expected = df.iloc[4:6]
        tm.assert_frame_equal(result, expected)

        # non-unique
        result = df.loc[[4, 5]]
        expected = df.take([4, 5, 4, 5])
        tm.assert_frame_equal(result, expected)

        def f():
            df.loc[[10]]

        self.assertRaises(KeyError, f)

        # partial missing
        result = df.loc[[10, 4]]
        expected = df.iloc[4:6]
        tm.assert_frame_equal(result, expected)
